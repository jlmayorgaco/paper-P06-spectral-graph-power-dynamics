include(joinpath(@__DIR__,"..","nonlinear_codesign_20261001","ReducedDAE.jl"))
using .ReducedDAE,LinearAlgebra,CSV,DataFrames,TOML,Sockets
const R=ReducedDAE
const OUT=joinpath(R.ROOT,"reports","retuning_boundary_20261002","search")
mkpath(OUT);BLAS.set_num_threads(1)
const CTX=R.N.design_context(R.ROOT)
const COUNT=Ref(0)
say(x...)=(println(x...);flush(stdout))
unpack(p)=(p[1:10],exp.(p[11:20]),exp.(p[21:30]))
function evaluate(p,grad,horizon,bus,delta)
    COUNT[]+=1;tag=lpad(COUNT[],4,'0')
    rho,kp,ki=unpack(p)
    m=R.model(CTX,rho,kp,ki;bus,delta,dc_convention=:physical_supply)
    mo=R.N.simple_mode_sensitivities(CTX,rho,kp,ki)
    modal=(real(mo.lambda)+.05)/.05
    if grad
        tr=R.tangent_simulate(m;horizon,dt=.0125,monitor_buses=collect(1:39))
        tr.limiter_fraction>0 || error("Smooth sensitivity crossed an actuator limit")
        pk=R.tangent_peaks(tr)
        g=vcat(modal,[pk[1].peak/.5-1,pk[2].peak/.5-1,(.9-tr.vmin)/.1,
                      (tr.vmax-1.1)/.1,(.002-tr.limiter_fraction)/.1])
        J=vcat(real.(vcat(mo.rho,mo.Kp.*kp,mo.Ki.*ki))'/.05,
            pk[1].gradient'/.5,pk[2].gradient'/.5,-tr.gvmin'/.1,
            tr.gvmax'/.1,-tr.limiter_gradient'/.1)
        data=Dict("parameters"=>p,"constraints"=>g,"horizon"=>horizon,"bus"=>bus,"delta"=>delta,
            "F"=>pk[1].peak,"R"=>pk[2].peak,"limiter_fraction"=>tr.limiter_fraction,
            "Ftime"=>pk[1].time,"Rtime"=>pk[2].time,"limiter_time"=>tr.limiter_time)
        open(io->TOML.print(io,data),joinpath(OUT,tag*"_gradient.toml"),"w")
        CSV.write(joinpath(OUT,tag*"_jacobian.csv"),DataFrame(J,:auto))
        say("GRAD ",tag," maxg=",maximum(g)," horizon=",horizon)
        return vcat(g,vec(J))
    end
    run=R.simulate(m;horizon,dt=.01,tol=2e-10,wall_limit=180.,maxiters=700000,rotating_frame=true)
    run.ok || error("Incomplete trajectory $(run.retcode), t=$(run.last_time)")
    met=R.metrics(m,run;dt=.01,monitor_buses=collect(1:39),
        savepath=joinpath(OUT,tag*"_trajectory.csv"))
    g=vcat(modal,[met.Fpeak_Hz/.5-1,met.Rpeak_Hz_s/.5-1,(.9-met.Vmin)/.1,
                  (met.Vmax-1.1)/.1,(.002-met.limiter_fraction)/.1])
    data=Dict{String,Any}(string(k)=>v for (k,v) in pairs(met))
    merge!(data,Dict("parameters"=>p,"constraints"=>g,"bus"=>bus,"delta"=>delta,
        "alpha"=>real(mo.lambda),"horizon"=>horizon,"tag"=>tag))
    open(io->TOML.print(io,data),joinpath(OUT,tag*"_evaluation.toml"),"w")
    say("EVAL ",tag," maxg=",maximum(g)," F=",met.Fpeak_Hz," alpha=",real(mo.lambda))
    g
end
server=listen(ip"127.0.0.1",0)
say("READY ",getsockname(server)[2])
sock=accept(server)
try
    while isopen(sock)
        line=readline(sock);line=="QUIT" && break
        try
            cmd,hh,bb,dd,vals=split(line,';')
            p=parse.(Float64,split(vals,','))
            g=evaluate(p,cmd=="GRAD",parse(Float64,hh),parse(Int,bb),parse(Float64,dd))
            println(sock,join(g,','));flush(sock)
        catch err
            say("ORACLE_ERROR ",sprint(showerror,err))
            println(sock,"ERR "*replace(sprint(showerror,err),'\n'=>' '));flush(sock)
        end
    end
finally
    close(sock);close(server)
end
