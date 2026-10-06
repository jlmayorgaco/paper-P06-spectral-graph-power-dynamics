include("validate_pd.jl")

# A complete short trajectory before the first observed limiter crossing.
# This is a diagnostic prefix, never a substitute for the 60 s validation.
function diagnose()
    designlabel=isempty(ARGS) ? "improved" : ARGS[1]
    designpath=joinpath(OUT,"candidate_"*designlabel*".toml")
    d=TOML.parsefile(designpath);designsha=bytes2hex(sha256(read(designpath)))
    fileprefix="bus8_prefix"*(designlabel=="improved" ? "" : "_"*designlabel)
    rho=Float64.(d["rho"]);kp=Float64.(d["Kp"]);ki=Float64.(d["Ki"])
    base=PDRef.frozen_baseline();nw=PDRef.build_architecture(base,rho,kp,ki)
    state=PDRef.trim_state(nw,base,rho,kp,ki);rows=NamedTuple[]
    for (label,tol,dt) in (("coarse",1e-8,.005),("fine",1e-10,.0025))
      active=event_network(nw,8,100.0)
      elapsed=@elapsed sol=SciMLBase.solve(SciMLBase.ODEProblem(active,state,(0.0,6.0)),Rodas5P();
        callback=get_callbacks(active),initializealg=SciMLBase.NoInit(),saveat=dt,
        abstol=tol,reltol=tol,maxiters=250000)
      SciMLBase.successful_retcode(sol.retcode) && sol.t[end]>=6.0-1e-8 || error("short run incomplete")
      times=collect(0.0:dt:6.0);nt=length(times);phase=zeros(10,nt);voltage=zeros(39,nt)
      prev=zeros(10);theta=zeros(10);initial=zeros(10)
      for (j,t) in enumerate(times)
        ss=NWState(sol,abs(t-1.0)<1e-12 ? t+1e-9 : t)
        for bus in 1:39
          ur=Float64(ss[VIndex(bus,:busbar₊u_r)]);ui=Float64(ss[VIndex(bus,:busbar₊u_i)])
          voltage[bus,j]=hypot(ur,ui)
          if bus>=30
            k=bus-29;raw=atan(ui,ur)
            if j==1;prev[k]=raw;theta[k]=raw;initial[k]=raw
            else;theta[k]+=mod(raw-prev[k]+pi,2pi)-pi;prev[k]=raw;end
            phase[k,j]=theta[k]-initial[k]
          end
        end
      end
      lag=Int(round(.5/dt));f=zeros(10,nt);r=zeros(10,nt)
      for j in 1:nt
        p1=j>lag ? phase[:,j-lag] : zeros(10)
        p2=j>2lag ? phase[:,j-2lag] : zeros(10)
        f[:,j]=(phase[:,j]-p1)/pi;r[:,j]=(phase[:,j]-2p1+p2)/(pi/2)
      end
      F=vec(maximum(abs.(f),dims=1));R=vec(maximum(abs.(r),dims=1))
      row=(;label,tol,dt,horizon_after_event_s=5.0,Fpeak_Hz=maximum(F),Rpeak_Hz_s=maximum(R),
        F_last_Hz=F[end],Vmin=minimum(voltage),Vmax=maximum(voltage),
        runtime_s=elapsed,retcode=string(sol.retcode),full_validation=false,candidate_sha256=designsha)
      push!(rows,row)
      data=DataFrame(time_after_event_s=times.-1,Fmax_Hz=F,Rmax_Hz_s=R,
        Vmin_pu=vec(minimum(voltage,dims=1)),Vmax_pu=vec(maximum(voltage,dims=1)),Vbus8_pu=voltage[8,:])
      CSV.write(joinpath(OUT,fileprefix*"_"*label*".csv"),data)
      println("BUS8_PREFIX ",row);flush(stdout)
    end
    CSV.write(joinpath(OUT,fileprefix*"_metrics.csv"),DataFrame(rows))
end
diagnose()
