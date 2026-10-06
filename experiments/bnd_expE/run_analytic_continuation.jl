using CSV, DataFrames, LinearAlgebra, TOML
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_design_e","PhysicalData.jl"))
include(joinpath(ROOT,"src","bnd_design_e","CollectiveModel.jl"))
include(joinpath(ROOT,"src","bnd_design_e","ClosureSpectrum.jl"))
include(joinpath(ROOT,"src","bnd_design_e","AnalyticContinuation.jl"))

function main()
domain=PhysicalData.freeze_design_domain(ROOT)
config=TOML.parsefile(domain.path)
net=CollectiveModel.frozen_network(ROOT)
buses=sort(collect(keys(net.sg)));n=length(buses)
nom=PhysicalData.nominal_pll_gains();Kp=fill(nom.Kp,n);Ki=fill(nom.Ki,n)
p=Float64.(domain.generators.SG_dispatch_initial_MW)
M=Matrix{Float64}(CSV.read(joinpath(ROOT,"reports","experiment_C","matrices","C0_M.csv"),DataFrame)[:,2:end])
Lg=Matrix{Float64}(CSV.read(joinpath(ROOT,"reports","experiment_C","matrices","C0_LG_primary.csv"),DataFrame)[:,2:end])
Qgraph=Symmetric(transpose(Lg)*(M\Lg))
Rinv=Diagonal(vcat(fill(nom.Kp^2,n),fill(nom.Ki^2,n)))
sigma=Float64(config["sigma_req_per_s"])
boundary=CSV.read(joinpath(ROOT,"reports","experiment_E","tables","TABLE_E09_uniform_boundary_sensitivities.csv"),DataFrame)
nrow(boundary)==1 || error("expected one analytically traced uniform margin crossing")
rho=fill(Float64(boundary.t_star[1]),n)
tables=joinpath(ROOT,"reports","experiment_E","tables")
branch_bus=get(ENV,"EXP_E_ANCHOR_BUS","")
if !isempty(branch_bus)
    b=parse(Int,branch_bus)
    row=only(eachrow(filter(r->Int(r.anchor_bus)==b,
        CSV.read(joinpath(tables,"TABLE_E12_one_SG_branch_summary.csv"),DataFrame))))
    isfinite(row.first_margin_crossing_rho) || error("branch bus $b lacks a margin crossing")
    rho.=1.0;rho[findfirst(==(b),buses)]=Float64(row.first_margin_crossing_rho)
end
path_file=joinpath(tables,isempty(branch_bus) ? "TABLE_E10_analytic_continuation_steps.csv" :
    "TABLE_E13_anchor$(branch_bus)_gain_continuation.csv")

function inspect(rho,Kp,Ki)
    m=CollectiveModel.mixed_jacobian(net,rho,Kp,Ki)
    sp=ClosureSpectrum.physical_spectrum(m.Ared;gauge_tolerance=config["gauge_tolerance"])
    return (model=m,sp=sp)
end

function write_path(rows)
    CSV.write(path_file,DataFrame(rows))
end

rows=NamedTuple[]
function record!(rows,k,rho,Kp,Ki,obs,direction,corrections,event)
    push!(rows,(step=k,event=event,replacement_MW=dot(p,rho),
        retained_SG_MW=dot(p,1 .- rho),critical_real=real(obs.sp.critical),
        critical_imag=imag(obs.sp.critical),spectral_abscissa=obs.sp.spectral_abscissa,
        gauge_detected=obs.sp.gauge_detected,dynamic_dimension=obs.model.n_dynamic,
        gain_corrections=corrections,
        rho=join(rho,";"),Kp=join(Kp,";"),Ki=join(Ki,";"),
        drho=join(direction,";")))
    write_path(rows)
end

step_start=0
if get(ENV,"EXP_E_RESUME","0")=="1" && isfile(path_file)
    old=CSV.read(path_file,DataFrame)
    rows=NamedTuple[NamedTuple(r) for r in eachrow(old)]
    lastrow=old[end,:]
    rho=parse.(Float64,split(String(lastrow.rho),';'))
    Kp=parse.(Float64,split(String(lastrow.Kp),';'))
    Ki=parse.(Float64,split(String(lastrow.Ki),';'))
    step_start=Int(lastrow.step)
    println("RESUMING_FROM_STEP: ",step_start)
end
obs=inspect(rho,Kp,Ki)
step_start==0 && isempty(rows) && record!(rows,0,rho,Kp,Ki,obs,zeros(n),0,
    isempty(branch_bus) ? "UNIFORM_BOUNDARY" : "ONE_ANCHOR_BOUNDARY")
stop_reason="MAX_STEPS"
for stepidx in step_start+1:80
    s=obs.sp.critical
    sens=ClosureSpectrum.pole_sensitivities(net,s,rho,Kp,Ki)
    Jrho=permutedims(real.(sens.ds_drho))
    JK=permutedims(vcat(real.(sens.ds_dKp),real.(sens.ds_dKi)))
    law=AnalyticContinuation.box_limited_direction(p,Jrho,JK,Rinv,Matrix(Qgraph),rho,Kp,Ki,config)
    drho=law.drho;dK=law.dK
    improvement=dot(p,drho)
    if improvement<=1e-8
        stop_reason="NO_POSITIVE_MW_TANGENT";break
    end
    if maximum(abs.(drho))<1e-12
        stop_reason="NO_FEASIBLE_TANGENT";break
    end
    alpha_box=minimum(vcat([drho[i]>1e-12 ? (1-rho[i])/drho[i] : Inf for i in 1:n],
        [drho[i]<-1e-12 ? -rho[i]/drho[i] : Inf for i in 1:n]))
    klo=vcat(fill(0.25*nom.Kp,n),fill(0.25*nom.Ki,n))
    khi=vcat(fill(4*nom.Kp,n),fill(4*nom.Ki,n))
    kcur=vcat(Kp,Ki)
    alpha_gain=minimum([dK[i]>1e-12 ? (khi[i]-kcur[i])/dK[i] :
        dK[i]<-1e-12 ? (klo[i]-kcur[i])/dK[i] : Inf for i in 1:2n])
    alpha=min(0.025,alpha_box,alpha_gain)
    alpha>1e-8 || (stop_reason=alpha==alpha_gain ? "GAIN_BOUND_TANGENT" : "RHO_BOUND_TANGENT";break)
    accepted=false;event="ACTIVE_BOUNDARY"
    for trial in 1:10
        rnew=clamp.(rho.+alpha.*drho,0.0,1.0)
        knewp=clamp.(Kp.+alpha.*dK[1:n],klo[1:n],khi[1:n])
        knewi=clamp.(Ki.+alpha.*dK[n+1:end],klo[n+1:end],khi[n+1:end])
        corrections=0
        good=true
        for _ in 1:8
            o=inspect(rnew,knewp,knewi)
            h=o.sp.spectral_abscissa+sigma
            if abs(h)<1e-8
                if o.sp.gauge_detected
                    accepted=true;obs=o;break
                end
            end
            ss=ClosureSpectrum.pole_sensitivities(net,o.sp.critical,rnew,knewp,knewi)
            Jgain=permutedims(vcat(real.(ss.ds_dKp),real.(ss.ds_dKi)))
            Rcorrect=Matrix(Rinv);corr=zeros(2n);correctable=false
            for _ in 1:2n+1
                auth=(Jgain*Rcorrect*transpose(Jgain))[1]
                if auth<1e-16 || !isfinite(auth);break;end
                corr=vec(-(Rcorrect*transpose(Jgain))*(h/auth))
                kval=vcat(knewp,knewi)
                blocked=[i for i in 1:2n if Rcorrect[i,i]>0 &&
                    (kval[i]+corr[i]<klo[i]-1e-10 || kval[i]+corr[i]>khi[i]+1e-10)]
                if isempty(blocked);correctable=true;break;end
                for i in blocked;Rcorrect[i,i]=0;end
            end
            if !correctable
                good=false;event="LOSS_OF_PLL_AUTHORITY";break
            end
            knewp.=clamp.(knewp.+corr[1:n],klo[1:n],khi[1:n])
            knewi.=clamp.(knewi.+corr[n+1:end],klo[n+1:end],khi[n+1:end])
            corrections+=1
        end
        if accepted && good
            rho=rnew;Kp=knewp;Ki=knewi
            event=alpha==alpha_box ? "RHO_BOUND_EVENT" : alpha==alpha_gain ? "GAIN_BOUND_EVENT" : "ACTIVE_BOUNDARY"
            record!(rows,stepidx,rho,Kp,Ki,obs,drho,corrections,event)
            println("CONTINUATION step=",stepidx," MW=",dot(p,rho),
                " abscissa=",obs.sp.spectral_abscissa," rho_min=",minimum(rho),
                " rho_max=",maximum(rho)," event=",event)
            break
        end
        alpha/=2
        alpha>1e-8 || break
    end
    if !accepted
        stop_reason=event=="ACTIVE_BOUNDARY" ? "CORRECTOR_OR_MODE_SWITCH_FAILURE" : event
        break
    end
end
println("CONTINUATION_STOP_REASON: ",stop_reason)
println("CONTINUATION_FINAL_MW: ",dot(p,rho))
println("CONTINUATION_FINAL_RHO: ",rho)
println("CONTINUATION_FINAL_KP: ",Kp)
println("CONTINUATION_FINAL_KI: ",Ki)
end

main()
