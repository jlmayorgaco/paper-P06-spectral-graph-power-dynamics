using CSV, DataFrames, LinearAlgebra, Statistics
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_design_e","PhysicalData.jl"))
include(joinpath(ROOT,"src","bnd_design_e","CollectiveModel.jl"))
include(joinpath(ROOT,"src","bnd_design_e","ClosureSpectrum.jl"))

domain=PhysicalData.freeze_design_domain(ROOT)
net=CollectiveModel.frozen_network(ROOT)
buses=sort(collect(keys(net.sg)));n=length(buses)
nom=PhysicalData.nominal_pll_gains()
Kp=fill(nom.Kp,n);Ki=fill(nom.Ki,n)
p=Float64.(domain.generators.SG_dispatch_initial_MW)
M=Matrix{Float64}(CSV.read(joinpath(ROOT,"reports","experiment_C","matrices","C0_M.csv"),DataFrame)[:,2:end])
Lg=Matrix{Float64}(CSV.read(joinpath(ROOT,"reports","experiment_C","matrices","C0_LG_primary.csv"),DataFrame)[:,2:end])
Qgraph=transpose(Lg)*(M\Lg)
symQ=(Qgraph+transpose(Qgraph))/2
minimum(eigvals(Symmetric(symQ)))>-1e-8 || error("graph placement metric is not PSD")
d0=(symQ+1e-6I(n))\p
d0./=maximum(d0)
println("UNCONSTRAINED_GRAPH_DIRECTION: ",d0)

function evaluate(t)
    rho=fill(t,n)
    model=CollectiveModel.mixed_jacobian(net,rho,Kp,Ki)
    sp=ClosureSpectrum.physical_spectrum(model.Ared)
    return (t=t,model=model,sp=sp)
end

rows=NamedTuple[]; evals=Dict{Float64,Any}()
for t in (0.0,0.05,0.1,0.2,0.4,0.6,0.8,0.9,0.95,0.99,0.999)
    e=evaluate(t);evals[t]=e
    push!(rows,(t=t,dimension=e.model.n_dynamic,gauge_detected=e.sp.gauge_detected,
        gauge_abs=e.sp.gauge_detected ? abs(e.sp.gauge_pole) : NaN,
        critical_real=real(e.sp.critical),critical_imag=imag(e.sp.critical),
        spectral_abscissa=e.sp.spectral_abscissa,
        feasible_margin=e.sp.spectral_abscissa<=-0.05))
    println("PATH t=",t," abscissa=",e.sp.spectral_abscissa," critical=",e.sp.critical,
        " gauge=",e.sp.gauge_detected)
end
tables=joinpath(ROOT,"reports","experiment_E","tables")
CSV.write(joinpath(tables,"TABLE_E08_uniform_analytic_path.csv"),DataFrame(rows))

crossings=NamedTuple[]
for j in 1:length(rows)-1
    a,b=rows[j],rows[j+1]
    fa=a.spectral_abscissa+0.05;fb=b.spectral_abscissa+0.05
    fa*fb<0 || continue
    lo,hi=a.t,b.t
    for _ in 1:28
        mid=(lo+hi)/2;em=evaluate(mid)
        if (em.sp.spectral_abscissa+0.05)*fa<=0
            hi=mid
        else
            lo=mid;fa=em.sp.spectral_abscissa+0.05
        end
    end
    tstar=(lo+hi)/2;e=evaluate(tstar)
    s=e.sp.critical
    sens=ClosureSpectrum.pole_sensitivities(net,s,fill(tstar,n),Kp,Ki)
    Jρ=real.(sens.ds_drho);JK=vcat(real.(sens.ds_dKp),real.(sens.ds_dKi))'
    Rinv=Diagonal(vcat(fill(nom.Kp^2,n),fill(nom.Ki^2,n)))
    authority=(JK*Rinv*transpose(JK))[1]
    qpll=authority>1e-16 ? (Jρ*transpose(Jρ))/authority : fill(Inf,n,n)
    Qeff=symQ+qpll+1e-6I(n)
    dρ=Qeff\p;dρ./=maximum(abs.(dρ))
    dK=authority>1e-16 ? -(Rinv*transpose(JK))*(dot(Jρ,dρ)/authority) : fill(NaN,2n)
    push!(crossings,(t_star=tstar,critical_real=real(s),critical_imag=imag(s),
        closure_sigma_min=sens.closure_sigma_min,closure_residual=sens.closure_residual,
        authority=authority,graph_PSD_minimum=minimum(eigvals(Symmetric(symQ))),
        d_rho=join(dρ,";"),d_Kp=join(dK[1:n],";"),d_Ki=join(dK[n+1:end],";")))
    println("BOUNDARY t=",tstar," pole=",s," closure_sigma=",sens.closure_sigma_min,
        " authority=",authority," direction=",dρ)
end
CSV.write(joinpath(tables,"TABLE_E09_uniform_boundary_sensitivities.csv"),DataFrame(crossings))
println("MARGIN_CROSSING_COUNT: ",length(crossings))
