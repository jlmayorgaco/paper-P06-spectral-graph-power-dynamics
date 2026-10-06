using CSV, DataFrames, LinearAlgebra
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_design_e","PhysicalData.jl"))
include(joinpath(ROOT,"src","bnd_design_e","CollectiveModel.jl"))
include(joinpath(ROOT,"src","bnd_design_e","ClosureSpectrum.jl"))

function main()
    domain=PhysicalData.freeze_design_domain(ROOT)
    net=CollectiveModel.frozen_network(ROOT)
    buses=sort(collect(keys(net.sg)));n=length(buses)
    nom=PhysicalData.nominal_pll_gains()
    Kp=fill(nom.Kp,n);Ki=fill(nom.Ki,n)
    p=Float64.(domain.generators.SG_dispatch_initial_MW)
    sigma=0.05
    function evaluate(j,t)
        rho=ones(n);rho[j]=t
        model=CollectiveModel.mixed_jacobian(net,rho,Kp,Ki)
        sp=ClosureSpectrum.physical_spectrum(model.Ared)
        return (rho=rho,model=model,sp=sp,margin=sp.spectral_abscissa+sigma)
    end
    branches=NamedTuple[];steps=NamedTuple[]
    for (j,b) in enumerate(buses)
        e0=evaluate(j,0.0)
        push!(steps,(anchor_bus=b,rho_anchor=0.0,margin=e0.margin,
            critical_real=real(e0.sp.critical),critical_imag=imag(e0.sp.critical),
            dimension=e0.model.n_dynamic))
        best_t=e0.margin<=0 ? 0.0 : NaN
        crossing_t=NaN;crossing_pole=NaN+NaN*im
        previous_t=0.0;previous_m=e0.margin
        for t in (0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.975,0.99,0.999)
            e=evaluate(j,t)
            push!(steps,(anchor_bus=b,rho_anchor=t,margin=e.margin,
                critical_real=real(e.sp.critical),critical_imag=imag(e.sp.critical),
                dimension=e.model.n_dynamic))
            if e.margin<=0;best_t=t;end
            if isfinite(best_t) && previous_m<=0 && e.margin>0 && !isfinite(crossing_t)
                lo,hi=previous_t,t
                for _ in 1:24
                    mid=(lo+hi)/2;em=evaluate(j,mid)
                    if em.margin<=0;lo=mid;else;hi=mid;end
                end
                crossing_t=(lo+hi)/2;ee=evaluate(j,crossing_t)
                crossing_pole=ee.sp.critical
            end
            previous_t=t;previous_m=e.margin
        end
        tstar=isfinite(crossing_t) ? crossing_t : best_t
        mw=isfinite(tstar) ? sum(p)-p[j]*(1-tstar) : NaN
        push!(branches,(anchor_bus=b,initial_SG_MW=p[j],
            feasible_at_full_SG=e0.margin<=0,
            first_margin_crossing_rho=crossing_t,
            largest_sampled_feasible_rho=best_t,
            first_crossing_pole_real=real(crossing_pole),
            first_crossing_pole_imag=imag(crossing_pole),
            replacement_MW_at_first_crossing=mw,
            status=isfinite(crossing_t) ? "LOCAL_BRANCH_MARGIN_ROOT" :
                isfinite(best_t) ? "FEASIBLE_SAMPLES_NO_CROSSING" : "NO_FEASIBLE_SAMPLE"))
        println("ANCHOR bus=",b," full_SG_margin=",e0.margin,
            " first_crossing_rho=",crossing_t," MW=",mw)
    end
    tables=joinpath(ROOT,"reports","experiment_E","tables")
    CSV.write(joinpath(tables,"TABLE_E11_one_SG_branch_steps.csv"),DataFrame(steps))
    CSV.write(joinpath(tables,"TABLE_E12_one_SG_branch_summary.csv"),DataFrame(branches))
    finite=[r for r in branches if isfinite(r.replacement_MW_at_first_crossing)]
    if !isempty(finite)
        best=argmax([r.replacement_MW_at_first_crossing for r in finite])
        println("BEST_ONE_SG_NOMINAL_BRANCH: ",finite[best])
    end
end

main()
