using CSV, DataFrames, LinearAlgebra
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","experiment_N")
include(joinpath(ROOT,"src","bnd_model_expN","PDExactDesignN.jl"))
using .PDExactDesignN

function continuous_pole(ctx,rho,kp,ki,reference)
    spec=PDExactDesignN.spectrum(ctx,rho,kp,ki)
    spec.lambda[argmin(abs.(spec.lambda .- reference))]
end

function main()
    ctx=PDExactDesignN.design_context(ROOT)
    rows=NamedTuple[]
    for label in ("random_1","random_2","random_3")
        path=joinpath(OUT,"validation_cases",label,"case_definition.csv")
        t=CSV.read(path,DataFrame)
        rho=Float64.(t.rho);kp=Float64.(t.Kp);ki=Float64.(t.Ki)
        deriv=PDExactDesignN.simple_mode_sensitivities(ctx,rho,kp,ki)
        base=deriv.lambda
        all_poles=PDExactDesignN.spectrum(ctx,rho,kp,ki).lambda
        gap=minimum(abs.(all_poles[abs.(all_poles .- base) .> 1e-10] .- base))
        for (kind,basevec,analytical) in (("rho",rho,deriv.rho),
                                          ("Kp",kp,deriv.Kp),
                                          ("Ki",ki,deriv.Ki))
            for i in 1:10
                # Scale h against the modal derivative: large sensitivities
                # need a smaller h to control O(h^2) curvature; small ones
                # need a larger h to suppress eigensolver roundoff. The
                # observed step sweep is archived by gradient_step_audit.jl.
                h=if kind=="rho"
                    a=abs(analytical[i])
                    a>0.5 ? 3e-5 : a>0.005 ? 3e-4 : 1e-3
                else
                    1e-5*max(1.0,abs(basevec[i]))
                end
                plus=copy(basevec);minus=copy(basevec)
                plus[i]+=h;minus[i]-=h
                p=kind=="rho" ? continuous_pole(ctx,plus,kp,ki,base) :
                    kind=="Kp" ? continuous_pole(ctx,rho,plus,ki,base) :
                                 continuous_pole(ctx,rho,kp,plus,base)
                m=kind=="rho" ? continuous_pole(ctx,minus,kp,ki,base) :
                    kind=="Kp" ? continuous_pole(ctx,rho,minus,ki,base) :
                                 continuous_pole(ctx,rho,kp,minus,base)
                fd=(p-m)/(2h)
                err=abs(fd-analytical[i])
                rel=err/max(abs(fd),abs(analytical[i]),1e-8)
                push!(rows,(;case=label,bus=i+29,parameter=kind,
                    lambda_real=real(base),lambda_imag=imag(base),
                    isolated_pole_gap=gap,pole_condition=deriv.condition,
                    analytical_real=real(analytical[i]),
                    analytical_imag=imag(analytical[i]),
                    finite_difference_real=real(fd),
                    finite_difference_imag=imag(fd),step=h,
                    absolute_error=err,relative_error=rel,
                    equilibrium_state_total_derivative=0.0,
                    status=(rel<1e-5 || err<1e-7) ? "PASS" : "FAIL_DERIVATIVES"))
            end
        end
        println("N_GRADIENT_CASE ",label," lambda=",base," gap=",gap,
            " condition=",deriv.condition);flush(stdout)
    end
    CSV.write(joinpath(OUT,"TABLE_N08_derivative_validation.csv"),DataFrame(rows))
    println("N_GRADIENT_DONE rows=",length(rows)," max_abs=",
        maximum(r.absolute_error for r in rows)," max_rel=",
        maximum(r.relative_error for r in rows)," failures=",
        count(r->r.status!="PASS",rows))
end

main()
