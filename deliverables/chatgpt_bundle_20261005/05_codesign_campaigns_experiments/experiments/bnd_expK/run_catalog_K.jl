using LinearAlgebra, CSV, DataFrames, Statistics, SHA, TOML
BLAS.set_num_threads(1)
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_design_k","BNDDesignK.jl"))
using .BNDDesignK
const OUT=joinpath(ROOT,"reports","experiment_K")
const TABLES=joinpath(OUT,"tables")
mkpath(TABLES)
ctx=design_context(ROOT)
ctx.buses==collect(30:39) || error("unexpected replaceable set")
abs(sum(ctx.power)-5402.761089978776)<1e-6 || error("SG MW integrity failure")

# A support is assessed at an interior reference eps=0.5 for every present SG.
# This gives one exact spectrum per architecture, not a feasibility certificate.
catalog=NamedTuple[]; branches=NamedTuple[]
for mask in 0:(2^length(ctx.buses)-1)
    e=support_epsilon(ctx,mask,0.5)
    ev=evaluate(ctx,e;vectors=true)
    active=ev.active
    push!(catalog,(support_mask=mask,support_buses=join(support_buses(ctx,mask),":"),
        support_size=count(>(0),e),dynamic_dimension=ev.model.n_dynamic,
        quotient_dimension=length(ev.lambda),reference_epsilon=0.5,
        reference_retained_sg_MW=ev.retained_mw,
        reference_alpha_s_inv=ev.alpha,
        reference_margin_feasible=ev.alpha<=-BNDDesignK.SIGMA_REQUIRED,
        active_pole_count=length(active),gauge_residual=ev.gauge_residual,
        continuous_status="REFERENCE_POINT_ONLY",
        infeasibility_certified=false,branch_exhausted=false))
    for j in active
        fam=pole_family(ev,j)
        push!(branches,(support_mask=mask,support_buses=join(support_buses(ctx,mask),":"),
            reference_epsilon=0.5,pole_index=j,real_s_inv=real(ev.lambda[j]),
            imag_s_inv=imag(ev.lambda[j]),spectral_gap_s_inv=ev.alpha-real(ev.lambda[j]),
            family=fam.family,sg_state_energy=fam.sg_energy,
            gfl_state_energy=fam.gfl_energy,
            classification_scope="participation_at_reference_point"))
    end
    mask%128==0 && (println("K catalog mask $mask / 1023");flush(stdout))
end
CSV.write(joinpath(TABLES,"TABLE_K01_support_catalog.csv"),DataFrame(catalog))
CSV.write(joinpath(TABLES,"TABLE_K02_branch_catalog.csv"),DataFrame(branches))

# Validate fixed-architecture ordinary eigenvalue derivatives at a regular
# all-bus interior point, retaining conjugate-sensitive complex comparisons.
v=evaluate(ctx,fill(0.5,length(ctx.buses));vectors=true)
j=v.active[1]
deriv=NamedTuple[]
for group in (:epsilon,:kp,:ki), i in eachindex(ctx.buses)
    h=group==:epsilon ? 1e-5 : 1e-4
    d=derivative_check(ctx,v,j,group,i;h=h)
    push!(deriv,(bus=ctx.buses[i],group=String(group),pole_real=real(v.lambda[j]),
        pole_imag=imag(v.lambda[j]),analytic_real=real(d.analytical),
        analytic_imag=imag(d.analytical),central_real=real(d.finite_difference),
        central_imag=imag(d.finite_difference),absolute_error=d.absolute_error,
        step=h,support_mask=1023))
end
CSV.write(joinpath(TABLES,"TABLE_K03_derivative_validation.csv"),DataFrame(deriv))

# A one-sided SG block has its own open-loop poles. The coefficient below is
# measured within the positive-epsilon architecture, never differentiated
# through the absent/present transition.
reentry=NamedTuple[]
for (i,b) in enumerate(ctx.buses)
    mu=eigvals(ctx.net.sg[b].J.A)
    ee=support_epsilon(ctx,1<<(i-1),1e-5)
    ev1=evaluate(ctx,ee)
    ee[i]=5e-6
    ev2=evaluate(ctx,ee)
    for (k,z) in enumerate(mu)
        p1=ev1.lambda[argmin(abs.(ev1.lambda.-z))]
        p2=ev2.lambda[argmin(abs.(ev2.lambda.-z))]
        coeff=(p1-p2)/5e-6
        push!(reentry,(bus=b,sg_open_loop_index=k,mu_real=real(z),mu_imag=imag(z),
            b_real=real(coeff),b_imag=imag(coeff),lambda_1e5_real=real(p1),
            lambda_1e5_imag=imag(p1),lambda_5e6_real=real(p2),
            lambda_5e6_imag=imag(p2),limit_match_error=abs(p2-z),
            status="ONE_SIDED_TWO_POINT_ESTIMATE_NOT_RIGOROUS_REMAINDER"))
    end
end
CSV.write(joinpath(TABLES,"TABLE_K04_reentry_branch_coefficients.csv"),DataFrame(reentry))

# Detect every feasible component seen on a preregistered log+linear mesh;
# refine observed sign-changing endpoints by deterministic bisection. This
# cannot exclude narrow unsampled components and is labelled accordingly.
intervals=NamedTuple[]
mesh=sort(unique(vcat(10.0.^range(-8,-2,length=25),collect(range(0.01,1.0,length=100)))))
function margin_at(i,x)
    ee=support_epsilon(ctx,1<<(i-1),x)
    evaluate(ctx,ee).alpha+BNDDesignK.SIGMA_REQUIRED
end
for (i,b) in enumerate(ctx.buses)
    vals=[margin_at(i,x) for x in mesh]
    inside=vals .<= 0
    runstart=nothing
    components=Tuple{Float64,Float64}[]
    for k in eachindex(mesh)
        if inside[k] && runstart===nothing
            lo=mesh[k]
            if k>1
                a=mesh[k-1]; z=lo
                for _ in 1:35
                    m=(a+z)/2
                    if margin_at(i,m)<=0; z=m else a=m end
                end
                lo=z
            end
            runstart=lo
        end
        if runstart!==nothing && (k==length(mesh) || !inside[k+1])
            hi=mesh[k]
            if k<length(mesh)
                a=hi; z=mesh[k+1]
                for _ in 1:35
                    m=(a+z)/2
                    if margin_at(i,m)<=0; a=m else z=m end
                end
                hi=a
            end
            push!(components,(runstart,hi));runstart=nothing
        end
    end
    if isempty(components)
        push!(intervals,(bus=b,component=0,epsilon_lo=NaN,epsilon_hi=NaN,
            alpha_min_sampled=minimum(vals)-BNDDesignK.SIGMA_REQUIRED,
            status="NO_FEASIBLE_MESH_POINT_NOT_PROVEN_EMPTY"))
    else
        for (k,(lo,hi)) in enumerate(components)
            push!(intervals,(bus=b,component=k,epsilon_lo=lo,epsilon_hi=hi,
                alpha_min_sampled=minimum(vals)-BNDDesignK.SIGMA_REQUIRED,
                status="OBSERVED_COMPONENT_BOUNDARIES_BISECTED"))
        end
    end
end
CSV.write(joinpath(TABLES,"TABLE_K05_single_bus_feasible_intervals.csv"),DataFrame(intervals))

println("K_CATALOG_SUPPORTS: ",nrow(DataFrame(catalog)))
println("K_REFERENCE_FEASIBLE_SUPPORTS: ",count(x->x.reference_margin_feasible,catalog))
println("K_MAX_DERIVATIVE_ABS_ERROR: ",maximum(x->x.absolute_error,deriv))
