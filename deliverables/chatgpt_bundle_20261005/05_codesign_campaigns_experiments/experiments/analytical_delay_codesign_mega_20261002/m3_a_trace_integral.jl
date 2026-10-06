using CSV, DataFrames, LinearAlgebra, TOML

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DC = include(joinpath(ROOT, "experiments", "delay_dressed_replacement_frontier_20261002", "DelayCharacteristic.jl"))
BLAS.set_num_threads(1)

# A separate implementation of the logarithmic-derivative contour integral.
# The input matrices remain Float64, so this is numerical evidence, not a
# directed-rounding certificate.
function diagonal_balance(A; sweeps=100)
    d = ones(Float64, size(A,1))
    for _ in 1:sweeps
        changed = false
        for i in eachindex(d)
            r = sum(abs(A[i,j])*d[j]/d[i] for j in eachindex(d) if i != j)
            c = sum(abs(A[j,i])*d[i]/d[j] for j in eachindex(d) if i != j)
            (r == 0 || c == 0) && continue
            f = 2.0^round(Int, 0.5*log2(r/c))
            if f != 1 && isfinite(d[i]*f)
                d[i] *= f
                changed = true
            end
        end
        !changed && break
    end
    d
end

function bound_radius(A,B,C,tau,gamma)
    nA = opnorm(A,2)
    cb = opnorm(C*B,2)
    ca = opnorm(C*A,2)
    bn = opnorm(B,2)
    ef = maximum(1+exp(-gamma*t) for t in tau)
    bound(r) = ef*(cb/r + ca*bn/(r*(r-nA)))
    lo = nA + max(1.0,1e-10*nA)
    hi = max(2lo,lo+1)
    while bound(hi) >= 1 && hi < 1e12
        hi *= 2
    end
    bound(hi) < 1 || error("no finite tail bound")
    for _ in 1:80
        mid = (lo+hi)/2
        bound(mid) < 1 ? (hi=mid) : (lo=mid)
    end
    hi + max(1e-6,1e-12*hi)
end

function trace_count(L,tau;gamma=-0.05,step=20.0,local_tolerance=5e-7,
        max_evaluations=200000,max_depth=15)
    Ac = L.A0+L.B*transpose(L.C)
    C = transpose(L.C)
    scale = diagonal_balance(Ac)
    A = (Ac .* transpose(scale)) ./ reshape(scale,:,1)
    B = L.B ./ reshape(scale,:,1)
    Cb = C .* transpose(scale)
    radius = bound_radius(A,B,Cb,tau,gamma)
    sf = schur(ComplexF64.(A))
    T,Z = sf.T,sf.Z
    left = Cb*Z
    right = adjoint(Z)*B
    reference_count = count(z->real(z)>gamma,diag(T))
    m = length(tau)
    ident = Matrix{ComplexF64}(I,m,m)
    cache = Dict{ComplexF64,Tuple{ComplexF64,Float64,ComplexF64}}()
    minsv = Inf
    depth_reached = 0

    function integrand(s)
        if haskey(cache,s)
            return cache[s][1]
        end
        length(cache) < max_evaluations || error("contour evaluation budget exceeded")
        Q = UpperTriangular(s*I-T)
        X = Q\right
        X2 = Q\X
        H = left*X
        H2 = left*X2
        e = exp.(-s.*tau)
        small = ident-Diagonal(e .- 1)*H
        slope = Diagonal(tau.*e)*H + Diagonal(e .- 1)*H2
        value = tr(small\slope)
        sv = minimum(svdvals(small))
        minsv = min(minsv,sv)
        isfinite(real(value)) && isfinite(imag(value)) && isfinite(sv) || error("nonfinite contour value")
        cache[s] = (value,sv,det(small))
        value
    end

    simpson(a,b,fa,fm,fb) = (b-a)*(fa+4fm+fb)/6
    function recurse(a,b,fa,fm,fb,coarse,depth)
        depth_reached = max(depth_reached,depth)
        mid = (a+b)/2
        q1,q3 = (a+mid)/2,(mid+b)/2
        fq1,fq3 = integrand(q1),integrand(q3)
        refined = simpson(a,mid,fa,fq1,fm)+simpson(mid,b,fm,fq3,fb)
        err = abs(refined-coarse)/15
        if err <= local_tolerance
            return refined+(refined-coarse)/15,err
        end
        depth < max_depth || error("contour quadrature unresolved at depth $(max_depth); err=$(err)")
        v1,e1 = recurse(a,mid,fa,fq1,fm,simpson(a,mid,fa,fq1,fm),depth+1)
        v2,e2 = recurse(mid,b,fm,fq3,fb,simpson(mid,b,fm,fq3,fb),depth+1)
        v1+v2,e1+e2
    end

    corners = ComplexF64[gamma-im*radius,radius-im*radius,radius+im*radius,
        gamma+im*radius,gamma-im*radius]
    integral = 0.0+0.0im
    errsum = 0.0
    for edge in 1:4
        a0,b0 = corners[edge],corners[edge+1]
        nseg = max(2,ceil(Int,abs(b0-a0)/step))
        for j in 0:nseg-1
            a = a0+(b0-a0)*j/nseg
            b = a0+(b0-a0)*(j+1)/nseg
            mid = (a+b)/2
            fa,fm,fb = integrand(a),integrand(mid),integrand(b)
            value,err = recurse(a,b,fa,fm,fb,simpson(a,b,fa,fm,fb),0)
            integral += value
            errsum += err
        end
        println("M3_A_EDGE_DONE edge=",edge," tau_ms=",1000*maximum(tau)," evals=",length(cache));flush(stdout)
    end
    phase_total = 0.0
    largest_phase_step = 0.0
    for edge in 1:4
        a0,b0 = corners[edge],corners[edge+1]
        nodes = Tuple{Float64,ComplexF64}[]
        for (s,(_,_,d)) in cache
            t = (s-a0)/(b0-a0)
            if abs(imag(t)) < 1e-10 && -1e-10 <= real(t) <= 1+1e-10
                push!(nodes,(real(t),d))
            end
        end
        sort!(nodes,by=first)
        for j in 1:length(nodes)-1
            step_angle = angle(nodes[j+1][2]/nodes[j][2])
            phase_total += step_angle
            largest_phase_step = max(largest_phase_step,abs(step_angle))
        end
    end
    phase_estimate = reference_count+phase_total/(2pi)
    estimate = reference_count+integral/(2pi*im)
    nearest = round(Int,real(estimate))
    (;gamma,radius,reference_count,estimate,phase_estimate,largest_phase_step,nearest,integer_distance=abs(estimate-nearest),
        quadrature_error_estimate=errsum/(2pi),min_sigma_small=minsv,
        evaluations=length(cache),max_refinement_depth=depth_reached,
        original_Ac_norm=opnorm(Ac,2),balanced_Ac_norm=opnorm(A,2))
end

function main()
    design_path = length(ARGS)>=2 ? ARGS[2] : joinpath(@__DIR__,"seed_uniform_875.toml")
    seed = TOML.parsefile(design_path)
    design_id = splitext(basename(design_path))[1]
    ctx = DC.R.N.design_context(ROOT)
    L = DC.linearization(ctx,Float64.(seed["rho"]),Float64.(seed["Kp"]),Float64.(seed["Ki"]))
    delays = isempty(ARGS) ? (20.0,40.0) : (parse(Float64,ARGS[1]),)
    rows = NamedTuple[]
    for tau_ms in delays
        tau = fill(tau_ms/1000,length(L.Ai))
        println("M3_A_START tau_ms=",tau_ms);flush(stdout)
        row = try
            r = trace_count(L,tau)
            (;design_id,tau_ms,status=(r.integer_distance<0.01 && r.quadrature_error_estimate<0.01 &&
                    abs(r.phase_estimate-r.nearest)<0.01 && r.largest_phase_step<pi/2 ?
                    "STABLE_NUMERICAL_INTEGER" : "INDETERMINATE"),
                trace_count_real=real(r.estimate),trace_count_imag=imag(r.estimate),phase_count=r.phase_estimate,
                largest_phase_step_rad=r.largest_phase_step,nearest_integer=r.nearest,
                integer_distance=r.integer_distance,quadrature_error_estimate=r.quadrature_error_estimate,
                min_sigma_small=r.min_sigma_small,contour_radius=r.radius,evaluations=r.evaluations,
                max_refinement_depth=r.max_refinement_depth,reference_pole_count=r.reference_count,
                original_Ac_norm=r.original_Ac_norm,balanced_Ac_norm=r.balanced_Ac_norm,error="")
        catch err
            (;design_id,tau_ms,status="INDETERMINATE",trace_count_real=NaN,trace_count_imag=NaN,
                phase_count=NaN,largest_phase_step_rad=NaN,
                nearest_integer=missing,integer_distance=NaN,quadrature_error_estimate=NaN,
                min_sigma_small=NaN,contour_radius=NaN,evaluations=missing,
                max_refinement_depth=missing,reference_pole_count=missing,
                original_Ac_norm=NaN,balanced_Ac_norm=NaN,error=sprint(showerror,err))
        end
        push!(rows,row)
        outfile=design_id=="seed_uniform_875" ? "M3_TRACE_INTEGRAL_ATTEMPT.csv" :
            "M3_TRACE_INTEGRAL_$(design_id).csv"
        output_path=joinpath(@__DIR__,outfile)
        prior=isfile(output_path) ? CSV.read(output_path,DataFrame) : DataFrame()
        if nrow(prior)>0
            if !(:design_id in propertynames(prior))
                prior.design_id=fill(design_id,nrow(prior))
            end
            prior=prior[.!((prior.design_id .== design_id) .& (prior.tau_ms .== tau_ms)),:]
            CSV.write(output_path,vcat(prior,DataFrame(rows);cols=:union))
        else
            CSV.write(output_path,DataFrame(rows))
        end
        println("M3_A_RESULT ",row);flush(stdout)
    end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
