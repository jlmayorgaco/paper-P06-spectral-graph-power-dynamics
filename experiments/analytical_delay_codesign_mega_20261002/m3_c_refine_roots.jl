using CSV, DataFrames, LinearAlgebra, TOML

module TraceCount
include("m3_a_trace_integral.jl")
end

const ROOT = normpath(joinpath(@__DIR__,"..",".."))
const DC = TraceCount.DC
BLAS.set_num_threads(1)

function small_factor(A,B,C,tau,s)
    X = (s*I-A)\B
    X2 = (s*I-A)\X
    H = C*X
    H2 = C*X2
    e = exp.(-s.*tau)
    F = Matrix{ComplexF64}(I,length(tau),length(tau))-Diagonal(e.-1)*H
    Fs = Diagonal(tau.*e)*H+Diagonal(e.-1)*H2
    F,Fs
end

function attempt(A,B,C,tau,s0)
    s=ComplexF64(s0)
    for k in 1:50
        F,Fs=small_factor(A,B,C,tau,s)
        sv=minimum(svdvals(F))
        if sv<1e-9
            return (;s,sv,iterations=k-1,converged=true)
        end
        deriv=try tr(F\Fs) catch; break end
        (isfinite(real(deriv)) && isfinite(imag(deriv)) && abs(deriv)>1e-12) || break
        ds=-inv(deriv)
        abs(ds)>20 && (ds*=20/abs(ds))
        s+=ds
        (isfinite(real(s)) && isfinite(imag(s)) && -100<real(s)<1000) || break
    end
    (;s,sv=NaN,iterations=50,converged=false)
end

function main()
    design_path=length(ARGS)>=2 ? ARGS[2] : joinpath(@__DIR__,"seed_uniform_875.toml")
    design_id=splitext(basename(design_path))[1]
    seed=TOML.parsefile(design_path)
    ctx=DC.R.N.design_context(ROOT)
    L=DC.linearization(ctx,Float64.(seed["rho"]),Float64.(seed["Kp"]),Float64.(seed["Ki"]))
    Ac=L.A0+L.B*transpose(L.C)
    scale=TraceCount.diagonal_balance(Ac)
    A=ComplexF64.((Ac.*transpose(scale))./reshape(scale,:,1))
    B=ComplexF64.(L.B./reshape(scale,:,1))
    C=ComplexF64.(transpose(L.C).*transpose(scale))
    tau_ms=isempty(ARGS) ? 40.0 : parse(Float64,ARGS[1])
    tau=fill(tau_ms/1000,length(L.Ai))
    starts=Tuple{Float64,ComplexF64}[]
    omega_grid=vcat(collect(0.0:0.25:150.0),collect(152.0:2.0:2000.0))
    for re in (-0.05,0.0,0.2,0.5,1.0,2.0,5.0,10.0,25.0)
        vals=Float64[]
        for om in omega_grid
            F,_=small_factor(A,B,C,tau,re+im*om)
            push!(vals,minimum(svdvals(F)))
        end
        push!(starts,(vals[1],re+0im))
        for j in 2:length(vals)-1
            vals[j]<vals[j-1] && vals[j]<vals[j+1] && vals[j]<1.0 &&
                push!(starts,(vals[j],re+im*omega_grid[j]))
        end
    end
    suffix=design_id=="seed_uniform_875" ? "" : "_$(design_id)"
    CSV.write(joinpath(@__DIR__,"M3_ROOT_DISCOVERY_SEEDS$(suffix).csv"),
        DataFrame([(;design_id,tau_ms,score,initial_real=real(s),initial_imag=imag(s)) for (score,s) in starts]))
    sort!(starts,by=first)
    roots=ComplexF64[]
    rows=NamedTuple[]
    for (score,s0) in first(starts,min(length(starts),200))
        r=attempt(A,B,C,tau,s0)
        r.converged || continue
        real(r.s)>-0.05 || continue
        any(abs(r.s-z)<1e-5 for z in roots) && continue
        push!(roots,r.s)
        D=DC.delta_matrix(L,r.s,tau)
        svdD=svd(D)
        residual=minimum(svdD.S)/max(opnorm(D),1)
        push!(rows,(;design_id,tau_ms,initial_real=real(s0),initial_imag=imag(s0),
            root_real=real(r.s),root_imag=imag(r.s),frequency_Hz=imag(r.s)/(2pi),
            small_factor_sigma_min=r.sv,full_characteristic_normalized_sigma_min=residual,
            iterations=r.iterations,scan_score=score,status="EXACT_CHARACTERISTIC_REFINED_NUMERICALLY"))
        println("M3_C_ROOT ",last(rows));flush(stdout)
    end
    CSV.write(joinpath(@__DIR__,"M3_REFINED_ROOTS$(suffix).csv"),DataFrame(rows))
    println("M3_C_TOTAL_POSITIVE_IMAG=",length(rows));flush(stdout)
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
