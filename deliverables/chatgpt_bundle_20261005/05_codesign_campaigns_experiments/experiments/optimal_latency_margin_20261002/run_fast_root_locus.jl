using CSV, DataFrames, LinearAlgebra, TOML

module MegaOracle
include(joinpath(@__DIR__,"..","analytical_delay_codesign_mega_20261002","m3_a_trace_integral.jl"))
end

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const MEGA=joinpath(ROOT,"experiments","analytical_delay_codesign_mega_20261002")
BLAS.set_num_threads(1)

function reduced_model(L)
    Ac=L.A0+L.B*transpose(L.C)
    scale=MegaOracle.diagonal_balance(Ac)
    A=(Ac.*transpose(scale))./reshape(scale,:,1)
    B=L.B./reshape(scale,:,1)
    C=transpose(L.C).*transpose(scale)
    f=schur(ComplexF64.(A))
    (;T=f.T,Z=f.Z,left=ComplexF64.(C*f.Z),right=ComplexF64.(adjoint(f.Z)*B),scale)
end

function small_factor(m,s,tau)
    Q=UpperTriangular(s*I-m.T)
    X=Q\m.right;X2=Q\X
    H=m.left*X;H2=m.left*X2
    e=exp.(-s.*tau)
    F=Matrix{ComplexF64}(I,length(tau),length(tau))-Diagonal(e.-1)*H
    Fs=Diagonal(tau.*e)*H+Diagonal(e.-1)*H2
    F,Fs,X,H,e
end

function refine(m,s0,tau;tol=1e-9,maxiter=40)
    s=ComplexF64(s0)
    for iteration in 1:maxiter
        F,Fs,X,H,e=small_factor(m,s,tau)
        decomp=svd(F);sv=decomp.S[end]
        if sv<tol
            z=decomp.V[:,end]
            vb=m.Z*(X*z);vb=m.scale.*vb;vb/=norm(vb)
            u=decomp.U[:,end]
            root_tau_sensitivity=-dot(u,(s*Diagonal(e)*H)*z)/dot(u,Fs*z)
            return (;s,sv,iteration,converged=true,vb,root_tau_sensitivity)
        end
        d=try tr(F\Fs) catch; break end
        abs(d)>1e-13 || break
        ds=-inv(d)
        abs(ds)>3 && (ds*=3/abs(ds))
        s+=ds
        (isfinite(real(s)) && isfinite(imag(s)) && -10000<real(s)<1000) || break
    end
    (;s,sv=NaN,iteration=maxiter,converged=false,vb=ComplexF64[],root_tau_sensitivity=NaN+NaN*im)
end

function run_design(ctx,design_id,design_path,seed_path)
    d=TOML.parsefile(design_path)
    L=MegaOracle.DC.linearization(ctx,Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"]))
    m=reduced_model(L)
    seeds=CSV.read(seed_path,DataFrame)
    rows=NamedTuple[]
    for (family_id,seed) in enumerate(eachrow(seeds))
        s=ComplexF64(seed.root_real+im*seed.root_imag)
        previous_vector=ComplexF64[]
        delay_step_ms=0.1
        for tau_ms in 40.0:-delay_step_ms:20.0
            tau=fill(tau_ms/1000,length(L.Ai))
            r=refine(m,s,tau)
            r.converged || begin
                push!(rows,(;design_id,family_id,tau_ms,root_real=real(r.s),root_imag=imag(r.s),
                    frequency_hz=imag(r.s)/(2pi),small_factor_sigma_min=NaN,
                    overlap_with_previous=NaN,uniform_delay_sensitivity_real=NaN,
                    uniform_delay_sensitivity_imag=NaN,status="TRACKING_FAILED"))
                break
            end
            overlap=isempty(previous_vector) ? 1.0 : abs(dot(previous_vector,r.vb))^2
            push!(rows,(;design_id,family_id,tau_ms,root_real=real(r.s),root_imag=imag(r.s),
                frequency_hz=imag(r.s)/(2pi),small_factor_sigma_min=r.sv,
                overlap_with_previous=overlap,
                uniform_delay_sensitivity_real=real(r.root_tau_sensitivity),
                uniform_delay_sensitivity_imag=imag(r.root_tau_sensitivity),
                status="LOCAL_EXACT_CHARACTERISTIC_CONTINUATION"))
            s=r.s-0.001*delay_step_ms*r.root_tau_sensitivity
            previous_vector=r.vb
        end
        CSV.write(joinpath(@__DIR__,"T08_FAST_MODE_PROVENANCE.csv"),DataFrame(rows))
        println("FAST_BRANCH_DONE ",design_id," family=",family_id,
            " last_tau_ms=",last(rows).tau_ms," last_status=",last(rows).status);flush(stdout)
    end
    rows
end

function main()
    ctx=MegaOracle.DC.R.N.design_context(ROOT)
    allrows=NamedTuple[]
    for (id,path,rootfile) in (("seed_875",joinpath(MEGA,"seed_uniform_875.toml"),
                                    joinpath(MEGA,"M3_REFINED_ROOTS.csv")),
                               ("best_zero_delay_88455",joinpath(MEGA,"M1_ZERO_DELAY_DESIGN.toml"),
                                    joinpath(MEGA,"M3_REFINED_ROOTS_M1_ZERO_DELAY_DESIGN.csv")))
        append!(allrows,run_design(ctx,id,path,rootfile))
        CSV.write(joinpath(@__DIR__,"T08_FAST_MODE_PROVENANCE.csv"),DataFrame(allrows))
    end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
