using CSV, DataFrames, LinearAlgebra, TOML

module Action
include(joinpath(@__DIR__, "..", "analytical_delay_codesign_mega_20261002", "q56_validate_action_space.jl"))
end

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(@__DIR__,"L1_ACTION_SPACE_VALIDATION.csv")
BLAS.set_num_threads(1)

function probe_values(row)
    [Float64(row[Symbol("Kp_",i)]) for i in 1:10],
    [Float64(row[Symbol("Ki_",i)]) for i in 1:10]
end

function main()
    ctx=Action.N.design_context(ROOT)
    probes=CSV.read(joinpath(@__DIR__,"L1_FROZEN_GAIN_PROBES.csv"),DataFrame)
    z=TOML.parsefile(joinpath(@__DIR__,"designs","Z_zero_delay_tuned.toml"))
    n=TOML.parsefile(joinpath(@__DIR__,"designs","N_nominal.toml"))
    rho=Float64.(z["rho"])
    rho==Float64.(n["rho"]) || error("not the same replacement vector")
    kpref=sqrt.(Float64.(z["Kp"]).*Float64.(n["Kp"]))
    kiref=sqrt.(Float64.(z["Ki"]).*Float64.(n["Ki"]))
    mr=Action.N.descriptor(ctx,rho,kpref,kiref)
    samples=((0.0,0.017),(37.38,5.26),(39.38,5.64),(40.0,5.26))
    rows=NamedTuple[]
    for row in eachrow(probes)
        id=String(row.id);kp,ki=probe_values(row)
        mp=Action.N.descriptor(ctx,rho,kp,ki)
        pairs=startswith(id,"random_") ? (samples[2],samples[3]) : samples
        for (tau_ms,freq_hz) in pairs
            s=-0.05+2pi*im*freq_hz;tau=fill(tau_ms/1000,10)
            p0=Action.descriptor_dde_pencil(ctx,rho,kpref,kiref,tau,s,mr)
            p1=Action.descriptor_dde_pencil(ctx,rho,kp,ki,tau,s,mp)
            nd=mr.n_dynamic;na=mr.n_algebraic
            U=zeros(ComplexF64,nd+na,10)
            H=zeros(ComplexF64,10,nd+na)
            for g in p0.geom
                U[1:nd,g.i]=(kp[g.i]-kpref[g.i]).*g.bp.+(ki[g.i]-kiref[g.i]).*g.bi
                H[g.i,:]=-exp(-s*tau[g.i]).*vcat(g.cx,g.cy)
            end
            difference=p1.T-p0.T
            reconstruction=U*H
            absolute_err=norm(difference-reconstruction)
            relative_err=absolute_err/max(norm(difference),eps())
            sig=svdvals(difference)
            numerical_rank=count(x->x>1e-9*maximum(sig),sig)
            gain_rank_upper_bound=rank(U;rtol=1e-9)
            absdet_ref,phase_ref=logabsdet(p0.T)
            absdet_full,phase_full=logabsdet(p1.T)
            small=Matrix{ComplexF64}(I,10,10)+H*(p0.T\U)
            absdet_small,phase_small=logabsdet(small)
            logdet_error=abs(absdet_full-(absdet_ref+absdet_small))
            phase_error=abs(angle(phase_full/(phase_ref*phase_small)))
            status=relative_err<1e-10 && numerical_rank<=10 && logdet_error<1e-8 && phase_error<1e-8 ?
                "EXACT_IDENTITY_NUMERICALLY_VERIFIED" : "FAIL_OR_ILL_CONDITIONED"
            push!(rows,(;design_id=id,tau_ms,freq_hz,full_descriptor_dimension=nd+na,
                fixed_rho_gain_action_columns=10,compressed_U_rank=gain_rank_upper_bound,
                numerical_difference_rank=numerical_rank,relative_reconstruction_error=relative_err,
                absolute_reconstruction_error=absolute_err,log_abs_determinant_lemma_error=logdet_error,
                determinant_phase_error_rad=phase_error,schur_err_ref=p0.schurerr,
                schur_err_probe=p1.schurerr,status,
                feasibility_scope=startswith(id,"random_") ? "GAIN_BOX_ONLY_EVENTS_NOT_SCREENED" :
                    "STORED_DESIGN_GAIN_BOX_ZERO_DELAY_EVENTS_EVALUATED_SEPARATELY"))
            CSV.write(OUT,DataFrame(rows))
            println("L1_ACTION ",id," tau_ms=",tau_ms," err=",relative_err," rank=",numerical_rank," status=",status);flush(stdout)
        end
    end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
