using CSV, DataFrames, LinearAlgebra, Random, TOML, Statistics

module Q56
include(joinpath(@__DIR__,"q56_validate_action_space.jl"))
end

const ROOT = normpath(joinpath(@__DIR__,"..",".."))
BLAS.set_num_threads(1)

function factor_difference(p0,p1,kp0,ki0,kp1,ki1,tau,s)
    dT = p1.T-p0.T
    n,ny = p0.n,p0.ny
    U = zeros(ComplexF64,n+ny,30)
    H = zeros(ComplexF64,30,n+ny)
    for g in p0.geom
        db = (kp1[g.i]-kp0[g.i]).*g.bp .+ (ki1[g.i]-ki0[g.i]).*g.bi
        U[1:n,g.i] = db
        H[g.i,:] .= -exp(-s*tau[g.i]).*vcat(g.cx,g.cy)
    end
    for i in 1:10
        bus = 29+i
        for q in 1:2
            row = n+2bus-2+q
            col = 10+2i-2+q
            U[row,col] = 1.0
            H[col,:] = dT[row,:]
        end
    end
    err = norm(dT-U*H)/max(norm(dT),eps())
    sv = svdvals(dT)
    numerical_rank = maximum(sv)==0 ? 0 : count(x->x>1e-9*maximum(sv),sv)
    small = Matrix{ComplexF64}(I,30,30)+H*(p0.T\U)
    ld0,ph0 = logabsdet(p0.T)
    ld1,ph1 = logabsdet(p1.T)
    ldsmall,phsmall = logabsdet(small)
    logabs_error = abs(ld1-ld0-ldsmall)
    phase_error = abs(angle(ph1/(ph0*phsmall)))
    (;err,numerical_rank,logabs_error,phase_error,dimension=n+ny)
end

function draw_design(rng,seed_rho,seed_kp,seed_ki,kind)
    rho = [0.80+0.15*rand(rng) for _ in 1:10]
    kp = [seed_kp[i]*exp(log(0.25)+rand(rng)*log(16.0)) for i in 1:10]
    ki = [seed_ki[i]*exp(log(0.25)+rand(rng)*log(16.0)) for i in 1:10]
    if kind == "GAIN_BOUND_VIOLATION"
        j = rand(rng,1:10)
        kp[j] = rand(rng,Bool) ? 0.20*seed_kp[j] : 4.20*seed_kp[j]
    end
    rho,kp,ki
end

function main()
    seed = TOML.parsefile(joinpath(@__DIR__,"seed_uniform_875.toml"))
    rho0,kp0,ki0 = Float64.(seed["rho"]),Float64.(seed["Kp"]),Float64.(seed["Ki"])
    ctx = Q56.N.design_context(ROOT)
    m0 = Q56.N.descriptor(ctx,rho0,kp0,ki0)
    patterns = CSV.read(joinpath(ROOT,"experiments","delay_dressed_replacement_frontier_20261002","FROZEN_DELAY_PATTERNS.csv"),DataFrame)
    ids = ["min_commutator","max_commutator","low_GSP_frequency"]
    tau_patterns = Dict(id=>parse.(Float64,split(only(patterns.tau_ms[patterns.delay_pattern_id.==id]),";")) ./1000 for id in ids)
    s = -0.05+2pi*im*0.5
    references = Dict(id=>Q56.descriptor_dde_pencil(ctx,rho0,kp0,ki0,tau_patterns[id],s,m0) for id in ids)
    rng = MersenneTwister(20261002)
    recrows,detrows,rankrows = NamedTuple[],NamedTuple[],NamedTuple[]
    for case_number in 0:40
        sample_id = case_number == 0 ? "baseline" : "sample_$(lpad(case_number,2,'0'))"
        kind = case_number == 0 ? "REPRODUCED_FEASIBLE_SEED" :
            (case_number <=20 ? "WITHIN_BOUNDS_EVENTS_UNVALIDATED" : "GAIN_BOUND_VIOLATION")
        rho,kp,ki = case_number == 0 ? (rho0,kp0,ki0) : draw_design(rng,rho0,kp0,ki0,kind)
        maxerr,maxlog,maxphase,maxrank = 0.0,0.0,0.0,0
        status = "COMPLETED"
        try
            m1 = Q56.N.descriptor(ctx,rho,kp,ki)
            for id in ids
                tau = tau_patterns[id]
                p0 = references[id]
                p1 = Q56.descriptor_dde_pencil(ctx,rho,kp,ki,tau,s,m1)
                r = factor_difference(p0,p1,kp0,ki0,kp,ki,tau,s)
                push!(recrows,(;sample_id,kind,delay_pattern_id=id,full_descriptor_dimension=r.dimension,
                    action_rank_upper_bound=30,numerical_difference_rank=r.numerical_rank,
                    reconstruction_relative_error=r.err))
                push!(detrows,(;sample_id,kind,delay_pattern_id=id,
                    determinant_logabs_error=r.logabs_error,determinant_phase_error_rad=r.phase_error,
                    reference_pencil_nonsingular=true))
                maxerr = max(maxerr,r.err)
                maxlog = max(maxlog,r.logabs_error)
                maxphase = max(maxphase,r.phase_error)
                maxrank = max(maxrank,r.numerical_rank)
            end
        catch err
            status = "MODEL_OR_FACTOR_FAILURE: "*sprint(showerror,err)
        end
        push!(rankrows,(;sample_id,kind,status,max_numerical_difference_rank=maxrank,
            max_reconstruction_relative_error=maxerr,max_determinant_logabs_error=maxlog,
            max_determinant_phase_error_rad=maxphase,fully_event_feasible=case_number==0))
        CSV.write(joinpath(@__DIR__,"M2_ACTION_SPACE_RANK.csv"),DataFrame(rankrows))
        CSV.write(joinpath(@__DIR__,"M2_RECONSTRUCTION_VALIDATION.csv"),DataFrame(recrows))
        CSV.write(joinpath(@__DIR__,"M2_DETERMINANT_IDENTITY.csv"),DataFrame(detrows))
        println("M2_CASE ",sample_id," ",status," maxerr=",maxerr);flush(stdout)
    end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
