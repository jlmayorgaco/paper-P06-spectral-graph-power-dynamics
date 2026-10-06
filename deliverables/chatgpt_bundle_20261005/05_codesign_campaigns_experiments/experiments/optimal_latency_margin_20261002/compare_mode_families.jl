using CSV, DataFrames, LinearAlgebra, TOML

module MegaOracle
include(joinpath(@__DIR__,"..","analytical_delay_codesign_mega_20261002","m3_a_trace_integral.jl"))
end

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const MEGA=joinpath(ROOT,"experiments","analytical_delay_codesign_mega_20261002")
BLAS.set_num_threads(1)

function mode(ctx,path,tau_ms,s)
    d=TOML.parsefile(path)
    L=MegaOracle.DC.linearization(ctx,Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"]))
    tau=fill(tau_ms/1000,length(L.Ai))
    D=MegaOracle.DC.delta_matrix(L,s,tau)
    x=svd(D).V[:,end]
    q=L.Q*x;q/=norm(q)
    port=[abs(dot(L.C[:,i],x))*norm(L.B[:,i]) for i in 1:10]
    port/=sum(port)
    (;q,port,s)
end

function main()
    ctx=MegaOracle.DC.R.N.design_context(ROOT)
    two=CSV.read(joinpath(@__DIR__,"T02_PRECISE_CROSSINGS.csv"),DataFrame)
    fixed=CSV.read(joinpath(@__DIR__,"T04_PRECISE_CROSSINGS.csv"),DataFrame)
    a=only(eachrow(two[two.design_id.=="seed_875",:]))
    b=only(eachrow(two[two.design_id.=="best_zero_delay_88455",:]))
    c=only(eachrow(fixed[fixed.design_id.=="rho_04",:]))
    modes=Dict(
        "seed_nominal"=>mode(ctx,joinpath(MEGA,"seed_uniform_875.toml"),
            a.local_root_crossing_ms,a.critical_root_real+im*a.critical_root_imag),
        "high_rho_nominal"=>mode(ctx,joinpath(@__DIR__,"fixed_gain_designs","rho_04.toml"),
            c.fixed_gain_tau_crit_ms,c.critical_root_real+2pi*im*c.critical_frequency_hz),
        "high_rho_zero_delay_tuned"=>mode(ctx,joinpath(MEGA,"M1_ZERO_DELAY_DESIGN.toml"),
            b.local_root_crossing_ms,b.critical_root_real+im*b.critical_root_imag))
    rows=NamedTuple[]
    for (left,right) in (("seed_nominal","high_rho_nominal"),
                          ("seed_nominal","high_rho_zero_delay_tuned"),
                          ("high_rho_nominal","high_rho_zero_delay_tuned"))
        x,y=modes[left],modes[right]
        mac=abs(dot(x.q,y.q))^2
        port_overlap=sum(sqrt.(x.port.*y.port))
        push!(rows,(;design_left=left,design_right=right,
            full_physical_right_vector_MAC=mac,
            PLL_port_score_Bhattacharyya_overlap=port_overlap,
            frequency_left_hz=imag(x.s)/(2pi),frequency_right_hz=imag(y.s)/(2pi),
            status="NUMERICAL_MODE_SHAPE_COMPARISON_NOT_BIORTHOGONAL_PARTICIPATION"))
    end
    CSV.write(joinpath(@__DIR__,"T03_FAMILY_OVERLAP.csv"),DataFrame(rows))
    foreach(x->println("FAMILY_OVERLAP ",x),rows)
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
