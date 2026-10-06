using CSV, DataFrames, LinearAlgebra, TOML

module MegaOracle
include(joinpath(@__DIR__,"..","analytical_delay_codesign_mega_20261002","m3_a_trace_integral.jl"))
end

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const MEGA=joinpath(ROOT,"experiments","analytical_delay_codesign_mega_20261002")
BLAS.set_num_threads(1)

function main()
    table=CSV.read(joinpath(@__DIR__,"T01_DELAY_TRANSITION_REPRODUCTION.csv"),DataFrame)
    ctx=MegaOracle.DC.R.N.design_context(ROOT)
    rows=NamedTuple[]
    for id in unique(table.design_id)
        path=id=="seed_875" ? joinpath(MEGA,"seed_uniform_875.toml") :
            joinpath(MEGA,"M1_ZERO_DELAY_DESIGN.toml")
        d=TOML.parsefile(path)
        L=MegaOracle.DC.linearization(ctx,Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"]))
        for point in eachrow(table[(table.design_id.==id).&(table.status.=="SAFE"),:])
            gamma=Float64(point.rightmost_real_part_candidate)+0.001
            println("RIGHTMOST_EXCLUSION ",id," tau_ms=",point.tau_ms," gamma=",gamma);flush(stdout)
            try
                r=MegaOracle.trace_count(L,fill(point.tau_ms/1000,length(L.Ai));
                    gamma,step=20.0)
                passed=r.nearest==0 && abs(r.phase_estimate)<0.01 &&
                    abs(r.estimate)<0.01 && r.largest_phase_step<pi/2
                push!(rows,(;design_id=String(id),tau_ms=point.tau_ms,
                    candidate_rightmost_real=point.rightmost_real_part_candidate,
                    exclusion_line_real=gamma,
                    roots_right_of_exclusion=passed ? 0 : r.nearest,
                    trace_count_real=real(r.estimate),phase_count=r.phase_estimate,
                    status=passed ? "RIGHTMOST_WITHIN_0P001_S_INV_NUMERICALLY" : "INDETERMINATE",
                    contour_evaluations=r.evaluations,error=""))
            catch err
                push!(rows,(;design_id=String(id),tau_ms=point.tau_ms,
                    candidate_rightmost_real=point.rightmost_real_part_candidate,
                    exclusion_line_real=gamma,roots_right_of_exclusion=missing,
                    trace_count_real=NaN,phase_count=NaN,status="INDETERMINATE",
                    contour_evaluations=missing,error=sprint(showerror,err)))
            end
            CSV.write(joinpath(@__DIR__,"T01_RIGHTMOST_EXCLUSION.csv"),DataFrame(rows))
            println("EXCLUSION_RESULT ",last(rows));flush(stdout)
        end
    end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
