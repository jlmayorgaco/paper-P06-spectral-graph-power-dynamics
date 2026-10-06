using CSV, DataFrames, LinearAlgebra, TOML

module MegaOracle
include(joinpath(@__DIR__,"..","analytical_delay_codesign_mega_20261002","m3_a_trace_integral.jl"))
end

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const MEGA=joinpath(ROOT,"experiments","analytical_delay_codesign_mega_20261002")
BLAS.set_num_threads(1)

function main()
    crossings=CSV.read(joinpath(@__DIR__,"T02_PRECISE_CROSSINGS.csv"),DataFrame)
    ctx=MegaOracle.DC.R.N.design_context(ROOT)
    rows=NamedTuple[]
    for item in eachrow(crossings)
        design_id=String(item.design_id)
        path=design_id=="seed_875" ? joinpath(MEGA,"seed_uniform_875.toml") :
            joinpath(MEGA,"M1_ZERO_DELAY_DESIGN.toml")
        d=TOML.parsefile(path)
        L=MegaOracle.DC.linearization(ctx,Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"]))
        for side in ("below","above")
            tau_ms=item.local_root_crossing_ms+(side=="below" ? -0.01 : 0.01)
            println("BOUNDARY_VERIFY ",design_id," ",side," ",tau_ms);flush(stdout)
            try
                r=MegaOracle.trace_count(L,fill(tau_ms/1000,length(L.Ai));
                    gamma=-0.05,step=10.0,local_tolerance=1e-7,max_evaluations=300000,max_depth=20)
                consistent=abs(real(r.estimate)-r.nearest)<0.005 &&
                    abs(r.phase_estimate-r.nearest)<0.005 && r.largest_phase_step<pi/2
                status=consistent ? (r.nearest==0 ? "SAFE" : "UNSAFE") : "INDETERMINATE"
                push!(rows,(;design_id,side,tau_ms,status,
                    roots_violating_margin=consistent ? r.nearest : missing,
                    trace_count_real=real(r.estimate),phase_count=r.phase_estimate,
                    max_phase_step=r.largest_phase_step,min_sigma_small=r.min_sigma_small,
                    quadrature_error_estimate=r.quadrature_error_estimate,
                    contour_evaluations=r.evaluations,error=""))
            catch err
                push!(rows,(;design_id,side,tau_ms,status="INDETERMINATE",
                    roots_violating_margin=missing,trace_count_real=NaN,phase_count=NaN,
                    max_phase_step=NaN,min_sigma_small=NaN,quadrature_error_estimate=NaN,
                    contour_evaluations=missing,error=sprint(showerror,err)))
            end
            CSV.write(joinpath(@__DIR__,"T02_BOUNDARY_VERIFICATION.csv"),DataFrame(rows))
            println("BOUNDARY_RESULT ",last(rows));flush(stdout)
        end
    end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
