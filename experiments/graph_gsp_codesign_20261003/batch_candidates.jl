include(joinpath(@__DIR__,"DelayedEvents.jl"))
include(joinpath(@__DIR__,"contour.jl"))
using .DelayedEvents, CSV, DataFrames
function batch(shard,nshards)
    ctx=DelayedEvents.R.N.design_context(DelayedEvents.R.ROOT)
    families=sort(String.(CSV.read(joinpath(@__DIR__,"TABLE_07_LOCAL_PREDICTORS.csv"),DataFrame).family))
    checkpoint=joinpath(@__DIR__,"BATCH_SHARD_$(shard).csv")
    summary=isfile(checkpoint) ? [NamedTuple(r) for r in eachrow(CSV.read(checkpoint,DataFrame))] : NamedTuple[]
    for (j,family) in enumerate(families)
        mod(j-1,nshards)+1==shard || continue
        any(r.family==family && r.fully_feasible for r in summary) && continue
        for fraction in ("1","0.5","0.25","0.125")
            id=family*"_f"*fraction;path=joinpath(@__DIR__,"proposals",id*".toml")
            any(r.design_id==id for r in summary) && continue
            println("CANDIDATE_START ",id);flush(stdout)
            spectralpath=joinpath(@__DIR__,"spectral",id*"_contour.csv")
            result=isfile(spectralpath) ? NamedTuple(only(eachrow(CSV.read(spectralpath,DataFrame)))) : check_contour(path;ctx,R=DelayedEvents.R)
            if result.status=="PASS_NUMERICAL"
                events=DelayedEvents.run_cases("eval",path;ctx)
                eventpass=all(x.pass for x in events)
                push!(summary,(;family,design_id=id,fraction=parse(Float64,fraction),spectrum=result.status,
                    event_pass=eventpass,fully_feasible=eventpass))
            else
                push!(summary,(;family,design_id=id,fraction=parse(Float64,fraction),spectrum=result.status,
                    event_pass=false,fully_feasible=false))
            end
            CSV.write(joinpath(@__DIR__,"BATCH_SHARD_$(shard).csv"),DataFrame(summary))
            last(summary).fully_feasible && break
        end
    end
end
batch(parse(Int,ARGS[1]),parse(Int,ARGS[2]))
