ENV["BND_SEARCH_SUBDIRECTORY"]="uniform_baseline"
include("search_oracle.jl")
rows=NamedTuple[]
if isfile(joinpath(DEST,"summary.csv"))
    append!(rows,NamedTuple.(eachrow(CSV.read(joinpath(DEST,"summary.csv"),DataFrame))))
    ids=[parse(Int,first(split(f,'_'))) for f in readdir(DEST) if occursin(r"^\d+_eval.csv$",f)]
    COUNT[]=isempty(ids) ? 0 : maximum(ids)
end
function check(rho)
    cached=findfirst(r->r.rho==rho,rows)
    cached!==nothing && return rows[cached].feasible
    p=vcat(fill(rho,10),fill(log(R.N.K0P),10),fill(log(R.N.K0I),10))
    alpha=R.N.spectrum(CTX,p[1:10],exp.(p[11:20]),exp.(p[21:30])).alpha
    if alpha > -CFG["nominal_decay_margin_s_inv"]
        c=Float64[(alpha+CFG["nominal_decay_margin_s_inv"])/CFG["nominal_decay_margin_s_inv"]]
        complete=false
    else
        c=evaluate(p);complete=all(c.<9.)
    end
    row=(;rho,alpha,feasible=maximum(c)<=0,complete_TDS=complete,max_constraint=maximum(c))
    push!(rows,row);CSV.write(joinpath(DEST,"summary.csv"),DataFrame(rows))
    println("UNIFORM_RESULT ",row);flush(stdout)
    row.feasible
end
grid=[.75,.8,.85,.875,.9,.925,.95,.975,.999]
feasible=[check(r) for r in grid]
idx=findlast(feasible)
if idx!==nothing && idx<length(grid)
    let lo=grid[idx],hi=grid[idx+1]
        for k in 1:8
            mid=(lo+hi)/2
            if check(mid);lo=mid;else;hi=mid;end
        end
    end
end
open(joinpath(DEST,"scope.txt"),"w") do io
    println(io,"Uniform replacement, nominal PLL gains, same physical model and constraints. Coarse scan followed by local transition refinement; no global scalar monotonicity or optimality theorem. Modal failures are rejected before time-domain simulation.")
end
