ENV["BND_SEARCH_SUBDIRECTORY"]="exchange_checkpoint"
include("search_oracle.jl")
input=CSV.read(joinpath(DEST,"parameters.csv"),DataFrame)
rows=NamedTuple[]
for row in eachrow(input)
    p=Float64[row[Symbol("p",i)] for i in 1:30]
    c=evaluate(p)
    push!(rows,(;id=row.id,increase_bus=row.increase_bus,decrease_bus=row.decrease_bus,
        MW=row.MW,tuned=row.tuned,feasible=maximum(c)<=0,
        constraints=join(c,";")))
    CSV.write(joinpath(DEST,"validation.csv"),DataFrame(rows))
    println("EXCHANGE ",row.id," max_constraint=",maximum(c));flush(stdout)
end
