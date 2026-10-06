include("solve.jl")
main()
# The final checkpoint remains explicitly exploratory. Audit in the same
# process so compilation is shared; no old candidate or freeze is overwritten.
ARGS[1]=joinpath(L.OUT,"FINAL_EXPLORATORY_CANDIDATE.toml")
include("second_order_audit.jl")
