ENV["BND_MODAL_BUNDLE_SIZE"]="6"
ENV["BND_SEARCH_SUBDIRECTORY"]="modal_bundle_audit"
include("search_oracle.jl")
d=TOML.parsefile(abspath(ARGS[1]));rho=Float64.(d["rho"]);kp=Float64.(d["Kp"]);ki=Float64.(d["Ki"])
modes=modal_bundle(rho,kp,ki)
rows=[(;rank=i,real=real(m.lambda),imag=imag(m.lambda),condition=m.condition) for (i,m) in enumerate(modes)]
CSV.write(joinpath(DEST,"modes.csv"),DataFrame(rows))
CSV.write(joinpath(DEST,"jacobian.csv"),DataFrame(reduce(vcat,[transpose(real.(vcat(m.rho,m.Kp.*kp,m.Ki.*ki))/.05) for m in modes]),:auto))
show(stdout,"text/plain",DataFrame(rows));println()
