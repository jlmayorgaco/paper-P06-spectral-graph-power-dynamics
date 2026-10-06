include("LocalOracle.jl")
using .LocalOracle, LinearAlgebra, TOML, CSV, DataFrames
const L=LocalOracle
d=TOML.parsefile(joinpath(L.OUT,"FEASIBLE_NUMERICAL_REPRESENTATIVE.toml"))
a=L.architecture(d["support"],d);x=L.encode(a.support,d);v=L.matrices(a,x)
As=Matrix(v.A+.05I);n=size(As,1);b=L.BETA
H=[As b*I;-b*I -As'];F=eigen(H);stable=findall(real.(F.values).<0)
println("HAMILTONIAN ",length(stable)," minreal=",minimum(abs,real.(F.values)));flush(stdout)
if length(stable)==n
    U=F.vectors[1:n,stable];V=F.vectors[n+1:2n,stable];X=real.(V/U);X=(X+X')/2
    println("X positivity=",eigmin(Symmetric(X))," norm=",opnorm(X)," invariant_condition=",cond(U));flush(stdout)
    for k in 1:8
        E=As'*X+X*As+b*(I+X*X);Ac=As+b*X
        println("RICCATI ",k," residual=",norm(E)," maximum=",eigmax(Symmetric(E))," closedalpha=",maximum(real.(eigvals(Ac))));flush(stdout)
        # Standard Newton correction is a Lyapunov solve, not an optimizer.
        global X=X+lyap(Ac',E);global X=(X+X')/2
    end
    CSV.write(joinpath(L.OUT,"RICCATI_X.csv"),DataFrame(X,:auto))
end
