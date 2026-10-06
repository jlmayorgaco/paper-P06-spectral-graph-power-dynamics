using CSV, DataFrames, TOML, LinearAlgebra
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const NEV=include(joinpath(ROOT,"experiments","delay_dressed_replacement_frontier_20261002","DDE_NEV.jl"))
const D=NEV.D;const R=D.R
function main()
    d=TOML.parsefile(joinpath(@__DIR__,"seed_uniform_875.toml"));ctx=R.N.design_context(ROOT)
    L=D.linearization(ctx,Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"]))
    rows=NamedTuple[];frequencies=[0.0,50.0,200.0,500.0,1000.0,2000.0,5000.0,9000.0]
    for tau_ms in (20.0,40.0)
        tau=fill(tau_ms/1000,length(L.Ai))
        for re0 in (-0.02,0.05),omega0 in frequencies
            s0=re0+im*omega0;M=D.delta_matrix(L,s0,tau);F=svd(M)
            root=NEV.newton_root(L,s0,tau,ComplexF64.(F.V[:,end]);tol=1e-9,maxiter=30,max_step=10.0)
            push!(rows,(;tau_ms,seed_real=re0,seed_imag=omega0,root_real=real(root.s),root_imag=imag(root.s),
                frequency_Hz=abs(imag(root.s))/(2pi),residual=root.residual,converged=root.converged,
                roots_right_of_margin=root.converged && real(root.s)>-0.05,
                solver="bordered_Newton_full_203_state_NEp",claim="local seeded root only"))
        end
        CSV.write(joinpath(@__DIR__,"TABLE_Q03_SEEDED_NEP_ROOTS.csv"),DataFrame(rows))
        println("SEEDED_ROOTS_DONE tau_ms=",tau_ms," converged=",count(r->r.tau_ms==tau_ms&&r.converged,rows),
            " right_of_margin=",count(r->r.tau_ms==tau_ms&&r.roots_right_of_margin,rows));flush(stdout)
    end
end
main()
