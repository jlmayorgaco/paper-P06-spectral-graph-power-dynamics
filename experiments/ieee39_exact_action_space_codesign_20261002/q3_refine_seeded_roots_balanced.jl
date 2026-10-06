using CSV,DataFrames,TOML,LinearAlgebra
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const NEV=include(joinpath(ROOT,"experiments","delay_dressed_replacement_frontier_20261002","DDE_NEV.jl"))
const D=NEV.D;const R=D.R
function balance_diagonal(A)
    d=ones(Float64,size(A,1));n=length(d)
    for _ in 1:100
        changed=false
        for i in 1:n
            r=0.0;c=0.0
            for j in 1:n
                i==j && continue
                r+=abs(A[i,j])*d[j]/d[i];c+=abs(A[j,i])*d[i]/d[j]
            end
            (r==0||c==0)&&continue
            f=2.0^round(Int,0.5*log2(r/c))
            if f!=1;d[i]*=f;changed=true;end
        end
        !changed&&break
    end
    d
end
function main()
    d=TOML.parsefile(joinpath(@__DIR__,"seed_uniform_875.toml"));ctx=R.N.design_context(ROOT)
    L=D.linearization(ctx,Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"]))
    C=transpose(L.C);Ac=L.A0+L.B*C;scale=balance_diagonal(Ac)
    A0b=(L.A0.*transpose(scale))./reshape(scale,:,1);Bb=L.B./reshape(scale,:,1);Cb=C.*transpose(scale)
    Aib=[Bb[:,j]*transpose(Cb[j,:]) for j in eachindex(L.Ai)];Ab=A0b+sum(Aib;init=zeros(size(A0b)))
    Lb=(;A=Ab,A0=A0b,Ai=Aib);rows=NamedTuple[]
    for tau_ms in (20.0,40.0),s0 in (-0.02+0im,0.05+0im)
        tau=fill(tau_ms/1000,length(Aib));M=D.delta_matrix(Lb,s0,tau);sv=svd(M)
        root=NEV.newton_root(Lb,s0,tau,ComplexF64.(sv.V[:,end]);tol=1e-12,maxiter=40,max_step=0.5)
        Mr=D.delta_matrix(Lb,root.s,tau);sr=minimum(svdvals(Mr));nr=norm(Mr*root.v)/max(1,norm(Mr)*norm(root.v))
        push!(rows,(;tau_ms,seed=string(s0),real=real(root.s),imag=imag(root.s),
            relative_sigma_min=sr/max(opnorm(Mr),eps()),normalized_residual=nr,
            solver_residual=root.residual,converged=root.converged,iterations=root.iterations,
            right_of_margin=real(root.s)>-0.05,diagonal_scale_ratio=maximum(scale)/minimum(scale)))
    end
    CSV.write(joinpath(@__DIR__,"TABLE_Q03_BALANCED_ROOT_REFINEMENT.csv"),DataFrame(rows));println(DataFrame(rows))
end
main()
