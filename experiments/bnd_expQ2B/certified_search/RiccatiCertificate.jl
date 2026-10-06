function riccati_certificate(A;beta=L.BETA,outdir=L.OUT)
    t=time();As=Matrix(A+.05I);n=size(A,1);bsolve=beta*(1+2e-8)
    H=[As bsolve*I;-bsolve*I -As'];F=eigen(H);stable=findall(real.(F.values).<0)
    length(stable)==n || error("Hamiltonian stable subspace count")
    U=F.vectors[1:n,stable];V=F.vectors[n+1:2n,stable]
    X=real.(V/U);Xb=BigFloat.((X+X')/2);Ab=BigFloat.(As);bb=BigFloat(bsolve)
    history=NamedTuple[]
    for k in 1:16
        Eb=Ab'*Xb+Xb*Ab+bb*(I+Xb*Xb)
        push!(history,(iteration=k,residual=Float64(norm(Eb))))
        Float64(norm(Eb))<1e-22 && break
        Ac=Float64.(Ab+bb*Xb)
        dx=lyap(Ac',Float64.(Eb));Xb+=BigFloat.((dx+dx')/2)
    end
    # Strict bounded-real LMI: positive X and negative residual prove the
    # resolvent gain below 1/beta for ALL real frequencies, including infinity.
    W=-(Ab'*Xb+Xb*Ab+BigFloat(beta)*(I+Xb*Xb))
    positivity=isposdef(LinearAlgebra.Symmetric(Xb));negativity=isposdef(LinearAlgebra.Symmetric(W))
    eigW=eigmin(LinearAlgebra.Symmetric(Float64.(W)));eigX=eigmin(LinearAlgebra.Symmetric(Float64.(Xb)))
    residual=Float64(norm(Ab'*Xb+Xb*Ab+bb*(I+Xb*Xb)))
    CSV.write(joinpath(outdir,"RICCATI_NEWTON_HISTORY.csv"),DataFrame(history))
    CSV.write(joinpath(outdir,"RICCATI_STORAGE_X_256BIT.csv"),DataFrame(string.(Xb),:auto))
    out=Dict("status"=>positivity&&negativity ? "NUMERICALLY_CERTIFIED_FULL_BAND_BOUNDED_REAL" : "INCONCLUSIVE",
        "beta_lower_tested"=>beta,"beta_Riccati_solve"=>bsolve,"X_positive"=>positivity,
        "strict_negative_LMI"=>negativity,"X_min_eigenvalue_Float64"=>eigX,
        "negative_LMI_min_eigenvalue_Float64"=>eigW,"Riccati_residual_big"=>residual,
        "precision_bits"=>precision(BigFloat),"rounding"=>"nearest; NOT outward-rounded",
        "formal"=>false,"Hamiltonian_min_abs_real"=>minimum(abs,real.(F.values)),
        "invariant_subspace_condition"=>cond(U),"elapsed_s"=>time()-t,
        "coverage"=>"Strict bounded-real storage LMI, entire real-frequency axis and infinite-frequency tail")
    open(joinpath(outdir,"RICCATI_FULL_BAND_CERTIFICATE.toml"),"w") do io;TOML.print(io,out);end
    println(out);flush(stdout)
    out
end
