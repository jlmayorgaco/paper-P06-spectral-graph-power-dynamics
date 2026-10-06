using SHA, TOML, CSV, DataFrames, LinearAlgebra
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","experiment_N")
const CANDIDATE=joinpath(OUT,"Z_N_NOMINAL_FINAL.toml")
bytes2hex(sha256(read(CANDIDATE)))==strip(read(CANDIDATE*".sha256",String)) ||
    error("candidate SHA mismatch")
include(joinpath(ROOT,"src","bnd_model_expN","PDExactDesignN.jl"))
using .PDExactDesignN

function main()
    z=TOML.parsefile(CANDIDATE)
    rho=Float64.(z["rho"]);kp=Float64.(z["Kp"]);ki=Float64.(z["Ki"])
    ctx=design_context(ROOT)
    sp=spectrum(ctx,rho,kp,ki)
    d=simple_mode_sensitivities(ctx,rho,kp,ki)
    λ=d.lambda
    m=sp.model;nport=size(ctx.net.y_static,1)
    T=ComplexF64.(ctx.net.y_static)
    Ts=zeros(ComplexF64,nport,nport)
    bybus=Dict{Int,Vector{Any}}()
    for b in m.blocks
        push!(get!(bybus,b.bus,Any[]),b)
        yi=(2b.bus-1):(2b.bus)
        R=(λ*I(size(b.J.A,1))-b.J.A)\Matrix{ComplexF64}(I,size(b.J.A,1),size(b.J.A,1))
        Y=b.J.D+b.J.C*R*b.J.B
        Ys=-b.J.C*R*R*b.J.B
        T[yi,yi].+=b.share.*Y
        Ts[yi,yi].+=b.share.*Ys
    end
    F=svd(T);v=F.V[:,end]
    pivot=argmax([norm(v[(2i-1):(2i)]) for i in 1:39])
    kk=collect((2pivot-1):(2pivot));rr=setdiff(collect(1:nport),kk)
    R=inv(T[rr,rr])
    Ckr=T[kk,rr];Crk=T[rr,kk]
    Gamma=-Ckr*R*Crk
    Teff=T[kk,kk]+Gamma
    FT=svd(Teff);u=FT.U[:,end];w=FT.V[:,end]
    function split(dT)
        direct=dT[kk,kk]
        self=-dT[kk,rr]*R*Crk-Ckr*R*dT[rr,kk]+
             Ckr*R*dT[rr,rr]*R*Crk
        direct,self
    end
    ds_direct,ds_self=split(Ts)
    denominator=dot(u,(ds_direct+ds_self)*w)
    paths=NamedTuple[]
    for l in 1:39, n in 1:39
        l==pivot && continue;n==pivot && continue
        il=collect((2l-1):(2l));im=collect((2n-1):(2n))
        rl=[findfirst(==(x),rr) for x in il]
        rm=[findfirst(==(x),rr) for x in im]
        G=-T[kk,il]*R[rl,rm]*T[im,kk]
        push!(paths,(;pivot_bus=pivot,from_bus=l,to_bus=n,
            pathway_frobenius=norm(G),trace_real=real(tr(G)),
            trace_imag=imag(tr(G))))
    end
    cancellation=sum(p.pathway_frobenius for p in paths)/max(norm(Gamma),eps())
    rows=NamedTuple[]
    params=[("rho",38)]
    append!(params,[("Kp",i) for i in 30:39])
    append!(params,[("Ki",i) for i in 30:39])
    for (kind,bus) in params
        yi=(2bus-1):(2bus)
        dlocal=zeros(ComplexF64,2,2)
        exact=0.0+0.0im
        if kind=="rho"
            bs=bybus[bus]
            sg=only(filter(x->x.kind=="SG",bs))
            gf=only(filter(x->x.kind=="GFL",bs))
            dlocal=port_transfer(gf.J,λ)-port_transfer(sg.J,λ)
            exact=d.rho[bus-29]
        else
            gf=only(filter(x->x.kind=="GFL",bybus[bus]))
            op=trim_gfl(ctx,bus)
            J=gf.J
            Rj=(λ*I(size(J.A,1))-J.A)\Matrix{ComplexF64}(I,size(J.A,1),size(J.A,1))
            dA=zeros(ComplexF64,size(J.A));dB=zeros(ComplexF64,size(J.B))
            θ=op.x[3];ur=op.u[1];ui=op.u[2]
            scale=kind=="Kp" ? 1/op.pars.pll_tau : 1.0
            row=kind=="Kp" ? 4 : 5
            dA[row,3]=scale*(-cos(θ)*ur-sin(θ)*ui)
            dB[row,1]=scale*(-sin(θ))
            dB[row,2]=scale*cos(θ)
            dlocal=gf.share.*(J.C*Rj*dA*Rj*J.B+J.C*Rj*dB)
            exact=kind=="Kp" ? d.Kp[bus-29] : d.Ki[bus-29]
        end
        dT=zeros(ComplexF64,nport,nport)
        dT[yi,yi].=dlocal
        direct,self=split(dT)
        dl_direct=-dot(u,direct*w)/denominator
        dl_self=-dot(u,self*w)/denominator
        total=dl_direct+dl_self
        push!(rows,(;pole_real=real(λ),pole_imag=imag(λ),
            pivot_bus=pivot,parameter=kind,parameter_bus=bus,
            direct_dlambda_real=real(dl_direct),
            self_energy_dlambda_real=real(dl_self),
            total_dlambda_real=real(total),
            exact_dlambda_real=real(exact),
            derivative_error=abs(total-exact),
            self_energy_share=abs(dl_self)/max(abs(dl_direct)+abs(dl_self),eps()),
            cancellation_ratio=cancellation,
            port_schur_sigma_min=FT.S[end],
            full_port_sigma_min=F.S[end]))
    end
    CSV.write(joinpath(OUT,"TABLE_N15_self_energy_explanation.csv"),DataFrame(rows))
    CSV.write(joinpath(OUT,"TABLE_N15_pairwise_pathways.csv"),DataFrame(paths))
    # The frozen closure stores complex Y_port as 2×2 real rectangular blocks.
    # Reconstruct the physical 39×39 complex graph operator before taking
    # real/imaginary parts; imag(y_static) would be identically zero.
    Yport=zeros(ComplexF64,39,39)
    rectangular_error=0.0
    for i in 1:39,j in 1:39
        block=ctx.net.y_static[(2i-1):(2i),(2j-1):(2j)]
        g=(block[1,1]+block[2,2])/2
        b=(block[2,1]-block[1,2])/2
        Yport[i,j]=g+im*b
        rectangular_error=max(rectangular_error,
            norm(block-[g -b;b g]))
    end
    LG=real.(Yport);LB=imag.(Yport)
    comm=LG*LB-LB*LG
    strength=[sum(abs(Yport[i,j]) for j in 1:39 if j!=i) for i in 1:39]
    rank38=1+count(x->x>strength[38],strength)
    CSV.write(joinpath(OUT,"TABLE_N15_graph_diagnostic.csv"),DataFrame([(
        commutator_relative_norm=norm(comm)/max(norm(LG)*norm(LB),eps()),
        rectangular_block_max_error=rectangular_error,
        bus38_graph_strength=strength[38],
        bus38_strength_rank_descending=rank38,
        gain_pattern="constant_at_box_corners",
        graph_strength_gain_correlation="undefined_constant_gain",
        interpretation="graph operators used only as diagnostics")]))
    println("BND pivot=",pivot," pole=",λ," gamma_norm=",norm(Gamma),
        " cancellation=",cancellation," max_derivative_error=",
        maximum(r.derivative_error for r in rows));flush(stdout)
end

main()
