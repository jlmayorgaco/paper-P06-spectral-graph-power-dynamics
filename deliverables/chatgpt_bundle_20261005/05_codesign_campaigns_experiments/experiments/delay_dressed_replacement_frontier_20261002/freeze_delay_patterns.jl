using LinearAlgebra,CSV,DataFrames,ForwardDiff,Random,SHA,Statistics
include(joinpath(@__DIR__,"..","nonlinear_codesign_20261001","ReducedDAE.jl"))
include(joinpath(@__DIR__,"..","..","src","pd39","PD39.jl"))
const R=ReducedDAE;const ROOT=R.ROOT;const OUT=@__DIR__
BLAS.set_num_threads(1)
ctx=R.N.design_context(ROOT);data=PD39.ieee39_data();nb=39
V0=ctx.net.voltage;theta0=angle.(V0);vm=abs.(V0)

# Exact full lossy active-power angle Jacobian for the frozen passive Ybus.
Yc=Matrix{ComplexF64}(undef,nb,nb)
for i in 1:nb,j in 1:nb
    Yc[i,j]=complex(ctx.net.y_static[2i-1,2j-1],ctx.net.y_static[2i,2j-1])
end
Ptheta=theta->begin
    v=vm.*complex.(cos.(theta),sin.(theta))
    real.(v.*conj.(Yc*v))
end
Jfull=ForwardDiff.jacobian(Ptheta,theta0)
ports=collect(30:39);interior=setdiff(1:nb,ports)
Jport=Jfull[ports,ports]-Jfull[ports,interior]*(Jfull[interior,interior]\Jfull[interior,ports])

# Separate valid lossless physical graph: gamma_e=b_e*V_i*V_j*cos(eta_e).
branches=data.branch;inc=zeros(Float64,nb,nrow(branches));weights=zeros(Float64,nrow(branches))
for (e,r) in enumerate(eachrow(branches))
    i=Int(r.src_bus);j=Int(r.dst_bus);tap=Float64(r.r_src)
    y=inv(complex(Float64(r.R),Float64(r.X)))
    bij=imag(-tap*y);eta=theta0[i]-theta0[j]
    inc[i,e]=1;inc[j,e]=-1
    weights[e]=bij*vm[i]*vm[j]*cos(eta)
end
LPfull=inc*Diagonal(weights)*inc'
Lport=LPfull[ports,ports]-LPfull[ports,interior]*(LPfull[interior,interior]\LPfull[interior,ports])
Lport=(Lport+Lport')/2

# Reference inertia is the original initialized SG inertia on buses 30:39,
# fixed before any replacement/delay result, to define a reproducible graph basis.
machines=data.machine;M=Float64[]
for b in ports
    r=only(eachrow(machines[machines.bus.==b,:]))
    push!(M,2Float64(r.H)*Float64(r.Sn)/100)
end
Mmat=Diagonal(M);Minvhalf=Diagonal(1 ./sqrt.(M));Lbar=Minvhalf*Lport*Minvhalf
eig=eigen(Symmetric(Lbar));lambda=eig.values;U=eig.vectors

symres=norm(Lport-Lport',Inf)/max(norm(Lport,Inf),eps())
rowsum=norm(Lport*ones(10),Inf)/max(norm(Lport,Inf),eps())
mineig=minimum(eigvals(Symmetric(Lport)))
minoff=minimum([Lport[i,j] for i in 1:10 for j in 1:10 if i!=j])
is_laplacian=symres<1e-10 && rowsum<1e-8 && mineig>=-1e-8 && minoff<=0 && maximum([Lport[i,j] for i in 1:10 for j in 1:10 if i!=j])<=1e-8

function scores(t)
    τ=Float64.(t)./1000
    T=Diagonal(τ);chi=norm(Lbar*T-T*Lbar)/max(norm(Lbar),eps())
    rough=dot(τ,Lport*τ);hat=U'*τ
    low=sum(max(lambda[k],0)*abs2(hat[k]) for k in eachindex(lambda))
    high=sum(max(lambda[k],0)^2*abs2(hat[k]) for k in eachindex(lambda))
    (;chi,rough,low,high)
end
function digest(t)
    bytes2hex(sha256(codeunits(join(repr.(Float64.(t)),","))))
end

base=collect(0:9).*(50/9);rng=MersenneTwister(20261002)
randoms=[randperm(rng,10) |> p->base[p] for _ in 1:100]
function search_permutations(base,rng,nsearch)
    bestchi=(Inf,nothing);worstchi=(-Inf,nothing);bestlow=(Inf,nothing);worsthigh=(-Inf,nothing)
    for _ in 1:nsearch
        t=base[randperm(rng,10)];q=scores(t)
        q.chi<bestchi[1] && (bestchi=(q.chi,t))
        q.chi>worstchi[1] && (worstchi=(q.chi,t))
        q.low<bestlow[1] && (bestlow=(q.low,t))
        q.high>worsthigh[1] && (worsthigh=(q.high,t))
    end
    (;bestchi,worstchi,bestlow,worsthigh)
end
found=search_permutations(base,rng,100_000)
bestchi=found.bestchi;worstchi=found.worstchi;bestlow=found.bestlow;worsthigh=found.worsthigh
selected=[("min_commutator",bestchi[2]),("max_commutator",worstchi[2]),
          ("low_GSP_frequency",bestlow[2]),("high_GSP_frequency",worsthigh[2])]
allpatterns=vcat(selected,[("random_$(lpad(i,3,'0'))",randoms[i]) for i in eachindex(randoms)])
patterns=NamedTuple[];metrics=NamedTuple[]
for (id,t) in allpatterns
    q=scores(t)
    push!(patterns,(;delay_pattern_id=id,tau_ms=join(Float64.(t),";"),mean_tau_ms=mean(t),std_tau_ms=std(t),max_tau_ms=maximum(t),
        vector_sha256=digest(t),multiset_sha256=digest(sort(Float64.(t)))))
    push!(metrics,(;delay_pattern_id=id,chi_tau=q.chi,graph_roughness=q.rough,low_frequency_score=q.low,
        high_frequency_score=q.high,mean_tau_ms=mean(t),std_tau_ms=std(t),max_tau_ms=maximum(t),
        laplacian_valid=is_laplacian,operator_name=is_laplacian ? "L_P_lossless_Kron_Laplacian" : "lossless_physical_graph_operator"))
end
CSV.write(joinpath(OUT,"FROZEN_DELAY_PATTERNS.csv"),DataFrame(patterns))
CSV.write(joinpath(OUT,"TABLE_D05_GSP_DELAY_METRICS.csv"),DataFrame(metrics))
CSV.write(joinpath(OUT,"GRAPH_LP_LOSSLESS_PORT.csv"),DataFrame(Lport,:auto))
CSV.write(joinpath(OUT,"GRAPH_LP_EXACT_ANGLE_JACOBIAN_PORT.csv"),DataFrame(Jport,:auto))
CSV.write(joinpath(OUT,"GRAPH_MASS_REFERENCE.csv"),DataFrame(bus=ports,M=M))
CSV.write(joinpath(OUT,"GRAPH_MODES.csv"),DataFrame(mode=1:10,eigenvalue=lambda))
CSV.write(joinpath(OUT,"GRAPH_GSP_BASIS.csv"),DataFrame(U,:auto))
summary=(;Lport_symmetry_residual=symres,Lport_row_sum_relative=rowsum,Lport_min_eigenvalue=mineig,
    Lport_max_offdiagonal=maximum([Lport[i,j] for i in 1:10 for j in 1:10 if i!=j]),
    Lport_is_laplacian=is_laplacian,exact_angle_jacobian_symmetry_residual=norm(Jport-Jport',Inf)/max(norm(Jport,Inf),eps()),
    exact_angle_vs_lossless_relative=norm(Jport-Lport)/max(norm(Jport),eps()),mass_definition="original IEEE39 initialized SG: 2*H*Sn/Sbase, frozen reference",
    seed=20261002,heuristic_random_candidates=100000,random_permutation_count=100)
open(joinpath(OUT,"GRAPH_OPERATOR_AUDIT.toml"),"w") do io
    for (k,v) in pairs(summary);println(io,k," = ",v isa Bool ? string(v) : v isa AbstractString ? "\"$v\"" : v);end
end
println(summary)
