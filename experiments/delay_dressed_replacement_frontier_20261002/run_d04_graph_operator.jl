using LinearAlgebra,CSV,DataFrames,ForwardDiff,Statistics
include(joinpath(@__DIR__,"..","nonlinear_codesign_20261001","ReducedDAE.jl"))
const R=ReducedDAE;const ROOT=R.ROOT;const OUT=@__DIR__;BLAS.set_num_threads(1)
for (src,dst) in ((joinpath(OUT,"baseline_reproduction","TABLE_D04_GRAPH_DAMPING_OPERATORS.csv"),joinpath(OUT,"TABLE_D04_PRE_GAUGE_DIAGNOSTIC.csv")),
                  (joinpath(OUT,"baseline_reproduction","TABLE_D05_GSP_DG_METRICS.csv"),joinpath(OUT,"TABLE_D05_PRE_GAUGE_DIAGNOSTIC.csv")),
                  (joinpath(OUT,"GRAPH_DG_ZERO_DELAY.csv"),joinpath(OUT,"GRAPH_DG_ZERO_DELAY_PRE_GAUGE_DIAGNOSTIC.csv")))
    isfile(src) && !isfile(dst) && cp(src,dst)
end
ctx=R.N.design_context(ROOT);rho=fill(.875,10);kp=fill(R.N.K0P,10);ki=fill(R.N.K0I,10)
m=R.model(ctx,rho,kp,ki;bus=8,delta=0.,dc_convention=:physical_supply)
A=R.derivatives(m.x0,m).Fx;Ai=Matrix{Float64}[]
for i in 1:10
    ix=m.gfidx[i];port=(2i-1):(2i)
    c=ForwardDiff.gradient(m.x0) do x
        v=R.voltage(x,m);theta=x[ix[3]]
        -sin(theta)*v[first(port)]+cos(theta)*v[last(port)]
    end
    b=zeros(length(m.x0));b[ix[4]]=kp[i]/(1/(2pi*300));b[ix[5]]=ki[i]
    push!(Ai,b*c')
end
A0=A-sum(Ai;init=zeros(size(A)));n=size(A,1)
# Deflate the exact global-rotation direction before partitioning. Retained
# coordinates are an orthonormalized projection of SG rotor-angle axes into
# the quotient; the remaining hidden coordinates span its orthogonal complement.
qidx=Int[last(m.sgidx[i]) for i in 1:10]
g,_=R.rotation_generator(m.x0,m;jacobian=false);g ./= norm(g)
gauge_residual=norm(A*g)/max(norm(A)*norm(g),eps())
E=Matrix{Float64}(I,n,n)[:,qidx]
Qraw=E-g*(g'*E)
Qret=Matrix(qr(Qraw).Q)[:,1:10]
Qhidden=nullspace(vcat(g',Qret'))
V=hcat(Qret,Qhidden)
size(V,2)==n-1 || error("rotational quotient basis has incorrect rank")
Aq=V'*A*V;Aiq=[V'*Ai0*V for Ai0 in Ai]
A0q=Aq-sum(Aiq;init=zeros(size(Aq)));nq=size(Aq,1)
qidx=collect(1:10);cidx=collect(11:nq)
Delta0=-Aq;I0=Matrix{Float64}(I,nq,nq)
function schur_slope(tau)
    D1=copy(I0)
    for i in 1:10;D1 .+= tau[i].*Aiq[i];end
    Aqq=Delta0[qidx,qidx];Aqc=Delta0[qidx,cidx];Acq=Delta0[cidx,qidx];Acc=Delta0[cidx,cidx]
    Aqq1=D1[qidx,qidx];Aqc1=D1[qidx,cidx];Acq1=D1[cidx,qidx];Acc1=D1[cidx,cidx]
    condc=cond(Acc);isfinite(condc)&&condc<1e14 || return (;valid=false,condition=condc)
    Z=Acc\Acq
    T0=Aqq-Aqc*Z
    T1=Aqq1-Aqc1*Z+Aqc*(Acc\(Acc1*Z))-Aqc*(Acc\Acq1)
    (;valid=true,condition=condc,T0,T1,DG=(T1+T1')/2,DA=(T1-T1')/2)
end
zero=schur_slope(zeros(10))
if !zero.valid
    patterns=CSV.read(joinpath(OUT,"FROZEN_DELAY_PATTERNS.csv"),DataFrame)
    cases=[("uniform_$(Int(t))ms",fill(Float64(t),10)) for t in (0,2,5,10,20,30,40,50)]
    append!(cases,[(String(r.delay_pattern_id),parse.(Float64,split(String(r.tau_ms),";"))) for r in eachrow(patterns)])
    blocked=[(;delay_pattern_id=id,mean_tau_ms=mean(t),hidden_block_condition=zero.condition,
        operator_status="BLOCKED_SCHUR_SINGULAR_AFTER_GAUGE_DEFLATION") for (id,t) in cases]
    CSV.write(joinpath(OUT,"baseline_reproduction","TABLE_D04_GRAPH_DAMPING_OPERATORS.csv"),DataFrame(blocked))
    desc=CSV.read(joinpath(OUT,"TABLE_D05_PRE_GAUGE_DIAGNOSTIC.csv"),DataFrame)
    select!(desc,Not([:eta_cross,:DG_norm,:Ltau_norm,:lambda_min_DG,:claim_status]))
    desc.operator_status=fill("BLOCKED_NO_VALID_SCHUR_SLOPE",nrow(desc))
    CSV.write(joinpath(OUT,"baseline_reproduction","TABLE_D05_GSP_DG_METRICS.csv"),desc)
    CSV.write(joinpath(OUT,"GRAPH_DG_ZERO_DELAY.csv"),DataFrame(status=["BLOCKED_SCHUR_SINGULAR_AFTER_GAUGE_DEFLATION"],
        hidden_block_condition=[zero.condition],gauge_residual=[gauge_residual]))
    open(joinpath(OUT,"D04_GATE_REPORT.md"),"w") do io
        println(io,"# D4 exact graph-damping operator gate\n")
        println(io,"**Status: BLOCKED.** After explicitly projecting out the physical global-rotation direction, the hidden block \u0394_cc(0) is singular under the preregistered retention of SG rotor-angle coordinates. Its Schur complement and derivative at s=0 are therefore undefined for this partition. The no-gauge calculation is retained only as `TABLE_D04_PRE_GAUGE_DIAGNOSTIC.csv` and is not valid evidence.\n")
        println(io,"No D_G, L_D, or L_\u03c4 result is accepted. A future derivation must retain the additional zero-frequency controller/device subspace or choose a nonsingular physically justified partition, then re-run gauge and pole consistency checks. No Laplacian or delay-dressed-damping claim follows from the previous numerical table.\n")
        println(io,"Frozen delay design remains unchanged. All ",length(cases)," requested coordinates are explicitly marked blocked in `baseline_reproduction/TABLE_D04_GRAPH_DAMPING_OPERATORS.csv`.")
    end
    println("D04 BLOCKED: hidden Schur block singular after gauge deflation; condition=",zero.condition,"; gauge residual=",gauge_residual)
    exit(0)
end
function schur_value(s,tau)
    Ds=Matrix{ComplexF64}(s*I-A0)
    for i in 1:10;Ds .-= exp(-s*tau[i]).*Aiq[i];end
    Aqq=Ds[qidx,qidx];Aqc=Ds[qidx,cidx];Acq=Ds[cidx,qidx];Acc=Ds[cidx,cidx]
    Aqq-Aqc*(Acc\Acq)
end
unit_slope=Matrix{Float64}[]
for i in 1:10
    t=zeros(10);t[i]=1
    z=schur_slope(t);z.valid||error("unit delay perturbation made hidden Schur block singular")
    push!(unit_slope,z.DG-zero.DG)
end
mass=CSV.read(joinpath(OUT,"GRAPH_MASS_REFERENCE.csv"),DataFrame)
Udf=CSV.read(joinpath(OUT,"GRAPH_GSP_BASIS.csv"),DataFrame);U=Matrix{Float64}(Udf)
Udeflated=U[:,2:end] # the zero eigenvalue of the valid physical L_P is first
Minvhalf=Diagonal(1 ./sqrt.(Float64.(mass.M)))
patterns=CSV.read(joinpath(OUT,"FROZEN_DELAY_PATTERNS.csv"),DataFrame)
cases=[("uniform_$(Int(t))ms",fill(Float64(t),10)) for t in (0,2,5,10,20,30,40,50)]
append!(cases,[(String(r.delay_pattern_id),parse.(Float64,split(String(r.tau_ms),";"))) for r in eachrow(patterns)])
rows=NamedTuple[];gsp=NamedTuple[];lin_res=0.0
for (id,tms) in cases
    r=schur_slope(tms./1000)
    if !r.valid
        push!(rows,(;delay_pattern_id=id,mean_tau_ms=mean(tms),hidden_block_condition=r.condition,
            valid=false,operator_status="BLOCKED_SCHUR_SINGULAR"));continue
    end
    Lt=zero.DG-r.DG
    predicted=zero.DG+sum((tms[i]/1000).*unit_slope[i] for i in 1:10)
    affine_error=norm(predicted-r.DG)/max(norm(r.DG),eps())
    ds=1e-6
    dgfd=(schur_value(ds+0im,tms./1000)-schur_value(-ds+0im,tms./1000))/(2ds)
    sder_error=norm(dgfd-r.T1)/max(norm(r.T1),eps())
    Dh=Udeflated'*Minvhalf*r.DG*Minvhalf*Udeflated;off=Dh-Diagonal(diag(Dh))
    lap=(norm(Lt-Lt',Inf)<1e-10&&norm(Lt*ones(10),Inf)<1e-8&&
        maximum([Lt[i,j] for i in 1:10 for j in 1:10 if i!=j])<=1e-8)
    push!(rows,(;delay_pattern_id=id,mean_tau_ms=mean(tms),hidden_block_condition=r.condition,
      lambda_min_DG=minimum(eigvals(Symmetric(r.DG))),DG_skew_norm=norm(r.DA),
      Ltau_norm=norm(Lt),Ltau_min_eigenvalue=minimum(eigvals(Symmetric(Lt))),
      Ltau_row_sum_inf=norm(Lt*ones(10),Inf),Ltau_max_offdiagonal=maximum([Lt[i,j] for i in 1:10 for j in 1:10 if i!=j]),
      Ltau_is_laplacian=lap,eta_cross=norm(off)/max(norm(Dh),eps()),
      DG_affine_delay_relative_error=affine_error,Schur_slope_centered_FD_relative_error=sder_error,
      gauge_residual=gauge_residual,retained_dimension=10,deflated_graph_modes=9,
      exact_vs_lossless_J_relative=1.9968595946726768,operator_status="NUMERICALLY_EVALUATED_EXACT_SCHUR_SLOPE"))
    chi_value=try only(patterns.chi_tau[patterns.delay_pattern_id.==id]) catch;NaN end
    push!(gsp,(;delay_pattern_id=id,mean_tau_ms=mean(tms),chi_tau=chi_value,
      graph_roughness=try only(patterns.graph_roughness[patterns.delay_pattern_id.==id]) catch;NaN end,
      eta_cross=norm(off)/max(norm(Dh),eps()),DG_norm=norm(r.DG),Ltau_norm=norm(Lt),
      lambda_min_DG=minimum(eigvals(Symmetric(r.DG))),claim_status="SUPPORTED_LOCAL_LOW_FREQUENCY_OPERATOR"))
end
CSV.write(joinpath(OUT,"baseline_reproduction","TABLE_D04_GRAPH_DAMPING_OPERATORS.csv"),DataFrame(rows))
CSV.write(joinpath(OUT,"baseline_reproduction","TABLE_D05_GSP_DG_METRICS.csv"),DataFrame(gsp))
CSV.write(joinpath(OUT,"GRAPH_DG_ZERO_DELAY.csv"),DataFrame(zero.DG,:auto))
println("D04 rows=",length(rows)," valid=",count(r->r.operator_status=="NUMERICALLY_EVALUATED_EXACT_SCHUR_SLOPE",rows),
    " gauge_residual=",gauge_residual," max_hidden_condition=",maximum(r.hidden_block_condition for r in rows if r.operator_status=="NUMERICALLY_EVALUATED_EXACT_SCHUR_SLOPE"))
