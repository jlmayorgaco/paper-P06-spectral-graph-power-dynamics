using CSV, DataFrames, LinearAlgebra, Statistics
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_design_e","PhysicalData.jl"))
include(joinpath(ROOT,"src","bnd_design_e","CollectiveModel.jl"))
include(joinpath(ROOT,"src","bnd_design_e","CollectiveClosure.jl"))

domain=PhysicalData.freeze_design_domain(ROOT)
net=CollectiveModel.frozen_network(ROOT)
buses=sort(collect(keys(net.sg)));n=length(buses)
nom=PhysicalData.nominal_pll_gains()
nearest(v,w)=maximum(minimum(abs(z-y) for y in w) for z in v)
spectral_haus(v,w)=max(nearest(v,w),nearest(w,v))

single_rows=NamedTuple[]
for b in (30,33,35,37)
    rho=[x==b ? 1.0 : 0.0 for x in buses]
    model=CollectiveModel.mixed_jacobian(net,rho,fill(nom.Kp,n),fill(nom.Ki,n))
    refname=b==33 ? "C33_regenerated_Ared.csv" : "C$(b)_Ared.csv"
    ref=Matrix{Float64}(CSV.read(joinpath(ROOT,"reports","experiment_C","matrices",refname),DataFrame)[:,2:end])
    d=spectral_haus(eigvals(model.Ared),eigvals(ref))
    push!(single_rows,(bus=b,analytic_dimension=model.n_dynamic,frozen_dimension=size(ref,1),
        spectral_Hausdorff=d,pass=d<1e-7))
end

port_rows=NamedTuple[]; closure_rows=NamedTuple[]
rho=[0.15+0.7*(j-1)/(n-1) for j in 1:n]
kp=[nom.Kp*(0.4+0.23*j) for j in 1:n]
ki=[nom.Ki*(0.5+0.13*j) for j in 1:n]
for s in (0.1+0.5im,0.2+2.0im,0.5+8.0im,1.0+20.0im,2.0+50.0im)
    for (j,b) in enumerate(buses)
        op=net.gfl[b].op
        J=CollectiveModel.AnalyticGFLPLL.jacobians(op.x,op.u,op.parameters;kp=kp[j],ki=ki[j])
        direct=CollectiveModel.port_admittance(J,s)
        f=CollectiveClosure.PLLLowRank.pll_factors(op)
        rank_audit=CollectiveClosure.PLLLowRank.factor_audit(f)
        wood=CollectiveClosure.PLLLowRank.woodbury_port(f,s,kp[j],ki[j])
        err=norm(wood-direct)/max(norm(direct),eps())
        push!(port_rows,(bus=b,s_real=real(s),s_imag=imag(s),Kp=kp[j],Ki=ki[j],
            A_Kp_rank=rank_audit.rank_Akp,A_Ki_rank=rank_audit.rank_Aki,
            combined_row_rank=rank_audit.joint_rank,
            A_factor_error=rank_audit.A_factor_error,B_factor_error=rank_audit.B_factor_error,
            Woodbury_relative_error=err))
    end
    x=CollectiveClosure.closure_at(net,s,rho,kp,ki)
    a=CollectiveClosure.determinant_audit(x)
    push!(closure_rows,(s_real=real(s),s_imag=imag(s),closure_dimension=x.closure_dimension,
        network_dimension=x.network_dimension,logabsdet_error=a.logabsdet_error,
        phase_error=a.phase_error,reconstruction_error=a.reconstruction_error,
        sigma_C=a.sigma_C,sigma_T=a.sigma_T))
end

tables=joinpath(ROOT,"reports","experiment_E","tables")
CSV.write(joinpath(tables,"TABLE_E04_single_bus_spectrum_audit.csv"),DataFrame(single_rows))
CSV.write(joinpath(tables,"TABLE_E05_PLL_Woodbury_audit.csv"),DataFrame(port_rows))
CSV.write(joinpath(tables,"TABLE_E06_collective_closure_audit.csv"),DataFrame(closure_rows))
errs=DataFrame(port_rows).Woodbury_relative_error
println("SINGLE_BUS_SPECTRAL_HAUSDORFF: ",[r.spectral_Hausdorff for r in single_rows])
println("WOODBURY_MEDIAN_RELATIVE_ERROR: ",median(errs))
println("WOODBURY_P95_RELATIVE_ERROR: ",quantile(errs,0.95))
println("CLOSURE_MAX_LOGDET_ERROR: ",maximum(r.logabsdet_error for r in closure_rows))
println("CLOSURE_MAX_PHASE_ERROR: ",maximum(r.phase_error for r in closure_rows))
println("CLOSURE_DIMENSION: ",first(closure_rows).closure_dimension)

neutral_rows=NamedTuple[]
for (name,kpf,kif) in (("nominal",ones(n),ones(n)),
    ("independent_pattern_a",[0.25+0.375*(j%3) for j in 1:n],[4.0-1.5*(j%3) for j in 1:n]),
    ("independent_pattern_b",[4.0-1.5*(j%3) for j in 1:n],[0.25+0.375*(j%3) for j in 1:n]))
    model=CollectiveModel.mixed_jacobian(net,ones(n),nom.Kp.*kpf,nom.Ki.*kif)
    A=model.Ared; λ=eigvals(A)
    order=sortperm(abs.(λ)); sv=svdvals(A)
    r0=zeros(size(A,1))
    for (j,b) in enumerate(buses)
        off=9*(j-1);x=net.gfl[b].op.x
        r0[off+3]=1;r0[off+6]=x[7];r0[off+7]=-x[6]
    end
    residual=norm(A*r0)/max(norm(A)*norm(r0),eps())
    push!(neutral_rows,(gain_pattern=name,closest_real=real(λ[order[1]]),closest_imag=imag(λ[order[1]]),
        second_real=real(λ[order[2]]),second_imag=imag(λ[order[2]]),
        sigma_min=sv[end],sigma_second=sv[end-1],rotation_null_residual=residual,
        next_physical_abscissa=maximum(real.(λ[setdiff(eachindex(λ),order[1:2])]))))
end
CSV.write(joinpath(tables,"TABLE_E07_all_GFL_neutrality.csv"),DataFrame(neutral_rows))
for r in neutral_rows
    println("ALL_GFL_NEUTRALITY ",r.gain_pattern," closest=",complex(r.closest_real,r.closest_imag),
        " second=",complex(r.second_real,r.second_imag)," rotation_residual=",r.rotation_null_residual,
        " next_abscissa=",r.next_physical_abscissa)
end
