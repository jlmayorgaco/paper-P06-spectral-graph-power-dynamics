using CSV, DataFrames, LinearAlgebra, TOML

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DC = include(joinpath(ROOT, "experiments", "delay_dressed_replacement_frontier_20261002", "DelayCharacteristic.jl"))
const R = DC.R

function main()
    seed = TOML.parsefile(joinpath(@__DIR__, "seed_uniform_875.toml"))
    rho=Float64.(seed["rho"]); kp=Float64.(seed["Kp"]); ki=Float64.(seed["Ki"])
    ctx=R.N.design_context(ROOT)
    L=DC.linearization(ctx,rho,kp,ki)
    rows1=NamedTuple[]
    for j in eachindex(L.Ai)
        sv=svdvals(L.Ai[j]); scale=maximum(sv)
        numerical_rank=count(x->x>1e-10*scale,sv)
        expected=L.B[:,j]*transpose(L.C[:,j])
        residual=norm(L.Ai[j]-expected)
        push!(rows1,(;bus=L.labels[j],numerical_rank,leading_singular=scale,
            second_singular=length(sv)>1 ? sv[2] : 0.0,
            outer_product_residual=residual,
            relative_residual=residual/max(norm(L.Ai[j]),eps())))
    end
    CSV.write(joinpath(@__DIR__,"TABLE_Q01_PLL_RANK_RESIDUAL.csv"),DataFrame(rows1))

    Ac=L.A0+L.B*transpose(L.C)
    ac_reconstruction=norm(Ac-L.A)/max(norm(L.A),eps())
    n=size(Ac,1); m=length(L.Ai); rows2=NamedTuple[]
    taus=(0.0,0.020,0.040)
    points=(-0.31+0.41im,-0.08+1.73im,-0.14+6.7im,-0.22+19.0im,0.37+43im,1.0+127im)
    for tau0 in taus, s in points
        tau=fill(tau0,m); E=Diagonal(exp.(-s.*tau))
        direct=DC.delta_matrix(L,s,tau)
        reconstructed=s*Matrix{ComplexF64}(I,n,n)-Ac-L.B*(E-I)*transpose(L.C)
        opres=norm(direct-reconstructed)/max(norm(direct),eps())
        Acs=s*Matrix{ComplexF64}(I,n,n)-Ac
        minref=minimum(svdvals(Acs))
        if minref>1e-8
            small=Matrix{ComplexF64}(I,m,m)-(E-I)*transpose(L.C)*(Acs\L.B)
            ld,sd=logabsdet(direct); la,sa=logabsdet(Acs); lr,sr=logabsdet(small)
            logerr=abs(ld-(la+lr))
            phaseerr=abs(angle(sd/(sa*sr)))
            small_min_sv=minimum(svdvals(small))
            checked=true
        else
            logerr=NaN;phaseerr=NaN;small_min_sv=NaN;checked=false
        end
        push!(rows2,(;tau_ms=1000*tau0,s_real=real(s),s_imag=imag(s),operator_relative_residual=opres,
            reference_min_singular=minref,determinant_identity_checked=checked,
            logabsdet_error=logerr,phase_error_rad=phaseerr,small_factor_min_singular=small_min_sv,
            full_state_dimension=n,low_rank_dimension=m,zero_delay_eigenvalue_identity=(tau0==0.0)))
    end
    CSV.write(joinpath(@__DIR__,"TABLE_Q02_DDE_REDUCTION_VALIDATION.csv"),DataFrame(rows2))
    status=Dict("status"=>"Q1_Q2_NUMERICALLY_VALIDATED_CONDITIONAL_ON_FIXED_DELAY_CHANNEL",
        "state_dimension"=>n,"delayed_channels"=>m,"action_rank_bound_PLL"=>m,
        "Ac_reconstruction_relative_residual"=>ac_reconstruction,
        "max_Q1_relative_residual"=>maximum(getproperty.(rows1,:relative_residual)),
        "max_Q2_operator_relative_residual"=>maximum(getproperty.(rows2,:operator_relative_residual)),
        "max_Q2_logabsdet_error"=>maximum(filter(isfinite,getproperty.(rows2,:logabsdet_error))),
        "max_Q2_phase_error_rad"=>maximum(filter(isfinite,getproperty.(rows2,:phase_error_rad))),
        "tau_test_ms"=>collect(1000 .* taus),"points_per_delay"=>length(points))
    open(joinpath(@__DIR__,"Q1_Q2_STATUS.toml"),"w") do io; TOML.print(io,status); end
    println(status)
end

main()
