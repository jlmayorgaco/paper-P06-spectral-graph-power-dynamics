using CSV, DataFrames, LinearAlgebra, TOML

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const DC=include(joinpath(ROOT,"experiments","delay_dressed_replacement_frontier_20261002","DelayCharacteristic.jl"))
const NEV=include(joinpath(ROOT,"experiments","delay_dressed_replacement_frontier_20261002","DDE_NEV.jl"))
const R=DC.R

function characteristic(L,Ai,s,tau)
    D=Matrix{ComplexF64}(s*I-L.A0)
    for j in eachindex(Ai);D.-=exp(-s*tau[j]).*Ai[j];end
    D
end

function main()
    seed=TOML.parsefile(joinpath(@__DIR__,"seed_uniform_875.toml"))
    rho=Float64.(seed["rho"]);kp0=Float64.(seed["Kp"]);ki0=Float64.(seed["Ki"])
    ctx=R.N.design_context(ROOT);L=DC.linearization(ctx,rho,kp0,ki0)
    sigma=0.05; buses=[30,32,34,36,38];freqs=[0.017,0.5,3.0];delays=[0.0,0.020,0.040]
    rows=NamedTuple[]; kpmin=ctx.kpmin[1];kpmax=ctx.kpmax[1];kimin=ctx.kimin[1];kimax=ctx.kimax[1]
    for tau0 in delays, bus in buses, f in freqs
        j=bus-29;sb=-sigma+2pi*im*f;tau=fill(tau0,length(L.Ai))
        Dminus=Matrix{ComplexF64}(sb*I-L.A0)
        for k in eachindex(L.Ai);k==j && continue;Dminus.-=exp(-sb*tau[k]).*L.Ai[k];end
        c=L.C[:,j];gp=transpose(c)*(Dminus\L.Bp[:,j]);gi=transpose(c)*(Dminus\L.Bi[:,j])
        G=[real(gp) real(gi);imag(gp) imag(gi)]
        sv=svdvals(G);authority=minimum(sv);condition=maximum(sv)/max(authority,eps())
        rhs=exp(sb*tau0);gains=try G\[real(rhs),imag(rhs)] catch; [NaN,NaN] end
        kps,kis=gains
        gain_solve_res=all(isfinite,gains) ? norm(G*gains-[real(rhs),imag(rhs)]) : Inf
        kp=copy(kp0);ki=copy(ki0);kp[j]=kps;ki[j]=kis
        Ai=[(kp[k].*L.Bp[:,k].+ki[k].*L.Bi[:,k])*transpose(L.C[:,k]) for k in eachindex(L.Ai)]
        Delta=characteristic(L,Ai,sb,tau)
        F=svd(Delta);v0=ComplexF64.(F.V[:,end]);nepres=minimum(F.S)/max(1,norm(Delta))
        Anew=L.A0+sum(Ai;init=zeros(size(L.A0)))
        Lnew=(;A=Anew,A0=L.A0,Ai)
        root=NEV.newton_root(Lnew,sb,tau,v0;tol=2e-11,maxiter=20,max_step=1.0)
        root_re_error=root.converged ? abs(real(root.s)+sigma) : NaN
        root_im_error=root.converged ? abs(imag(root.s)-2pi*f) : NaN
        bounds=all(isfinite,gains) && kpmin<=kps<=kpmax && kimin<=kis<=kimax
        push!(rows,(;bus,tau_ms=1000*tau0,target_frequency_Hz=f,target_real_s_inv=-sigma,
            g_p_real=real(gp),g_p_imag=imag(gp),g_I_real=real(gi),g_I_imag=imag(gi),
            sigma_min_G=authority,condition_G=condition,Kp_closed_form=kps,Ki_closed_form=kis,
            gains_within_frozen_bounds=bounds,gain_equation_residual=gain_solve_res,
            full_NEp_relative_residual=nepres,root_converged=root.converged,
            root_real_boundary_error=root_re_error,root_imag_boundary_error=root_im_error,
            root_NEp_residual=root.residual,claim_scope="conditional target-root branch; not a global nearest-pole certificate"))
    end
    CSV.write(joinpath(@__DIR__,"TABLE_Q04_PI_CLOSED_FORM_VALIDATION.csv"),DataFrame(rows))
    status=Dict("status"=>"CONDITIONAL_PI_LAW_VALIDATED_ON_FULL_CHARACTERISTIC",
        "tests"=>length(rows),"buses"=>buses,"tau_ms"=>1000 .*delays,"target_frequency_Hz"=>freqs,
        "max_gain_equation_residual"=>maximum(getproperty.(rows,:gain_equation_residual)),
        "max_full_NEp_residual"=>maximum(getproperty.(rows,:full_NEp_relative_residual)),
        "max_real_boundary_error"=>maximum(filter(isfinite,getproperty.(rows,:root_real_boundary_error))),
        "max_imag_boundary_error"=>maximum(filter(isfinite,getproperty.(rows,:root_imag_boundary_error))),
        "min_PI_authority"=>minimum(getproperty.(rows,:sigma_min_G)),
        "max_PI_condition"=>maximum(getproperty.(rows,:condition_G)),
        "all_target_roots_converged"=>all(getproperty.(rows,:root_converged)))
    open(joinpath(@__DIR__,"Q4_STATUS.toml"),"w") do io; TOML.print(io,status); end
    println(status)
end

main()
