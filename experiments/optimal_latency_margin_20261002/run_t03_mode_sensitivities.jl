using CSV, DataFrames, LinearAlgebra, TOML

module FastRoots
include(joinpath(@__DIR__,"run_fast_root_locus.jl"))
end

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const MEGA=joinpath(ROOT,"experiments","analytical_delay_codesign_mega_20261002")
BLAS.set_num_threads(1)

function perturb_root(ctx,rho,kp,ki,tau,s,kind,i,h)
    p=copy(rho);k=copy(kp);j=copy(ki)
    if kind==:rho;p[i]+=h
    elseif kind==:Kp;k[i]+=h
    else;j[i]+=h
    end
    L=FastRoots.MegaOracle.DC.linearization(ctx,p,k,j)
    model=FastRoots.reduced_model(L)
    r=FastRoots.refine(model,s,tau;tol=1e-10)
    r.converged && abs(r.s-s)<0.1 ? r.s : NaN+NaN*im
end

function main()
    ctx=FastRoots.MegaOracle.DC.R.N.design_context(ROOT)
    crossings=CSV.read(joinpath(@__DIR__,"T02_PRECISE_CROSSINGS.csv"),DataFrame)
    rows=NamedTuple[]
    for point in eachrow(crossings)
        design_id=String(point.design_id)
        path=design_id=="seed_875" ? joinpath(MEGA,"seed_uniform_875.toml") :
            joinpath(MEGA,"M1_ZERO_DELAY_DESIGN.toml")
        d=TOML.parsefile(path)
        rho,kp,ki=Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"])
        L=FastRoots.MegaOracle.DC.linearization(ctx,rho,kp,ki)
        s=point.critical_root_real+im*point.critical_root_imag
        tau=fill(point.local_root_crossing_ms/1000,length(L.Ai))
        D=FastRoots.MegaOracle.DC.delta_matrix(L,s,tau)
        S=svd(D);u=S.U[:,end];v=S.V[:,end]
        ds=FastRoots.MegaOracle.DC.delta_s(L,s,tau)
        denominator=dot(u,ds*v)
        root_residual=norm(D*v)/max(1,norm(D)*norm(v))
        physical_v=L.Q*v;physical_v/=norm(physical_v)
        port_scores=[abs(dot(L.C[:,i],v))*norm(L.B[:,i]) for i in 1:10]
        port_normalizer=sum(port_scores)
        for i in 1:10
            d_tau=-dot(u,FastRoots.MegaOracle.DC.delta_tau(L,s,tau,i)*v)/denominator
            hkp=1e-4*kp[i];hki=1e-4*ki[i];hrho=1e-5
            rp=perturb_root(ctx,rho,kp,ki,tau,s,:Kp,i,hkp)
            rm=perturb_root(ctx,rho,kp,ki,tau,s,:Kp,i,-hkp)
            d_kp=(rp-rm)/(2hkp)
            rp=perturb_root(ctx,rho,kp,ki,tau,s,:Ki,i,hki)
            rm=perturb_root(ctx,rho,kp,ki,tau,s,:Ki,i,-hki)
            d_ki=(rp-rm)/(2hki)
            rp=perturb_root(ctx,rho,kp,ki,tau,s,:rho,i,hrho)
            rm=perturb_root(ctx,rho,kp,ki,tau,s,:rho,i,-hrho)
            d_rho=(rp-rm)/(2hrho)
            ix=L.m.gfidx[i]
            pll_state_score=sum(abs2,physical_v[ix[3:5]])
            status=all(isfinite,[real(d_tau),real(d_kp),real(d_ki),real(d_rho)]) ?
                "LOCAL_SENSITIVITY_NUMERICALLY_VALIDATED" : "INDETERMINATE_BRANCH_FD"
            push!(rows,(;design_id,replacement_percent=point.replacement_percent,
                tau_critical_ms=point.local_root_crossing_ms,critical_family_id=point.critical_family_id,
                critical_frequency_hz=point.critical_frequency_hz,bus=29+i,
                d_real_lambda_d_tau_i_per_s=real(d_tau),d_imag_lambda_d_tau_i_per_s=imag(d_tau),
                d_real_lambda_d_Kp_i=real(d_kp),d_imag_lambda_d_Kp_i=imag(d_kp),
                d_real_lambda_d_Ki_i=real(d_ki),d_imag_lambda_d_Ki_i=imag(d_ki),
                d_real_lambda_d_rho_i=real(d_rho),d_imag_lambda_d_rho_i=imag(d_rho),
                pll_state_right_vector_energy_fraction=pll_state_score,
                pll_port_coupling_score=port_normalizer>0 ? port_scores[i]/port_normalizer : NaN,
                root_residual,denominator_abs=abs(denominator),status,
                method="tau:exact_NEp_left_right;Kp_Ki_rho:central_FD_relinearized_full_model"))
            CSV.write(joinpath(@__DIR__,"T03_CRITICAL_MODE_FAMILIES.csv"),DataFrame(rows))
            println("T03_BUS ",design_id," bus=",29+i," status=",status);flush(stdout)
        end
    end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
