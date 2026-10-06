using CSV, DataFrames, LinearAlgebra, TOML

module DC
include(joinpath(@__DIR__, "..", "delay_dressed_replacement_frontier_20261002", "DelayCharacteristic.jl"))
end

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
BLAS.set_num_threads(1)

function main()
    ctx=DC.DelayCharacteristic.R.N.design_context(ROOT)
    baseline=CSV.read(joinpath(@__DIR__,"L0_REPRODUCTION.csv"),DataFrame)
    rows=NamedTuple[]
    for row in eachrow(baseline)
        id=String(row.design_id)
        d=TOML.parsefile(joinpath(@__DIR__,"designs",id=="Z" ?
            "Z_zero_delay_tuned.toml" : "N_nominal.toml"))
        rho,kp,ki=Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"])
        L=DC.DelayCharacteristic.linearization(ctx,rho,kp,ki)
        tau=fill(row.tau_crit_local_ms/1000,10)
        s=ComplexF64(row.critical_root_real+im*row.critical_root_imag)
        D=DC.DelayCharacteristic.delta_matrix(L,s,tau)
        decomp=svd(D);u=decomp.U[:,end];v=decomp.V[:,end]
        denom=dot(u,DC.DelayCharacteristic.delta_s(L,s,tau)*v)
        dtau=sum(DC.DelayCharacteristic.delta_tau(L,s,tau,i) for i in 1:10)
        root_tau=-dot(u,dtau*v)/denom
        real(root_tau)>0 || error("nonpositive delay slope for $id")
        for i in 1:10
            dkp=exp(-s*tau[i])*dot(u,L.Bp[:,i])*dot(L.C[:,i],v)/denom
            dki=exp(-s*tau[i])*dot(u,L.Bi[:,i])*dot(L.C[:,i],v)/denom
            push!(rows,(;design_id=id,bus=29+i,critical_delay_ms=row.tau_crit_local_ms,
                critical_frequency_hz=row.critical_delay_frequency_hz,
                alpha_tau_s_inv_per_s=real(root_tau),
                analytic_d_real_lambda_d_Kp=real(dkp),
                analytic_d_real_lambda_d_Ki=real(dki),
                analytic_d_tau_margin_ms_d_Kp=-1000*real(dkp)/real(root_tau),
                analytic_d_tau_margin_ms_d_Ki=-1000*real(dki)/real(root_tau),
                analytic_d_tau_margin_ms_d_logKp=-1000*kp[i]*real(dkp)/real(root_tau),
                analytic_d_tau_margin_ms_d_logKi=-1000*ki[i]*real(dki)/real(root_tau),
                denominator_abs=abs(denom),normalized_root_residual=norm(D*v)/max(1,norm(D)*norm(v)),
                method="EXACT_NEP_DERIVATIVE_FOR_FIXED_RHO_AND_SIMPLE_ROOT"))
        end
        CSV.write(joinpath(@__DIR__,"L2_BASELINE_ANALYTIC_GRADIENTS.csv"),DataFrame(rows))
        println("L2_GRADIENT_DONE ",id," alpha_tau=",real(root_tau));flush(stdout)
    end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
