using CSV, DataFrames, LinearAlgebra, TOML

module DC
include(joinpath(@__DIR__, "..", "delay_dressed_replacement_frontier_20261002", "DelayCharacteristic.jl"))
end

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
BLAS.set_num_threads(1)

function main()
    ctx=DC.DelayCharacteristic.R.N.design_context(ROOT)
    z=TOML.parsefile(joinpath(@__DIR__,"designs","Z_zero_delay_tuned.toml"))
    n=TOML.parsefile(joinpath(@__DIR__,"designs","N_nominal.toml"))
    rho=Float64.(z["rho"])
    scan=CSV.read(joinpath(@__DIR__,"L2_FAMILY_SCAN.csv"),DataFrame)
    rows=NamedTuple[]
    for eta in (0.0,0.5,0.9,1.0)
        kp=exp.((1-eta).*log.(Float64.(z["Kp"])).+eta.*log.(Float64.(n["Kp"])))
        ki=exp.((1-eta).*log.(Float64.(z["Ki"])).+eta.*log.(Float64.(n["Ki"])))
        L=DC.DelayCharacteristic.linearization(ctx,rho,kp,ki)
        selected=scan[(abs.(scan.eta.-eta).<1e-8).&(scan.rank_by_local_crossing.<=2),:]
        for sr in eachrow(selected)
            tau=fill(Float64(sr.local_crossing_ms)/1000,10)
            s=ComplexF64(sr.critical_real+2pi*im*sr.frequency_hz)
            D=DC.DelayCharacteristic.delta_matrix(L,s,tau)
            S=svd(D);u=S.U[:,end];v=S.V[:,end]
            denominator=dot(u,DC.DelayCharacteristic.delta_s(L,s,tau)*v)
            slope=-dot(u,sum(DC.DelayCharacteristic.delta_tau(L,s,tau,i) for i in 1:10)*v)/denominator
            real(slope)>0 || error("invalid delay slope eta=$eta rank=$(sr.rank_by_local_crossing)")
            for i in 1:10
                dkp=exp(-s*tau[i])*dot(u,L.Bp[:,i])*dot(L.C[:,i],v)/denominator
                dki=exp(-s*tau[i])*dot(u,L.Bi[:,i])*dot(L.C[:,i],v)/denominator
                push!(rows,(;eta,mode_rank=sr.rank_by_local_crossing,family_label=sr.family_label,
                    local_crossing_ms=sr.local_crossing_ms,bus=29+i,
                    d_margin_ms_d_logKp=-1000*kp[i]*real(dkp)/real(slope),
                    d_margin_ms_d_logKi=-1000*ki[i]*real(dki)/real(slope),
                    d_margin_ms_d_Kp=-1000*real(dkp)/real(slope),
                    d_margin_ms_d_Ki=-1000*real(dki)/real(slope),
                    alpha_tau_s_inv_per_s=real(slope),root_residual=norm(D*v)/max(1,norm(D)*norm(v)),
                    status="EXACT_SIMPLE_ROOT_LOCAL_DERIVATIVE;ENVELOPE_NONSMOOTH_AT_SWITCH"))
            end
        end
        CSV.write(joinpath(@__DIR__,"L2_MULTIMODE_GRADIENTS.csv"),DataFrame(rows))
        println("L2_MULTIMODE ",eta," rows=",length(rows));flush(stdout)
    end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
