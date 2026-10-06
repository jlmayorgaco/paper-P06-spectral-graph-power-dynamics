using CSV, DataFrames, LinearAlgebra, TOML

module DC
include(joinpath(@__DIR__,"..","delay_dressed_replacement_frontier_20261002","DelayCharacteristic.jl"))
end

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
BLAS.set_num_threads(1)

function main()
    length(ARGS)==1 || error("usage: evaluate_gradients.jl design.toml")
    path=abspath(ARGS[1]);id=splitext(basename(path))[1]
    folder=joinpath(@__DIR__,"evaluations",id)
    rootrows=CSV.read(joinpath(folder,"ROOTS.csv"),DataFrame)
    d=TOML.parsefile(path)
    rho,kp,ki=Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"])
    ctx=DC.DelayCharacteristic.R.N.design_context(ROOT)
    L=DC.DelayCharacteristic.linearization(ctx,rho,kp,ki)
    rows=NamedTuple[]
    for rr in eachrow(rootrows[1:min(3,nrow(rootrows)),:])
        tau=fill(Float64(rr.local_crossing_ms)/1000,10)
        s=ComplexF64(rr.critical_real+im*rr.critical_imag)
        D=DC.DelayCharacteristic.delta_matrix(L,s,tau)
        S=svd(D);u=S.U[:,end];v=S.V[:,end]
        den=dot(u,DC.DelayCharacteristic.delta_s(L,s,tau)*v)
        dsdt=-dot(u,sum(DC.DelayCharacteristic.delta_tau(L,s,tau,i) for i in 1:10)*v)/den
        real(dsdt)>0 || error("nonpositive delay derivative for $id root $(rr.rank)")
        for i in 1:10
            factor=exp(-s*tau[i])*dot(L.C[:,i],v)/den
            dlkp=factor*dot(u,L.Bp[:,i])
            dlki=factor*dot(u,L.Bi[:,i])
            push!(rows,(;candidate_id=id,mode_rank=rr.rank,
                local_crossing_ms=rr.local_crossing_ms,bus=29+i,
                d_margin_ms_d_logKp=-1000*kp[i]*real(dlkp)/real(dsdt),
                d_margin_ms_d_logKi=-1000*ki[i]*real(dlki)/real(dsdt),
                alpha_tau_s_inv_per_s=real(dsdt),
                normalized_root_residual=norm(D*v)/max(1,norm(D)*norm(v)),
                status="EXACT_SIMPLE_ROOT_LOCAL_GRADIENT"))
        end
    end
    CSV.write(joinpath(folder,"GRADIENTS.csv"),DataFrame(rows))
    println("EVALUATE_GRADIENTS_DONE ",id," modes=",maximum(getproperty.(rows,:mode_rank)));flush(stdout)
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
