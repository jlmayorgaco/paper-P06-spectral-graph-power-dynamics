using LinearAlgebra, ForwardDiff, CSV, DataFrames, TOML, SHA, Dates
include(joinpath(@__DIR__,"..","delay_dressed_replacement_frontier_20261002","DelayCharacteristic.jl"))
const DC=DelayCharacteristic
const R=DC.R
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
BLAS.set_num_threads(1)
function write_matrix(path,a)
    CSV.write(path,DataFrame(a,:auto))
end
function main()
    ctx=R.N.design_context(ROOT)
    inputs=[("N",joinpath(ROOT,"experiments","latency_robust_pll_codesign_20261003","designs","N_nominal.toml")),
        ("Z",joinpath(ROOT,"experiments","latency_robust_pll_codesign_20261003","designs","Z_zero_delay_tuned.toml")),
        ("T",joinpath(ROOT,"experiments","latency_robust_pll_closure_20261003","designs","full20_step15_medium.toml"))]
    gates=NamedTuple[]
    for (id,path) in inputs
        d=TOML.parsefile(path);rho,kp,ki=Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"])
        m=R.model(ctx,rho,kp,ki;dc_convention=:physical_supply)
        f=zeros(length(m.x0));R.rhs!(f,m.x0,m,0.)
        A=R.derivatives(m.x0,m).Fx
        Acheck=ForwardDiff.jacobian(m.x0) do x
            dx=similar(x);R.rhs!(dx,x,m,0.);dx
        end
        C=ForwardDiff.jacobian(m.x0) do x
            v=R.voltage(x,m)
            [-sin(x[m.gfidx[i][3]])*v[2i-1]+cos(x[m.gfidx[i][3]])*v[2i] for i=1:10]
        end
        B=zeros(length(m.x0),10)
        for i=1:10
            B[m.gfidx[i][4],i]=kp[i]/m.gp[i].pll_tau
            B[m.gfidx[i][5],i]=ki[i]
        end
        old=DC.linearization(ctx,rho,kp,ki)
        ev=eigvals(A);deleteat!(ev,argmin(abs.(ev)))
        evold=eigvals(old.A)
        eigdiff=max(maximum(minimum(abs.(z.-evold)) for z in ev),maximum(minimum(abs.(z.-ev)) for z in evold))
        out=joinpath(@__DIR__,"models",id);mkpath(out)
        saved=joinpath(out,"design.toml")
        if isfile(saved)
            sha256(read(saved))==sha256(read(path)) || error("frozen design copy differs")
        else
            cp(path,saved)
        end
        for (name,mat) in (("A",A),("A0",A-B*C),("B",B),("C",C),("x0",reshape(m.x0,:,1)))
            write_matrix(joinpath(out,name*".csv"),mat)
        end
        ports=DataFrame(bus=30:39,state_index=[ix[end-1]-1 for ix in m.sgidx],
            M=[2*m.sp[i].inertia*(1-rho[i])*m.sp[i].rating_mva for i=1:10],
            rho=rho,Kp=kp,Ki=ki,P0=ctx.power,
            pll_angle_index=[ix[3]-1 for ix in m.gfidx],pll_frequency_index=[ix[4]-1 for ix in m.gfidx])
        CSV.write(joinpath(out,"ports.csv"),ports)
        CSV.write(joinpath(out,"spectrum_tau0.csv"),DataFrame(real=real.(ev),imag=imag.(ev)))
        push!(gates,(;design=id,equilibrium_inf=norm(f,Inf),jacobian_relative=norm(A-Acheck)/norm(Acheck),
            eigenvalue_max_set_difference=eigdiff,alpha=maximum(real.(ev)),
            GFL_MW=dot(ctx.power,rho),SG_MW=dot(ctx.power,1 .-rho),
            input_sha256=bytes2hex(sha256(read(path)))))
        println("EXPORTED ",last(gates));flush(stdout)
    end
    CSV.write(joinpath(@__DIR__,"TABLE_01_PARITY.csv"),DataFrame(gates))
end
main()
