include(joinpath(@__DIR__,"..","nonlinear_codesign_20261001","ReducedDAE.jl"))
using .ReducedDAE, LinearAlgebra, CSV, DataFrames, TOML
const R=ReducedDAE
const OUT=joinpath(R.ROOT,"reports","retuning_boundary_20261002")
mkpath(OUT);BLAS.set_num_threads(1)
println("CONTEXT_START");flush(stdout)
const CTX=R.N.design_context(R.ROOT)
println("CONTEXT_READY");flush(stdout)
rows=NamedTuple[]
for rho in (.875,.90,.925,.95,.975,.99,.995,.999)
    for bus in (8,16,29),delta in (-100.,100.)
        m=R.model(CTX,fill(rho,10),fill(R.N.K0P,10),fill(R.N.K0I,10);
            bus,delta,dc_convention=:physical_supply)
        vv=R.voltage(m.x0,m;allbus=true)
        phase=atan.(vv[2:2:end],vv[1:2:end])
        jump=mod.(phase-angle.(CTX.net.voltage).+pi,2pi).-pi
        voltage=hypot.(vv[1:2:end],vv[2:2:end])
        row=(;rho,bus,delta,F0=maximum(abs,jump)/pi,R0=maximum(abs,jump)/(pi/2),
            jump_bus=argmax(abs.(jump)),Vmin=minimum(voltage),Vmax=maximum(voltage))
        push!(rows,row);println(row);flush(stdout)
    end
end
CSV.write(joinpath(OUT,"preflight_jump.csv"),DataFrame(rows))
for filename in ("reports/analytic_iteration_20261001/refined/joint_final.toml",
                 "reports/nonlinear_codesign_20261001/candidate_final_physical.toml")
    d=TOML.parsefile(joinpath(R.ROOT,filename))
    for factor in (1.,.75,.5,1.5,2.)
        kp=clamp.(factor*d["Kp"],CTX.kpmin,CTX.kpmax)
        ki=clamp.(factor*d["Ki"],CTX.kimin,CTX.kimax)
        ss=R.N.spectrum(CTX,d["rho"],kp,ki)
        println("MODAL ",filename," factor=",factor," alpha=",ss.alpha);flush(stdout)
    end
end
