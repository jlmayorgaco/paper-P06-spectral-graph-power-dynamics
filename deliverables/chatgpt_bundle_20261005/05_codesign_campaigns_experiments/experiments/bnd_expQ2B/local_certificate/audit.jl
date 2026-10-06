include("LocalOracle.jl")
using .LocalOracle, LinearAlgebra, TOML, CSV, DataFrames, SHA
const L=LocalOracle
seedpath=length(ARGS)>0 ? abspath(ARGS[1]) : joinpath(L.ROOT,"reports","codesign_validation_20261001","candidate_repaired.toml")
d=TOML.parsefile(seedpath)
a=L.architecture(d["support"],d);x=L.encode(a.support,d)
println("EXACT_GRADIENT_AUDIT support=",a.support," dimensions=",length(x));flush(stdout)
t=@elapsed v=L.evaluate(a,x;derivatives=true)
println("SEED J=",v.J," alpha=",v.alpha," beta=",v.beta," F=",maximum(p.value for p in v.fp),
    " Finf=",maximum(abs.(v.dc))," R=",maximum(p.value for p in v.rp)," seconds=",t);flush(stdout)
rows=NamedTuple[];matrows=NamedTuple[]
# Condition-aware centered differences on an interior nearby point.
z=clamp.(x,.025,.975);z[1:length(a.support)].*=1.08
v=L.evaluate(a,z;derivatives=true)
for j in eachindex(z),h in (1e-4,2e-5,2e-4,j<=length(a.support) ? 5e-6 : .002)
    zp=copy(z);zm=copy(z);zp[j]+=h;zm[j]-=h
    vp=L.evaluate(a,zp);vm=L.evaluate(a,zm)
    fd=(vp.g-vm.g)/(2h)
    for (label,indices) in (("modal",1:1),("beta",2:2),("steady",3:12),("peak",13:22),("rocof",23:32))
        abserr=maximum(abs.(fd[indices]-v.Jac[indices,j]))
        scale=max(maximum(abs.(fd[indices])),maximum(abs.(v.Jac[indices,j])),1e-7)
        push!(rows,(;coordinate=j,step=h,family=label,absolute_error=abserr,relative_error=abserr/scale,
            derivative_magnitude=scale))
    end
    bp=L.matrices(a,zp);bm=L.matrices(a,zm)
    for name in (:A,:B,:C,:D,:jump)
        f=(getproperty(bp,name)-getproperty(bm,name))/(2h);g=getproperty(v.ds[j],name)
        push!(matrows,(;coordinate=j,step=h,matrix=String(name),relative_error=norm(f-g)/max(norm(f),norm(g),1e-10)))
    end
    println("GRAD coordinate=",j," h=",h," maxabs=",maximum(abs.(fd-v.Jac[:,j])));flush(stdout)
end
CSV.write(joinpath(L.OUT,"TABLE_exact_gradient_audit.csv"),DataFrame(rows))
CSV.write(joinpath(L.OUT,"TABLE_exact_matrix_derivative_audit.csv"),DataFrame(matrows))
open(joinpath(L.OUT,"DERIVATIVE_AUDIT.toml"),"w") do io
    TOML.print(io,Dict("support"=>a.support,"analytic_derivatives"=>true,"phase_units"=>"radians",
        "no_PowerDynamics_in_evaluator"=>true,"source_candidate_sha"=>bytes2hex(sha256(read(seedpath))),
        "model_sha"=>"e2f104608f1eb1705f0ad2c07764e7beab0a5df4deb4d4ce4f359b1e79947f0a"))
end
println(combine(groupby(DataFrame(rows),[:family,:step]),:absolute_error=>maximum,:relative_error=>maximum))
scales=Dict("modal"=>.05,"beta"=>L.BETA,"steady"=>.5,"peak"=>.5,"rocof"=>.5)
gates=NamedTuple[]
for group in groupby(DataFrame(rows),[:family,:coordinate])
    ok=any((group.relative_error .<1e-5) .| (group.absolute_error.*scales[group.family[1]] .<1e-7))
    push!(gates,(;family=group.family[1],coordinate=group.coordinate[1],pass=ok,
        best_relative=minimum(group.relative_error),best_physical_absolute=minimum(group.absolute_error)*scales[group.family[1]]))
end
CSV.write(joinpath(L.OUT,"TABLE_gradient_gate.csv"),DataFrame(gates))
println("GRADIENT_GATE ",all(r.pass for r in gates)," passed=",count(r->r.pass,gates),"/",length(gates));flush(stdout)
