include("LocalOracle.jl")
using .LocalOracle, TOML, LinearAlgebra, CSV, DataFrames
const L=LocalOracle
d=TOML.parsefile(abspath(ARGS[1]));a=L.architecture(d["support"],d);x=L.encode(a.support,d)
z=clamp.(x,.025,.975);z[1:length(a.support)].*=1.08
v=L.evaluate(a,z;derivatives=true,fine=true);rows=NamedTuple[]
for h in (5e-5,1e-4,2e-4,5e-4,1e-3)
    j=17;zp=copy(z);zm=copy(z);zp[j]+=h;zm[j]-=h
    vp=L.evaluate(a,zp;fine=true);vm=L.evaluate(a,zm;fine=true)
    for k in 1:10
        fd=(vp.g[k+12]-vm.g[k+12])/(2h);an=v.Jac[k+12,j]
        push!(rows,(;coordinate=j,step=h,bus=k+29,analytic=an,finite_difference=fd,
            physical_absolute_error=.5abs(fd-an),relative_error=abs(fd-an)/max(abs(fd),abs(an),1e-7),
            time=v.fp[k].time,time_plus=vp.fp[k].time,time_minus=vm.fp[k].time))
    end
end
CSV.write(joinpath(L.OUT,"TABLE_peak_step_diagnostic.csv"),DataFrame(rows))
println(combine(groupby(DataFrame(rows),:step),:physical_absolute_error=>maximum,:relative_error=>maximum))
println("peak_times=",[p.time for p in v.fp]);println("condition_eigenvectors=",cond(v.V))
