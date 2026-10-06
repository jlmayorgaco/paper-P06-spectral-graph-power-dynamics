rows=NamedTuple[]
for j in (1,2,3,4,5,7,17,18,24), h in (1e-5,1e-4)
    xp=copy(x);xm=copy(x);xp[j]+=h;xm[j]-=h
    vp=L.evaluate(a,xp);vm=L.evaluate(a,xm)
    fd=(vp.g-vm.g)/(2h)
    for k in (1,2,3,12,13,22)
        push!(rows,(;coordinate=j,constraint=k,h,analytic=v.Jac[k,j],fd=fd[k],
            abserr=abs(v.Jac[k,j]-fd[k]),relerr=abs(v.Jac[k,j]-fd[k])/max(abs(fd[k]),abs(v.Jac[k,j]),1e-8)))
    end
end
CSV.write(joinpath(L.OUT,"GRADIENTS_AT_STALLED_POINT.csv"),DataFrame(rows))
println(combine(groupby(DataFrame(rows),[:constraint,:h]),:abserr=>maximum,:relerr=>maximum))
