include("oracle.jl")
y=parse.(Float64,split(strip(read(joinpath(OUT,"seed_faces_margin.csv"),String)),','))
J=dynamic_jacobian(y);rows=NamedTuple[]
for col in (1,4,6,8,9,10,11,15,20,21,25,29),h in (2e-5,2e-6,2e-7)
  yp=copy(y);ym=copy(y);yp[col]=min(1-1e-12,y[col]+h);ym[col]=max(0.0,y[col]-h)
  a=oracle(yp);b=oracle(ym);fd=(a.g[1:2]-b.g[1:2])/(yp[col]-ym[col])
  for row in 1:2
    push!(rows,(;coordinate=col,step=h,constraint=row==1 ? "alpha" : "beta",
      sensitivity=J[row,col],finite_difference=fd[row],
      absolute_error=abs(J[row,col]-fd[row]),
      relative_error=abs(J[row,col]-fd[row])/max(abs(J[row,col]),abs(fd[row]),1e-6)))
  end
end
CSV.write(joinpath(OUT,"dynamic_derivative_step_audit.csv"),DataFrame(rows))
println("DERIVATIVE_AUDIT ",combine(groupby(DataFrame(rows),[:constraint,:step]),:relative_error=>maximum))
