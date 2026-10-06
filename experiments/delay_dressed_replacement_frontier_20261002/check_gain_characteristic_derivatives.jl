using LinearAlgebra,CSV,DataFrames
include(joinpath(@__DIR__,"DelayCharacteristic.jl"))
const D=DelayCharacteristic;const R=D.R;BLAS.set_num_threads(1)
ctx=R.N.design_context(R.ROOT);rho=fill(.875,10);kp=fill(R.N.K0P,10);ki=fill(R.N.K0I,10)
L=D.linearization(ctx,rho,kp,ki);rows=NamedTuple[]
for (parameter,vector) in (("Kp",kp),("Ki",ki)), i in (1,3,4)
    step=.01*vector[i];pp=copy(vector);pm=copy(vector);pp[i]+=step;pm[i]-=step
    Lp=parameter=="Kp" ? D.linearization(ctx,rho,pp,ki) : D.linearization(ctx,rho,kp,pp)
    Lm=parameter=="Kp" ? D.linearization(ctx,rho,pm,ki) : D.linearization(ctx,rho,kp,pm)
    Afd=(Lp.A-Lm.A)/(2step)
    Aexpect=parameter=="Kp" ? L.Bp[:,i]*L.C[:,i]' : L.Bi[:,i]*L.C[:,i]'
    A0fd=(Lp.A0-Lm.A0)/(2step)
    push!(rows,(;parameter,bus=29+i,relative_A_derivative_error=norm(Afd-Aexpect)/max(norm(Afd),eps()),
      relative_A0_change=norm(A0fd)/max(norm(Afd),eps()),state_change=norm(Lp.m.x0-Lm.m.x0)/(2step),
      quotient_change=norm(Lp.Q-Lm.Q)/(2step),Afd_norm=norm(Afd),Aexpect_norm=norm(Aexpect)))
end
CSV.write(joinpath(@__DIR__,"baseline_reproduction","TABLE_D02_GAIN_MATRIX_DIAGNOSTIC.csv"),DataFrame(rows))
println(DataFrame(rows))
