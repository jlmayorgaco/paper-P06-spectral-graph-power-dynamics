using LinearAlgebra,CSV,DataFrames,TOML
include(joinpath(@__DIR__,"DelayCharacteristic.jl"))
const D=DelayCharacteristic;const R=D.R
ROOT=R.ROOT;OUT=joinpath(@__DIR__,"baseline_reproduction");BLAS.set_num_threads(1)
ctx=R.N.design_context(ROOT)
designs=[("baseline",fill(.875,10),fill(R.N.K0P,10),fill(R.N.K0I,10)),
    ("joint_final",Float64.(TOML.parsefile(joinpath(@__DIR__,"inputs","prior_joint_input.toml"))["rho"]),
     Float64.(TOML.parsefile(joinpath(@__DIR__,"inputs","prior_joint_input.toml"))["Kp"]),
     Float64.(TOML.parsefile(joinpath(@__DIR__,"inputs","prior_joint_input.toml"))["Ki"]))]
rows=NamedTuple[];poles=NamedTuple[]
for (name,rho,kp,ki) in designs
    L=D.linearization(ctx,rho,kp,ki);vals=eigen(ComplexF64.(L.A));dtau=zeros(length(L.Ai))
    errs=[D.residual(L,vals.values[j],dtau,vals.vectors[:,j]) for j in eachindex(vals.values)]
    reconstruction=norm(L.A-(L.A0+sum(L.Ai;init=zeros(size(L.A)))))/max(1,norm(L.A))
    push!(rows,(;design=name,state_dimension=size(L.A,1),delayed_PLL_channels=length(L.Ai),
      A_reconstruction_relative_error=reconstruction,max_tau0_NEV_residual=maximum(errs),
      tau0_operator="Delta(s,0)=sI-A exactly",
      spectral_abscissa=maximum(real,vals.values),status="EXACT_IDENTITY_AT_TAU_ZERO"))
    for j in eachindex(vals.values)
        z=vals.values[j]
        push!(poles,(;design=name,delay_pattern_id="uniform_tau_0",tau_vector_ms=join(fill(0.0,length(L.Ai)),";"),
          pole_id=j,real=real(z),imag=imag(z),frequency_hz=abs(imag(z))/(2pi),
          damping_ratio=-real(z)/abs(z),residual=errs[j],solver="tau_zero_matrix_identity/eig",root_tracking_id="ode_$j"))
    end
end
CSV.write(joinpath(OUT,"TABLE_D01_EXACT_DDE_POLES.csv"),DataFrame(poles))
CSV.write(joinpath(OUT,"TABLE_D01_TAU0_IDENTITY.csv"),DataFrame(rows))
println(DataFrame(rows))
