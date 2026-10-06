# Read-only checks of the current measurement and Jacobian implementations.
# No PowerDynamics or design optimization is loaded; upstream files are not edited.
using LinearAlgebra, TOML, SHA
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))

module LinearSecurity
using LinearAlgebra
function reduced_step_model(ctx,m,rho,bus;kwargs...)
    (;A=reshape([-1.0],1,1),B=[0.0],C_bus=zeros(10,1),D_bus=zeros(10))
end
function load_input_vector(ctx,m,bus;kwargs...)
    y=zeros(78); y[60]=-0.02; y
end
end
const FWPATH=joinpath(ROOT,"src","bnd_expQ2B","FiniteWindow.jl")
include(FWPATH)
ctx=(net=(voltage=ones(ComplexF64,39),),)
m=(Gy=Matrix{Float64}(I,78,78),)
t=collect(-100:400).*0.01
sig=FiniteWindow.phase_step(ctx,m,ones(10),t;load_bus=16,
    disturbance_MW=1.0,gauge_vector=x->zeros(1))
metric=only(FiniteWindow.window_metrics(sig;windows=(0.5,),dt_s=0.01,
    horizon_s=4.0,disturbance_MW=1.0))
physical_phase_jump=0.02

# Evaluate the actual constraint_jacobian function with an affine DC oracle.
# Rows 3/4 are active and the peak is exactly the steady plateau.
module JacobianHarness
const FREQUENCY_LIMIT=0.5
const ROCOF_LIMIT=0.5
decode(ctx,s,y)=(y[1:length(s)],zeros(10),zeros(10))
_bounds(ctx,s)=(zeros(length(s)+20),ones(length(s)+20))
_dc_frequency(ctx,s,y,d)=3.0*y[1]
end
const CDPATH=joinpath(ROOT,"src","bnd_expQ2B","CoreDesign.jl")
source=read(CDPATH,String)
start=first(findfirst("function constraint_jacobian(",source))
stop=first(findnext("function _solve_qp(",source,start))-1
include_string(JacobianHarness,source[start:stop],CDPATH)
e=(g=[-10.0,-10.0,0.0,0.0],rho=ones(10),Fpeak=0.6,Finf=0.6)
j=JacobianHarness.constraint_jacobian(nothing,[30],fill(0.2,21),e)

d=TOML.parsefile(joinpath(ROOT,"reports","experiment_Q2B","CORE",
    "Z_Q2B_CORE_d100_restored.toml"))
kp0=10pi; ki0=kp0^2/4
y=vcat(d["epsilon"],(d["Kp"].-0.25kp0)./(3.75kp0),
    (d["Ki"].-0.25ki0)./(3.75ki0))
lo=vcat(fill(1e-5,10),zeros(20)); hi=vcat(fill(1-1e-5,10),ones(20))

results=Dict(
    "scope"=>"Synthetic source-level reproduction; not an IEEE39 candidate rerun",
    "FiniteWindow_sha256"=>bytes2hex(sha256(read(FWPATH))),
    "CoreDesign_sha256"=>bytes2hex(sha256(read(CDPATH))),
    "phase_jump_expected_rad"=>physical_phase_jump,
    "phase_jump_returned_as_rad"=>sig.phase_jump_rad[1],
    "phase_jump_ratio_returned_to_expected"=>sig.phase_jump_rad[1]/physical_phase_jump,
    "frequency_jump_expected_Hz"=>physical_phase_jump/(2pi*0.5),
    "frequency_jump_returned_Hz"=>metric.F_peak_Hz,
    "rocof_jump_expected_Hz_s"=>physical_phase_jump/(2pi*0.5^2),
    "rocof_jump_returned_Hz_s"=>metric.R_peak_Hz_s,
    "phase_jump_units_match"=>isapprox(sig.phase_jump_rad[1],physical_phase_jump;rtol=1e-12),
    "dc_gradient"=>j.J[3,1],
    "plateau_peak_gradient"=>j.J[4,1],
    "plateau_gradients_match"=>isapprox(j.J[3,1],j.J[4,1];atol=1e-10),
    "current_candidate_minimum_normalized_box_distance"=>minimum(min.(y.-lo,hi.-y)),
    "box_omission_explains_current_candidate_KKT_residual"=>false,
    "SG_support_count"=>2^10,
    "three_state_architecture_strata_count"=>3^10,
)
open(joinpath(@__DIR__,"isolated_check_results.toml"),"w") do io
    TOML.print(io,results)
end
TOML.print(stdout,results)
