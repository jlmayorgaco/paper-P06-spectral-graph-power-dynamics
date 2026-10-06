include("ReducedDAE.jl")
using .ReducedDAE, LinearAlgebra, CSV, DataFrames, TOML, SHA
const R=ReducedDAE
const OUT=joinpath(R.OUT,"full_gfl_sector_contract")
ctx=R.N.design_context(R.ROOT)
dd=TOML.parsefile(joinpath(R.OUT,"candidate_final_physical.toml"))
m=R.model(ctx,Float64.(dd["rho"]),Float64.(dd["Kp"]),Float64.(dd["Ki"]);dc_convention=:physical_supply)
vv=R.voltage(m.x0,m)
witness_path=joinpath(OUT,"closure_witnesses.csv")
rows=CSV.read(witness_path,DataFrame)
maxerr=0.;maxusage=0.;violations=0
for row in eachrow(rows)
    donor=Int(row.donor_bus)-29;receiver=Int(row.receiver_bus)-29
    x=copy(m.x0);ix=m.gfidx[donor];op=R.N.trim_gfl(ctx,Int(row.donor_bus))
    theta=op.x[3]+row.delta
    id=op.x[8]+row.i_d;iq=op.pars.iset_q+row.i_q
    x[ix[1]]+=row.gamma_q;x[ix[2]]+=row.gamma_d
    x[ix[3]]=theta;x[ix[4]]=row.omega;x[ix[5]]=row.xi
    x[ix[6]]=sin(theta)*id+cos(theta)*iq
    x[ix[7]]=cos(theta)*id-sin(theta)*iq
    x[ix[8]]+=row.v_dc_i;x[ix[9]]+=row.v_dc
    voltage=R.voltage(x,m)
    measured=voltage[2receiver-1:2receiver]-vv[2receiver-1:2receiver]
    expected=[row.predicted_delta_ur,row.predicted_delta_ui]
    global maxerr=max(maxerr,norm(measured-expected,Inf))
    usage=norm(measured)/row.receiver_input_radius_pu
    global maxusage=max(maxusage,usage)
    global violations+=usage>1
end
@assert maxerr<1e-10
@assert violations>0
joint_path=joinpath(OUT,"joint_closure_witness.csv")
joint=CSV.read(joint_path,DataFrame)
xjoint=copy(m.x0)
for row in eachrow(joint)
    donor=Int(row.donor_bus)-29;ix=m.gfidx[donor];op=R.N.trim_gfl(ctx,Int(row.donor_bus))
    theta=op.x[3]+row.delta;id=op.x[8]+row.i_d;iq=op.pars.iset_q+row.i_q
    xjoint[ix[1]]+=row.gamma_q;xjoint[ix[2]]+=row.gamma_d
    xjoint[ix[3]]=theta;xjoint[ix[4]]=row.omega;xjoint[ix[5]]=row.xi
    xjoint[ix[6]]=sin(theta)*id+cos(theta)*iq;xjoint[ix[7]]=cos(theta)*id-sin(theta)*iq
    xjoint[ix[8]]+=row.v_dc_i;xjoint[ix[9]]+=row.v_dc
end
target=Int(joint.receiver_bus[1])-29
joint_delta=(R.voltage(xjoint,m)-vv)[2target-1:2target]
eps_target=first(filter(r->r.receiver_bus==target+29,eachrow(rows))).receiver_input_radius_pu
joint_usage=norm(joint_delta)/eps_target
@assert joint_usage>1
out=Dict("status"=>"INDEPENDENT_NONLINEAR_DAE_CONFIRMS_LOCAL_CONTRACT_CLOSURE_FAILURE",
    "cases"=>nrow(rows),"violating_pairs"=>violations,"max_input_usage"=>maxusage,
    "max_absolute_voltage_prediction_error"=>maxerr,
    "witness_csv_sha256"=>bytes2hex(sha256(read(witness_path))),
    "joint_witness_csv_sha256"=>bytes2hex(sha256(read(joint_path))),
    "joint_input_usage"=>joint_usage,"joint_delta_voltage"=>joint_delta,
    "ReducedDAE_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"ReducedDAE.jl")))),
    "script_sha256"=>bytes2hex(sha256(read(@__FILE__))),
    "scope"=>"Algebraic input-assumption failure inside product of local ellipsoids, not instability")
open(joinpath(OUT,"julia_network_witness_audit.toml"),"w") do io;TOML.print(io,out);end
TOML.print(stdout,out)
