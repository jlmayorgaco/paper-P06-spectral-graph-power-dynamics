using CSV, DataFrames, LinearAlgebra, NetworkDynamics, PowerDynamics
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_model_expN","PDReferenceN.jl"))
using .PDReferenceN

function main()
    println("N_TRIM_BASELINE");flush(stdout)
    base=frozen_baseline()
    k0=PDReferenceN.NOMINAL_KP; i0=PDReferenceN.NOMINAL_KI
    cases=[(0.10,k0,i0,"rho_010_nominal"),
           (0.50,k0,i0,"rho_050_nominal"),
           (0.90,k0,i0,"rho_090_nominal"),
           (0.99,k0,i0,"rho_099_nominal"),
           (0.10,0.25k0,4i0,"rho_010_Kpmin_Kimax"),
           (0.90,4k0,0.25i0,"rho_090_Kpmax_Kimin")]
    rows=NamedTuple[]
    for bus in 30:39
        println("N_TRIM_BUILD bus=",bus);flush(stdout)
        seed=zeros(10);seed[bus-29]=0.5
        nw=build_architecture(base,seed)
        for (rho,kp,ki,label) in cases
            r=zeros(10);r[bus-29]=rho
            kv=fill(k0,10);iv=fill(i0,10)
            kv[bus-29]=kp;iv[bus-29]=ki
            try
                s=trim_state(nw,base,r,kv,iv)
                residue=residual_audit(nw,s)
                powers=direct_power_audit(s,base,r)
                loc=only(eachrow(powers[powers.bus .== bus,:]))
                tol=maximum((residue.maximum,loc.max_P_error_pu,
                    loc.max_Q_error_pu,abs(loc.KCL_P_error_MW)/100,
                    abs(loc.KCL_Q_error_Mvar)/100))
                status=residue.finite && residue.maximum<1e-10 &&
                    loc.max_P_error_pu<1e-9 && loc.max_Q_error_pu<1e-9 &&
                    loc.bounds_status=="PASS" ? "PASS" : "TRIM_INFEASIBLE"
                push!(rows,(;bus,rho,Kp=kp,Ki=ki,label,status,
                    residual_max=residue.maximum,residual_norm=residue.norm,
                    P_error_pu=loc.max_P_error_pu,Q_error_pu=loc.max_Q_error_pu,
                    KCL_P_error_MW=loc.KCL_P_error_MW,
                    KCL_Q_error_Mvar=loc.KCL_Q_error_Mvar,
                    SG_tau_m_pu=loc.SG_tau_m_pu,bounds_status=loc.bounds_status,
                    error=""))
                println("N_TRIM_ROW bus=",bus," rho=",rho," status=",status,
                    " max=",tol);flush(stdout)
            catch e
                msg=sprint(showerror,e)
                push!(rows,(;bus,rho,Kp=kp,Ki=ki,label,status="TRIM_INFEASIBLE",
                    residual_max=NaN,residual_norm=NaN,P_error_pu=NaN,Q_error_pu=NaN,
                    KCL_P_error_MW=NaN,KCL_Q_error_Mvar=NaN,SG_tau_m_pu=NaN,
                    bounds_status="UNKNOWN",error=msg))
                println("N_TRIM_ERROR bus=",bus," rho=",rho," ",msg);flush(stdout)
            end
        end
    end
    path=joinpath(ROOT,"reports","experiment_N","TABLE_N04_trim_validation.csv")
    CSV.write(path,DataFrame(rows))
    println("N_TRIM_DONE rows=",length(rows)," failed=",count(row->row.status!="PASS",rows))
end

main()
