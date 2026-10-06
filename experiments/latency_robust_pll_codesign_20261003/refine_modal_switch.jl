using CSV, DataFrames, LinearAlgebra, TOML

module Scan
include(joinpath(@__DIR__,"run_l2_family_scan.jl"))
end

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
BLAS.set_num_threads(1)

function write_design(path,rho,kp,ki,eta)
    open(path,"w") do io
        println(io,"rho = [",join(rho,", "),"]")
        println(io,"Kp = [",join(kp,", "),"]")
        println(io,"Ki = [",join(ki,", "),"]")
        println(io,"eta = ",eta)
        println(io,"label = \"modal_switch_on_Z_to_N_gain_interpolation\"")
    end
end

function main()
    ctx=Scan.Roots.MegaOracle.DC.R.N.design_context(ROOT)
    z=TOML.parsefile(joinpath(@__DIR__,"designs","Z_zero_delay_tuned.toml"))
    n=TOML.parsefile(joinpath(@__DIR__,"designs","N_nominal.toml"))
    rho=Float64.(z["rho"])
    scan=CSV.read(joinpath(@__DIR__,"L2_FAMILY_SCAN.csv"),DataFrame)
    base=scan[(abs.(scan.eta.-0.8).<1e-8).&(scan.rank_by_local_crossing.<=2),:]
    rz=only(eachrow(base[base.family_label.=="Z_template_like",:]))
    rn=only(eachrow(base[base.family_label.=="N_template_like",:]))
    sz=ComplexF64(rz.critical_real+2pi*im*rz.frequency_hz)
    sn=ComplexF64(rn.critical_real+2pi*im*rn.frequency_hz)
    function crossings(eta)
        kp=exp.((1-eta).*log.(Float64.(z["Kp"])).+eta.*log.(Float64.(n["Kp"])))
        ki=exp.((1-eta).*log.(Float64.(z["Ki"])).+eta.*log.(Float64.(n["Ki"])))
        L=Scan.Roots.MegaOracle.DC.linearization(ctx,rho,kp,ki)
        m=Scan.Roots.reduced_model(L)
        cz=Scan.root_crossing(m,sz;initial=39.15)
        cn=Scan.root_crossing(m,sn;initial=39.15)
        (cz===nothing || cn===nothing) && error("root continuation failed eta=$eta")
        (;cz,cn,L,kp,ki)
    end
    lo,hi=0.8,0.9
    f_lo=crossings(lo);f_hi=crossings(hi)
    f_lo.cz.tau_ms<f_lo.cn.tau_ms && f_hi.cz.tau_ms>f_hi.cn.tau_ms ||
        error("mode-family ordering not bracketed")
    rows=NamedTuple[]
    while hi-lo>1e-5
        eta=(lo+hi)/2;r=crossings(eta)
        push!(rows,(;eta,Z_crossing_ms=r.cz.tau_ms,N_crossing_ms=r.cn.tau_ms,
            signed_gap_ms=r.cz.tau_ms-r.cn.tau_ms,Z_frequency_hz=imag(r.cz.s)/(2pi),
            N_frequency_hz=imag(r.cn.s)/(2pi)))
        if r.cz.tau_ms<r.cn.tau_ms;lo=eta;else;hi=eta;end
    end
    CSV.write(joinpath(@__DIR__,"L4_MODE_SWITCH_BISECTION.csv"),DataFrame(rows))
    eta=(lo+hi)/2;r=crossings(eta)
    crossing=(r.cz.tau_ms+r.cn.tau_ms)/2
    write_design(joinpath(@__DIR__,"designs","mode_switch_eta.toml"),rho,r.kp,r.ki,eta)
    vZ=Scan.modevector(r.L,r.cz.s,r.cz.tau_ms)
    vN=Scan.modevector(r.L,r.cn.s,r.cn.tau_ms)
    mac=abs(dot(vZ,vN))^2
    crows=NamedTuple[]
    for (side,offset) in (("below",-0.03),("above",0.03))
        tau_ms=crossing+offset
        println("L4_SWITCH_COUNT_START ",side," ",tau_ms);flush(stdout)
        try
            q=Scan.Roots.MegaOracle.trace_count(r.L,fill(tau_ms/1000,10);gamma=-0.05,step=10.0)
            good=abs(real(q.estimate)-q.nearest)<0.01 && abs(imag(q.estimate))<0.01 &&
                abs(q.phase_estimate-q.nearest)<0.01 && q.largest_phase_step<pi/2 &&
                q.quadrature_error_estimate<0.01
            push!(crows,(;side,tau_ms,status=good ? (q.nearest==0 ? "SAFE" : "UNSAFE") : "INDETERMINATE",
                root_count=good ? q.nearest : missing,trace_real=real(q.estimate),
                phase_count=q.phase_estimate,quadrature_error=q.quadrature_error_estimate,error=""))
        catch err
            push!(crows,(;side,tau_ms,status="INDETERMINATE",root_count=missing,
                trace_real=NaN,phase_count=NaN,quadrature_error=NaN,error=sprint(showerror,err)))
        end
        CSV.write(joinpath(@__DIR__,"L4_MODE_SWITCH_COUNTS.csv"),DataFrame(crows))
    end
    result=(;eta_lower=lo,eta_upper=hi,eta_estimate=eta,
        Z_family_crossing_ms=r.cz.tau_ms,N_family_crossing_ms=r.cn.tau_ms,
        crossing_gap_ms=abs(r.cz.tau_ms-r.cn.tau_ms),
        Z_frequency_hz=imag(r.cz.s)/(2pi),N_frequency_hz=imag(r.cn.s)/(2pi),
        physical_right_vector_MAC=mac,
        status=all(x.status in ("SAFE","UNSAFE") for x in crows) &&
            crows[1].root_count==0 && crows[2].root_count>=4 ?
            "NUMERICAL_TWO_FAMILY_SWITCH_WITH_FULL_CONTOUR_BRACKET" :
            "INDETERMINATE_CONTOUR")
    CSV.write(joinpath(@__DIR__,"L4_MODE_SWITCH_POINT.csv"),DataFrame([result]))
    println("L4_SWITCH_RESULT ",result)
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
