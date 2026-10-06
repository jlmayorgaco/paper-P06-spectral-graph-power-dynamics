module Oracle
include(joinpath(@__DIR__,"..","analytical_delay_codesign_mega_20261002","m3_a_trace_integral.jl"))
end
using .Oracle, CSV, DataFrames, TOML, LinearAlgebra
function physical_linearization(ctx,rho,kp,ki;R=Oracle.DC.R)
    m=R.model(ctx,rho,kp,ki;dc_convention=:physical_supply)
    de=R.derivatives(m.x0,m);n=length(m.x0)
    B=zeros(n,10);C=zeros(10,n)
    for i=1:10
        ix=m.gfidx[i];s,c=sincos(m.x0[ix[3]])
        C[i,:].=-s.*de.vx[2i-1,:]+c.*de.vx[2i,:]
        C[i,ix[3]]+=-c*de.v[2i-1]-s*de.v[2i]
        B[ix[4],i]=kp[i]/m.gp[i].pll_tau;B[ix[5],i]=ki[i]
    end
    gauge,_=R.rotation_generator(m.x0,m;jacobian=false)
    gauge./=norm(gauge)
    norm(de.Fx*gauge)/norm(de.Fx)<1e-10 || error("GAUGE_PARITY_FAILED")
    norm(C*gauge)<1e-10 || error("DELAY_CHANNEL_GAUGE_FAILED")
    Q=nullspace(reshape(gauge,1,:));Bq=Q'*B;Cq=C*Q
    (;A0=Q'*(de.Fx-B*C)*Q,B=Bq,C=Matrix(Cq'))
end
function check_contour(path;ctx=Oracle.DC.R.N.design_context(Oracle.ROOT),R=Oracle.DC.R)
    started=time()
    d=TOML.parsefile(path);id=splitext(basename(path))[1]
    L=physical_linearization(ctx,Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"]);R)
    out=joinpath(@__DIR__,"spectral");mkpath(out)
    r=try
        Oracle.trace_count(L,fill(Float64(d["tau"]),10);step=10.,max_depth=22)
    catch err
        row=(;design_id=id,tau_ms=1000*d["tau"],status="INDETERMINATE",error=sprint(showerror,err),contour_seconds=time()-started)
        CSV.write(joinpath(out,id*"_contour.csv"),DataFrame([row]))
        println("CONTOUR_UNRESOLVED ",row);flush(stdout);return row
    end
    resolved=r.integer_distance<.01 && r.quadrature_error_estimate<.01 && abs(r.phase_estimate-r.nearest)<.01 && r.largest_phase_step<pi/2
    row=(;design_id=id,tau_ms=1000*d["tau"],status=resolved ? (r.nearest==0 ? "PASS_NUMERICAL" : "FAIL_ROOTS_RIGHT_OF_MARGIN") : "INDETERMINATE",
        count_real=real(r.estimate),count_imag=imag(r.estimate),phase_count=r.phase_estimate,count=r.nearest,
        integer_distance=r.integer_distance,quadrature_error=r.quadrature_error_estimate,largest_phase_step=r.largest_phase_step,
        min_sigma=r.min_sigma_small,radius=r.radius,evaluations=r.evaluations,contour_seconds=time()-started)
    CSV.write(joinpath(out,id*"_contour.csv"),DataFrame([row]))
    println("CONTOUR_RESULT ",row);flush(stdout);row
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    check_contour(abspath(ARGS[1]))
end
