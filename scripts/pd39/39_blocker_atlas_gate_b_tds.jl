using CSV
using DataFrames
using FFTW
using NetworkDynamics
using OrdinaryDiffEqRosenbrock
using PowerDynamics
using SciMLBase
using Statistics

include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))
using .PD39

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const RESULTS = joinpath(ROOT, "results")
const CONDITIONS = CSV.read(joinpath(RESULTS, "PD39_HOLDOUT_CONDITIONS.csv"), DataFrame)
const DEFAULT_AUDIT = CSV.read(joinpath(RESULTS, "PD39_BLOCKER_DEFAULT_PATH_AUDIT.csv"), DataFrame)
const SELECTED = [
    (blocker_id = "H1", portfolio = "35;36", control = "32;36"),
    (blocker_id = "H2", portfolio = "37;38", control = "36;38"),
    (blocker_id = "H3", portfolio = "30;32;33", control = "30;33;35"),
]
const LOAD_BUS = 39
const PULSE = 0.01

parse_portfolio(p::AbstractString) = p == "none" ? Int[] : parse.(Int, split(p, ";"))
function condition_row(id)
    id == "nominal" && return (id=id, delta_pll=0.0, delta_xf=0.0, delta_cc=0.0,
        delta_load_p=0.0, delta_load_q=0.0, delta_ibr_p=0.0, branches=zeros(Float64,46))
    r=only(filter(x->String(x.condition)==id,eachrow(CONDITIONS)))
    b=[Float64(r[Symbol("delta_line_$(lpad(i,2,'0'))")]) for i in 1:46]
    (id=id,delta_pll=Float64(r.delta_pll),delta_xf=Float64(r.delta_xf),delta_cc=Float64(r.delta_cc),
        delta_load_p=Float64(r.delta_load_p),delta_load_q=Float64(r.delta_load_q),
        delta_ibr_p=Float64(r.delta_ibr_p),branches=b)
end
kwargs(c)=(controller_delta=(c.delta_pll,c.delta_xf,c.delta_cc),
    load_delta=(c.delta_load_p,c.delta_load_q),ibr_delta=c.delta_ibr_p,branch_delta=c.branches)
function proper_predecessors(p)
    bs=parse_portfolio(p)
    [join([bs[j] for j in eachindex(bs) if j != i],";") for i in eachindex(bs)]
end
function cases()
    out=NamedTuple[]
    for h in SELECTED
        push!(out,(case_id=h.blocker_id*"_blocker",role="blocker",portfolio=h.portfolio))
        for (i,p) in enumerate(proper_predecessors(h.portfolio))
            push!(out,(case_id=h.blocker_id*"_pred$(i)",role="immediate_predecessor",portfolio=p))
        end
        push!(out,(case_id=h.blocker_id*"_control",role="stable_matched_control",portfolio=h.control))
    end
    out
end

function build_tds_network_holdout(portfolio,c; pulse=PULSE)
    nw=build_confirmatory_network(parse_portfolio(portfolio);kwargs(c)...,bounds=:primary)
    vertices,edges=copy_network_components(nw)
    defaults=get_defaults_dict(vertices[LOAD_BUS])
    ps=only([s for s in keys(defaults) if occursin("Pset",string(s))])
    qs=only([s for s in keys(defaults) if occursin("Qset",string(s))])
    p0,q0=defaults[ps],defaults[qs]
    affect=(u,p,ctx)->begin
        factor=ctx.t < 1.1 ? 1+pulse : 1.0
        p[ps]=p0*factor
        p[qs]=q0*factor
    end
    cb=PresetTimeComponentCallback([1.0,1.1],ComponentAffect(affect,(),(ps,qs)))
    set_callback!(vertices[LOAD_BUS],cb)
    out=Network(vertices,edges); set_jac_prototype!(out); out
end

function voltage(s)
    try
        ur=Float64(s[VIndex(LOAD_BUS,:busbar₊u_r)])
        ui=Float64(s[VIndex(LOAD_BUS,:busbar₊u_i)])
        return hypot(ur,ui)
    catch
        return NaN
    end
end
function dominant_frequency(t,y)
    idx=findall(i->t[i]>=1.1 && isfinite(y[i]),eachindex(t))
    length(idx)<8 && return NaN
    z=Float64.(y[idx]); z .-= mean(z)
    n=length(z); dt=median(diff(Float64.(t[idx])))
    a=abs.(fft(z)); f=(0:n-1)./(n*dt)
    keep=findall(i->f[i]>=0.05 && f[i]<=5.0,1:n)
    isempty(keep) ? NaN : f[keep[argmax(a[keep])]]
end
function logfit(t,y)
    idx=findall(i->t[i]>=1.1 && isfinite(y[i]) && y[i]>1e-10,eachindex(t))
    length(idx)<5 && return (rate=NaN,r2=NaN)
    x=Float64.(t[idx]); z=log.(Float64.(y[idx]))
    xm,zm=mean(x),mean(z)
    slope=sum((x.-xm).*(z.-zm))/max(sum((x.-xm).^2),eps(Float64))
    intercept=zm-slope*xm; pred=intercept .+ slope.*x
    sse=sum((z.-pred).^2); sst=sum((z.-zm).^2)
    (rate=slope,r2=sst<=eps(Float64) ? NaN : 1-sse/sst)
end
function linear_row(p,cid)
    rs=filter(r->String(r.portfolio)==p && String(r.condition)==cid,eachrow(DEFAULT_AUDIT))
    isempty(rs) ? nothing : only(rs)
end
function one(case,cid)
    lr=linear_row(case.portfolio,cid)
    base=(case_id=case.case_id,role=case.role,portfolio=case.portfolio,condition=cid,status="failed",
        error_type="",error_message="",pulse=PULSE,linear_alpha=lr===nothing ? NaN : Float64(lr.native_alpha),
        linear_frequency_hz=lr===nothing ? NaN : Float64(lr.native_frequency_hz),
        nonlinear_rate_s_inv=NaN,nonlinear_frequency_hz=NaN,fit_r2=NaN,max_frequency_spread_hz=NaN,
        max_voltage_deviation_pu=NaN,voltage_min_pu=NaN,voltage_max_pu=NaN,
        frequency_settling_s=NaN,voltage_settling_s=NaN,sign_consistent=false,trace="")
    try
        nw=build_tds_network_holdout(case.portfolio,condition_row(cid))
        eq=initialize_equilibrium(nw;sparse=false)
        eq.powerflow_finite&&eq.state_finite&&eq.fixed_point ||
            return merge(base,(error_type="tds_equilibrium_not_qualified",
                error_message="default TDS initialization gate failed"))
        vref=voltage(eq.state)
        prob=SciMLBase.ODEProblem(nw,eq.state,(0.0,20.0))
        sol=SciMLBase.solve(prob,OrdinaryDiffEqRosenbrock.Rodas5P();
            callback=get_callbacks(nw),initializealg=SciMLBase.NoInit(),saveat=0.02,
            abstol=1e-8,reltol=1e-8)
        t=Float64.(sol.t); states=[NetworkDynamics.NWState(sol,ti) for ti in t]
        v=[voltage(s) for s in states]
        spread=Float64[]
        for s in states
            f=PD39._frequency_values(s)
            push!(spread,isempty(f) ? NaN : maximum(f)-minimum(f))
        end
        vd=abs.(v.-vref)
        combined=[max(isfinite(spread[i]) ? spread[i] : 0.0,
            isfinite(vd[i]) ? PD39.BASE_FREQ*vd[i] : 0.0) for i in eachindex(t)]
        fit=logfit(t,combined)
        trace_name="PD39_BLOCKER_TDS_TRACE_$(case.case_id)_$(cid).csv"
        CSV.write(joinpath(RESULTS,trace_name),DataFrame(time=t,frequency_spread_hz=spread,
            voltage_pu=v,voltage_deviation_pu=vd,combined_observable=combined))
        sfinite=filter(isfinite,spread); vfinite=filter(isfinite,v)
        return merge(base,(status="ok",nonlinear_rate_s_inv=fit.rate,
            nonlinear_frequency_hz=dominant_frequency(t,vd),fit_r2=fit.r2,
            max_frequency_spread_hz=isempty(sfinite) ? NaN : maximum(sfinite),
            max_voltage_deviation_pu=isempty(vfinite) ? NaN : maximum(vd[isfinite.(v)]),
            voltage_min_pu=isempty(vfinite) ? NaN : minimum(vfinite),
            voltage_max_pu=isempty(vfinite) ? NaN : maximum(vfinite),
            frequency_settling_s=PD39._settling_time(t,spread,1.1),
            voltage_settling_s=PD39._settling_time(t,vd,1.1),
            sign_consistent=isfinite(lr===nothing ? NaN : Float64(lr.native_alpha)) &&
                isfinite(fit.rate) && sign(Float64(lr.native_alpha))==sign(fit.rate),trace=trace_name))
    catch err
        merge(base,(error_type="tds_exception",error_message=sprint(showerror,err)))
    end
end

rows=NamedTuple[]
for c in cases(), cid in ("C12","nominal")
    println("B ",c.case_id," / ",cid); flush(stdout)
    push!(rows,one(c,cid))
    CSV.write(joinpath(RESULTS,"PD39_BLOCKER_TDS_VALIDATION.csv"),DataFrame(rows))
end
println("wrote TDS rows=",length(rows))
