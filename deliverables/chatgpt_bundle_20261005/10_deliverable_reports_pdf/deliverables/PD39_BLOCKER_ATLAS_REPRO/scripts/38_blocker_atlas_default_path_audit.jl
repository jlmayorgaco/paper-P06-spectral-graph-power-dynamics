using CSV
using DataFrames
using LinearAlgebra
using NetworkDynamics
using PowerDynamics

include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))
using .PD39

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const RESULTS = joinpath(ROOT, "results")
const CONDITIONS = CSV.read(joinpath(RESULTS, "PD39_HOLDOUT_CONDITIONS.csv"), DataFrame)
const SELECTED = [
    (blocker_id = "H1", condition_id = "C12", portfolio = "35;36", control = "32;36"),
    (blocker_id = "H2", condition_id = "C12", portfolio = "37;38", control = "36;38"),
    (blocker_id = "H3", condition_id = "C12", portfolio = "30;32;33", control = "30;33;35"),
]

parse_portfolio(p::AbstractString) = p == "none" ? Int[] : parse.(Int, split(p, ";"))
portfolio_string(p) = isempty(p) ? "none" : join(sort(Int.(collect(p))), ";")
function condition_row(id)
    id == "nominal" && return (id=id, delta_pll=0.0, delta_xf=0.0, delta_cc=0.0,
        delta_load_p=0.0, delta_load_q=0.0, delta_ibr_p=0.0, branches=zeros(Float64,46))
    r = only(filter(x -> String(x.condition) == id, eachrow(CONDITIONS)))
    b = [Float64(r[Symbol("delta_line_$(lpad(i,2,'0'))")]) for i in 1:46]
    (id=id, delta_pll=Float64(r.delta_pll), delta_xf=Float64(r.delta_xf),
        delta_cc=Float64(r.delta_cc), delta_load_p=Float64(r.delta_load_p),
        delta_load_q=Float64(r.delta_load_q), delta_ibr_p=Float64(r.delta_ibr_p), branches=b)
end
kwargs(c) = (controller_delta=(c.delta_pll,c.delta_xf,c.delta_cc),
    load_delta=(c.delta_load_p,c.delta_load_q), ibr_delta=c.delta_ibr_p, branch_delta=c.branches)
function proper_subsets(p)
    bs=parse_portfolio(p); out=String[]
    for mask in 0:(2^length(bs)-2)
        push!(out,portfolio_string([bs[i] for i in eachindex(bs) if (mask>>(i-1))&1==1]))
    end
    unique(out)
end
function roles()
    d=Dict{String,Vector{String}}()
    function add!(p,x); d[p]=unique(vcat(get(d,p,String[]),[x])); end
    for h in SELECTED
        add!(h.portfolio,"selected_blocker_$(h.blocker_id)")
        for s in proper_subsets(h.portfolio); add!(s,"proper_subset_of_$(h.blocker_id)"); end
        add!(h.control,"stable_matched_control_for_$(h.blocker_id)")
    end
    d
end
function residual(state)
    nw=extract_nw(state); du=zeros(Float64,length(uflat(state)))
    nw(du,uflat(state),pflat(state),state.t); maximum(abs,du)
end
function vrange(state)
    v=Float64[]
    for b in 1:39
        try
            ur=Float64(state[VIndex(b,:busbar₊u_r)]); ui=Float64(state[VIndex(b,:busbar₊u_i)])
            isfinite(ur)&&isfinite(ui)&&push!(v,hypot(ur,ui))
        catch
        end
    end
    isempty(v) ? (minimum=NaN,maximum=NaN) : (minimum=minimum(v),maximum=maximum(v))
end
function critical(vals)
    z=ComplexF64.(vals); z=filter(x->isfinite(real(x))&&isfinite(imag(x))&&abs(x)>1e-8,z)
    isempty(z)&&error("no non-gauge eigenvalue"); q=z[argmax(real.(z))]
    (alpha=real(q),value=q,frequency=abs(imag(q))/(2pi))
end
function one(p,c)
    base=(portfolio=p,condition=c.id,status="failed",error_type="",error_message="",
        equilibrium_residual=NaN,voltage_min_pu=NaN,voltage_max_pu=NaN,
        fixed_point=false,native_alpha=NaN,dense_alpha=NaN,descriptor_alpha=NaN,
        native_frequency_hz=NaN,dense_frequency_hz=NaN,descriptor_frequency_hz=NaN,
        eigenvalue_native="",eigenvalue_dense="",eigenvalue_descriptor="",
        dense_minus_native=NaN,descriptor_minus_native=NaN,jacobian_condition=NaN,
        smallest_singular_value=NaN,g_z_condition=NaN,g_z_smallest_singular_value=NaN,
        eigenvector_condition=NaN,mode_family="",participation="",path="initialize_equilibrium_default")
    try
        nw=build_confirmatory_network(parse_portfolio(p);kwargs(c)...,bounds=:primary)
        eq=initialize_equilibrium(nw;sparse=false)
        v=vrange(eq.state)
        (!eq.powerflow_finite||!eq.state_finite||!eq.fixed_point) &&
            return merge(base,(error_type="equilibrium_not_qualified",error_message="default path gate failed",
                equilibrium_residual=residual(eq.state),voltage_min_pu=v.minimum,voltage_max_pu=v.maximum))
        sys=linearize_network(eq.state); red=reduce_dae(sys)
        nat=stability_audit(eq.state); den=critical(eigvals(Matrix(red.A)))
        des=if sys.M isa UniformScaling
            (critical(Matrix(sys.A))...,status="uniform_scaling_equivalent")
        else
            (critical(eigvals(Matrix(sys.A),Matrix(sys.M)))...,status="generalized_dense")
        end
        m=modal_report(eq.state); sv=svdvals(Matrix(red.A))
        return merge(base,(status="ok",equilibrium_residual=residual(eq.state),
            voltage_min_pu=v.minimum,voltage_max_pu=v.maximum,fixed_point=eq.fixed_point,
            native_alpha=nat.max_real,dense_alpha=den.alpha,descriptor_alpha=des.alpha,
            native_frequency_hz=m.critical_frequency_hz,dense_frequency_hz=den.frequency,
            descriptor_frequency_hz=des.frequency,eigenvalue_native=string(m.critical_eigenvalue),
            eigenvalue_dense=string(den.value),eigenvalue_descriptor=string(des.value),
            dense_minus_native=den.alpha-nat.max_real,descriptor_minus_native=des.alpha-nat.max_real,
            jacobian_condition=m.jacobian_condition,smallest_singular_value=m.smallest_singular_value,
            g_z_condition=m.g_z_condition,g_z_smallest_singular_value=m.g_z_smallest_singular_value,
            eigenvector_condition=m.eigenvector_condition,mode_family=m.critical_mode_family,
            participation=m.critical_participation))
    catch err
        merge(base,(error_type="default_path_exception",error_message=sprint(showerror,err)))
    end
end

rs=NamedTuple[]; rr=roles()
for p in sort(collect(keys(rr))), cid in ("C12","nominal")
    println("default ",p," / ",cid); flush(stdout)
    push!(rs,merge(one(p,condition_row(cid)),(roles=join(rr[p],"|"),)))
end
CSV.write(joinpath(RESULTS,"PD39_BLOCKER_DEFAULT_PATH_AUDIT.csv"),DataFrame(rs))
println("wrote default-path rows=",length(rs))
