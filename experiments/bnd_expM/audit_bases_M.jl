using CSV, DataFrames, NetworkDynamics, PowerDynamics
include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))
using .PD39

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const OUT = joinpath(ROOT, "reports", "experiment_M")
mkpath(joinpath(OUT, "tables"))
const model = PD39.PD39Model
const expected_ω = 2pi * 60

function defval(component, suffix)
    d = get_defaults_dict(component)
    vals = [(string(k), v) for (k, v) in d if endswith(string(k), suffix)]
    isempty(vals) && return missing
    length(vals) == 1 || error("ambiguous $suffix: $vals")
    return Float64(only(vals)[2])
end

function baserow(case, bus, name, kind, comp, expected_v)
    s = defval(comp, "systembase₊Sbase")
    ω = defval(comp, "systembase₊ωbase")
    frame = defval(comp, "systembase₊ωframe")
    v = defval(comp, "busbar₊Vbase")
    pass = !ismissing(s) && !ismissing(ω) && !ismissing(frame) &&
        isapprox(s, 100.0; atol=1e-12) && isapprox(ω, expected_ω; atol=1e-10) &&
        isapprox(frame, 1.0; atol=1e-12) &&
        (ismissing(expected_v) || (!ismissing(v) && isapprox(v, expected_v; atol=1e-12)))
    return (; case, bus, component=name, component_type=kind, Sbase=s,
        fbase_equivalent=ismissing(ω) ? missing : ω/(2pi), omega_base=ω,
        omega_frame=frame, Vbase=v, expected_Sbase=100.0,
        expected_omega_base=expected_ω, expected_Vbase=expected_v, pass)
end

function main()
    data = model.ieee39_data()
    nw = PD39.baseline_network()
    rows = NamedTuple[]
    for bus in 1:39
        row = data.bus[findfirst(==(bus), data.bus.bus), :]
        push!(rows, baserow("all_SG", bus, "bus$bus", String(row.category), nw[VIndex(bus)], Float64(row.base_kv)))
    end
    for edge in 1:length(data.branch.src_bus)
        push!(rows, baserow("all_SG", missing, "line$edge", "PiLine_fault", nw[EIndex(edge)], missing))
    end
    template = model.simple_gfldc_template()
    push!(rows, baserow("GFL_template", 38, "SimpleGFLDC", "compiled_GFL", template, missing))
    mixed = model.weighted_replacement_network(nw, 38, 0.5, 1.0)
    push!(rows, baserow("mixed_bus38", 38, "bus38", "SG_GFL", mixed[VIndex(38)], 16.5))
    CSV.write(joinpath(OUT, "tables", "TABLE_M01_component_bases.csv"), DataFrame(rows))
    for r in rows
        !r.pass && println("BASE_FAIL ", r.component, " ", r.omega_base)
    end
    println("BASE_ROWS=", length(rows), " PASS=", count(r->r.pass, rows))
    println("GFL_TEMPLATE_FBASE_HZ=", rows[end-1].fbase_equivalent)
end

main()
