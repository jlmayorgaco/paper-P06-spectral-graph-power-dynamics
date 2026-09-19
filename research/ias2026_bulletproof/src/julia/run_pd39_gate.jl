"""PowerDynamics IEEE-39 Gate A.

This is a reproducible package-validation gate, not SG->GFL parity. It
exercises the official IEEE-39 tutorial model over several nonlinear solver
tolerances and deterministic initial guesses, then checks both dynamic
initialization paths, residuals, eigenvalue-path consistency, and Jacobian
conditioning. The POWERDYNAMICS_VALIDATED label is emitted only when all of
those checks pass.
"""

using PowerDynamics
using PowerDynamics.Library
using ModelingToolkitBase
using NetworkDynamics
using LinearAlgebra

example_dir = joinpath(pkgdir(PowerDynamics), "docs", "examples")
include(joinpath(example_dir, "ieee39_part1.jl"))

formula = @initformula :ZIPLoad₊Vset = sqrt(:busbar₊u_r^2 + :busbar₊u_i^2)
set_initformula!(nw[VIndex(31)], formula)
set_initformula!(nw[VIndex(39)], formula)

num_buses = length(nw.im.g.fadjlist)
num_branches = nw.im.g.ne
@assert num_buses == 39
@assert num_branches == 46

const TOLS = [1e-8, 1e-10, 1e-12]
const RESIDUAL_LIMIT = 1e-8
const SPECTRUM_LIMIT = 1e-6

function state_residual(net, state)
    du = zeros(Float64, length(uflat(state)))
    net(du, uflat(state), pflat(state), 0.0)
    return maximum(abs, du)
end

function make_initial_state(pfnw, kind::Symbol)
    state = kind === :default ? NWState(pfnw) : NWState(pfnw; guess=true)
    if kind === :perturbed
        for i in eachindex(uflat(state))
            if isfinite(uflat(state)[i])
                uflat(state)[i] += 0.01 * sin(i)
            end
        end
    end
    return state
end

function spectrum_key(z)
    return (real(z), imag(z))
end

function sorted_spectrum(state)
    return sort(collect(jacobian_eigenvals(state)); by=spectrum_key)
end

function spectrum_delta(a, b)
    length(a) == length(b) || return Inf
    isempty(a) && return 0.0
    return maximum(abs.(a .- b))
end

function csvquote(x)
    s = string(x)
    return occursin(',', s) || occursin('"', s) ? "\"" * replace(s, "\"" => "\"\"") * "\"" : s
end

function run_gate()
    campaign_root = normpath(joinpath(@__DIR__, "..", ".."))
    raw_dir = joinpath(campaign_root, "raw", "powerdynamics")
    mkpath(raw_dir)

pfnw = powerflow_model(nw)
rows = NamedTuple[]
reference_spectrum = nothing
reference_pfs = nothing
high_accuracy_pfs = nothing
all_pass = true
path_count = 0

for tol in TOLS
    for guess_kind in (:default, :guess, :perturbed)
        path_count += 1
        path_name = string("nonmutating/", guess_kind)
        pf_ok = false
        init_ok = false
        pf_resid = Inf
        dyn_resid = Inf
        spectral_delta = Inf
        jac_cond = Inf
        jac_sigma_min = 0.0
        eig_count = 0
        spectral_abscissa = Inf
        message = ""
        try
            pfs0 = make_initial_state(pfnw, guess_kind)
            pfs = solve_powerflow(nw; pfnw, pfs0, verbose=false,
                                  use_guesses=true, tol=tol,
                                  abstol=tol, reltol=tol, maxiters=10000)
            pf_resid = state_residual(pfnw, pfs)
            pf_ok = isfinite(pf_resid) && pf_resid <= RESIDUAL_LIMIT
            s0 = initialize_from_pf(nw; pfnw, pfs=pfs, verbose=false,
                                    subverbose=false, check=:none,
                                    tol=tol, nwtol=tol)
            dyn_resid = state_residual(nw, s0)
            sys = linearize_network(s0)
            spectrum = sorted_spectrum(s0)
            reference_spectrum = isnothing(reference_spectrum) ? spectrum : reference_spectrum
            spectral_delta = spectrum_delta(spectrum, reference_spectrum)
            jac_matrix = Matrix(sys.A)
            jac_cond = cond(jac_matrix)
            jac_sigma_min = minimum(svdvals(jac_matrix))
            eig_count = length(spectrum)
            spectral_abscissa = maximum(real, spectrum)
            init_ok = isfinite(dyn_resid) && dyn_resid <= RESIDUAL_LIMIT &&
                      isfinite(spectral_delta) && spectral_delta <= SPECTRUM_LIMIT &&
                      isfinite(jac_cond) && isfinite(jac_sigma_min)
            if isnothing(reference_pfs)
                reference_pfs = pfs
            end
            if tol == minimum(TOLS)
                high_accuracy_pfs = pfs
            end
            message = "ok"
        catch err
            message = sprint(showerror, err)
        end
        row_ok = pf_ok && init_ok
        all_pass &= row_ok
        push!(rows, (path=path_name, tolerance=tol, pf_ok=pf_ok,
                     init_ok=init_ok, row_status=row_ok ? "PASS" : "FAIL",
                     pf_residual=pf_resid, dynamic_residual=dyn_resid,
                     spectrum_delta=spectral_delta, jacobian_condition=jac_cond,
                     jacobian_sigma_min=jac_sigma_min, eigenvalue_count=eig_count,
                     spectral_abscissa=spectral_abscissa, message=message))
    end
end

# Exercise the mutating componentwise initialization path independently.
mutating_ok = false
mutating_resid = Inf
mutating_delta = Inf
mutating_message = ""
try
    @assert !isnothing(high_accuracy_pfs)
    mutating_state = initialize_from_pf!(nw; pfnw, pfs=high_accuracy_pfs,
                                         verbose=false, subverbose=false,
                                         check=:none, tol=minimum(TOLS),
                                         nwtol=minimum(TOLS))
    mutating_resid = state_residual(nw, mutating_state)
    mutating_spectrum = sorted_spectrum(mutating_state)
    mutating_delta = spectrum_delta(mutating_spectrum, reference_spectrum)
    mutating_ok = isfinite(mutating_resid) && mutating_resid <= RESIDUAL_LIMIT &&
                   isfinite(mutating_delta) && mutating_delta <= SPECTRUM_LIMIT
    mutating_message = "ok"
catch err
    mutating_message = sprint(showerror, err)
end
all_pass &= mutating_ok

csv_path = joinpath(raw_dir, "gate_a_paths.csv")
open(csv_path, "w") do io
    headers = ["path", "tolerance", "pf_ok", "init_ok", "row_status",
               "pf_residual", "dynamic_residual", "spectrum_delta",
               "jacobian_condition", "jacobian_sigma_min", "eigenvalue_count",
               "spectral_abscissa", "message"]
    println(io, join(headers, ','))
    for row in rows
        println(io, join(map(h -> csvquote(getproperty(row, Symbol(h))), headers), ','))
    end
    mutating_values = ["mutating_componentwise", missing, missing, missing,
                       mutating_ok ? "PASS" : "FAIL", missing,
                       mutating_resid, mutating_delta, missing, missing,
                       missing, missing, mutating_message]
    println(io, join(map(csvquote, mutating_values), ','))
end

report_path = joinpath(raw_dir, "pd39_equilibrium_gate.md")
open(report_path, "w") do io
    println(io, "# PowerDynamics IEEE-39 Gate A")
    println(io)
    println(io, "status: ", all_pass ? "POWERDYNAMICS_VALIDATED" : "STOPPED_BY_GATE")
    println(io, "gate: A")
    println(io, "scope: official PowerDynamics tutorial model; not SG->GFL parity")
    println(io, "package_version: 5.0.0")
    println(io, "julia_version: ", VERSION)
    println(io, "buses: ", num_buses)
    println(io, "branches: ", num_branches)
    println(io, "nonmutating_paths: ", path_count)
    println(io, "tolerances: ", join(TOLS, ", "))
    println(io, "initial_guess_paths: default, guess, perturbed")
    println(io, "residual_limit: ", RESIDUAL_LIMIT)
    println(io, "spectrum_consistency_limit: ", SPECTRUM_LIMIT)
    println(io, "mutating_componentwise_status: ", mutating_ok ? "PASS" : "FAIL")
    println(io, "mutating_componentwise_residual: ", mutating_resid)
    println(io, "mutating_componentwise_spectrum_delta: ", mutating_delta)
    println(io, "paths_csv: raw/powerdynamics/gate_a_paths.csv")
    println(io, "source: PowerDynamics official docs/examples/ieee39_part1.jl")
    if !isempty(mutating_message) && mutating_message != "ok"
        println(io, "mutating_componentwise_message: ", mutating_message)
    end
end

    println(all_pass ? "POWERDYNAMICS_GATE_A_PASS" : "POWERDYNAMICS_GATE_A_STOPPED",
            " paths=", path_count + 1, " residual_limit=", RESIDUAL_LIMIT,
            " spectrum_limit=", SPECTRUM_LIMIT)
end

run_gate()
