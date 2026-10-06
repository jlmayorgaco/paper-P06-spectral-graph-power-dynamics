using NetworkDynamics
using PowerDynamics

export EquilibriumAudit, initialize_equilibrium

struct EquilibriumAudit
    pfs
    state
    powerflow_finite::Bool
    state_finite::Bool
    fixed_point::Bool
end

"Run the documented PowerDynamics two-step PF + dynamic initialization."
function initialize_equilibrium(nw; pfs = nothing, sparse = false, check = :error, tol = 1e-8)
    local pf_state
    if isnothing(pfs)
        pf_state = solve_powerflow(nw; verbose = false, sparse = sparse)
    else
        pf_state = pfs
    end
    state = initialize_from_pf!(
        nw;
        pfs = pf_state,
        verbose = false,
        sparsepf = sparse,
        check = check,
    )
    return EquilibriumAudit(
        pf_state,
        state,
        all(isfinite, uflat(pf_state)),
        all(isfinite, uflat(state)),
        isfixpoint(state; tol = tol),
    )
end
