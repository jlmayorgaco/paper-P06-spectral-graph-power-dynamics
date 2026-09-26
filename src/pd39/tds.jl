using OrdinaryDiffEqRosenbrock
using SciMLBase

export simulate_equilibrium

"Run an unforced holdout trajectory from a qualified equilibrium."
function simulate_equilibrium(nw, s0; tspan = (0.0, 5.0), saveat = 0.01)
    prob = SciMLBase.ODEProblem(nw, s0, tspan)
    sol = SciMLBase.solve(prob, OrdinaryDiffEqRosenbrock.Rodas5P(); saveat = saveat)
    return sol
end
