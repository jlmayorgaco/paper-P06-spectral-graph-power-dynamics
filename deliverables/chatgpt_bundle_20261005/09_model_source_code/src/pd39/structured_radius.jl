export uncertainty_scenarios, scenario_id

"""
Fixed controller-uncertainty scenarios for the robust screen.

The factors are applied to the stock SimpleGFLDC construction before any
portfolio is initialized. The nominal case is included explicitly; the
eight corners are a preregistered Cartesian product of ±20% PLL bandwidth,
filter reactance, and current-controller bandwidth.
"""
function uncertainty_scenarios()
    scenarios = [(id = "nominal", pll_scale = 1.0, filter_scale = 1.0, current_control_scale = 1.0)]
    for pll_scale in (0.8, 1.2), filter_scale in (0.8, 1.2), current_control_scale in (0.8, 1.2)
        push!(scenarios, (
            id = scenario_id(pll_scale, filter_scale, current_control_scale),
            pll_scale = pll_scale,
            filter_scale = filter_scale,
            current_control_scale = current_control_scale,
        ))
    end
    return scenarios
end

scenario_id(pll, filter, current) = "pll$(pll)_xf$(filter)_cc$(current)"
