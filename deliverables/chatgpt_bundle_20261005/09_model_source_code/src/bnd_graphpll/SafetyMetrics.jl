module SafetyMetrics

export gate_status, hard_current_limit_status

gate_status() = (experiment="F7",status="NOT_RUN",
    reason="No frozen graph-controller design entered PowerDynamics validation.",
    hard_current_limit="NOT_MODELED")

hard_current_limit_status() = "NOT_MODELED: SimpleGFLDC has no hard current limiter or validated Imax parameter."

end
