module GraphPLLContinuation

export gate_status, require_graph_design

gate_status() = (experiment="F6",status="BLOCKED_UPSTREAM_GATE",
    upstream="F3/F4/F5",reason="No graph-structured gain vector or coefficient path is frozen for continuation.")

require_graph_design() = error("F6 graph replacement continuation is stopped: F2 is NON_MODAL and F4/F5 have no candidate.")

end
