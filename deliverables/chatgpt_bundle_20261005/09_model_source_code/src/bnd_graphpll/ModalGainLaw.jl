module ModalGainLaw

export gate_status, require_modal_gate

gate_status() = (experiment="F3",status="BLOCKED_UPSTREAM_GATE",
    upstream="F2",upstream_status="NON_MODAL",
    reason="The actual GFL/network characteristic does not separate by the F0 graph eigenvalue.")

function require_modal_gate()
    error("F3 is stopped: F2 classified the detailed model NON_MODAL; no Kp*(nu) or Ki*(nu) law is defined.")
end

end
