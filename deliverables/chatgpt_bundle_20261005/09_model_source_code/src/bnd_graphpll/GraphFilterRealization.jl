module GraphFilterRealization

export gate_status, require_gain_law

gate_status() = (experiment="F4-A/B",status="BLOCKED_UPSTREAM_GATE",
    upstream="F3",reason="No admissible spectral gain law exists after the F2 NON_MODAL stop.")

require_gain_law() = error("F4 spectral controller realization is stopped: F3 produced no hp(nu) or hi(nu).")

end
