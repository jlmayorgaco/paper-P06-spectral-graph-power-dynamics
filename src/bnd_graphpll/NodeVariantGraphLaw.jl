module NodeVariantGraphLaw

export gate_status, require_gain_targets

gate_status() = (experiment="F4-C",status="BLOCKED_UPSTREAM_GATE",
    upstream="F3/F4-A",reason="There are no frozen ideal modal target gains for the node-local regression.")

require_gain_targets() = error("F4-C is stopped: graph-generated local gains need a frozen F3/F4 target.")

end
