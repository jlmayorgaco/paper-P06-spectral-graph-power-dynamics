module GraphPLLPowerDynamicsValidation

using SHA

export frozen_candidate_hash, require_frozen_candidate, validation_gate_status

frozen_candidate_hash(path::AbstractString) = bytes2hex(sha256(read(path)))

"""Enforce the analytical freeze/hash gate before a validation adapter is called."""
function require_frozen_candidate(path::AbstractString,expected_sha256::AbstractString)
    actual=frozen_candidate_hash(path)
    lowercase(actual)==lowercase(expected_sha256) ||
        error("PowerDynamics validation refused: candidate hash does not match the frozen analytical hash")
    return actual
end

validation_gate_status() = (status="BLOCKED_NO_FROZEN_CANDIDATE",
    reason="F2 stopped the graph spectral program; no F3-F6 candidate was frozen. PowerDynamics is not imported by this guard or the analytical pipeline.")

function validate_frozen_candidate(args...;kwargs...)
    error("F7 validation is not entered: no frozen analytical graph-PLL candidate exists. Pass the hash gate and supply a separate validation adapter when one exists.")
end

end
