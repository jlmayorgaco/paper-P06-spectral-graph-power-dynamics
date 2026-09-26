using LinearAlgebra

export dominant_modes, spectral_summary

"Return the leading modes ordered by decreasing real part."
function dominant_modes(audit::StabilityAudit; n = 10)
    perm = sortperm(audit.nontrivial_eigenvalues; by = real, rev = true)
    return audit.nontrivial_eigenvalues[perm[1:min(n, length(perm))]]
end

"Convert a stability audit into a compact, serializable summary."
function spectral_summary(audit::StabilityAudit)
    return (
        mode_count = length(audit.eigenvalues),
        gauge_mode_count = length(audit.gauge_eigenvalues),
        nontrivial_mode_count = length(audit.nontrivial_eigenvalues),
        max_real = audit.max_real,
        dynamic_margin = audit.dynamic_margin,
        stable = audit.stable,
        finite = audit.finite,
        dominant_modes = dominant_modes(audit),
    )
end
