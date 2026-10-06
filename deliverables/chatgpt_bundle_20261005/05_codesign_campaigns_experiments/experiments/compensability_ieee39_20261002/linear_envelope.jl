using LinearAlgebra, CSV, DataFrames, TOML, SHA
BLAS.set_num_threads(1)

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "src", "bnd_model_expN", "PDExactDesignN.jl"))
using .PDExactDesignN
const REPORTS = joinpath(ROOT, "reports", "nonlinear_codesign_20261001")
const OUT = joinpath(ROOT, "reports", "poster", "ias2026", "compensability_ieee39_20261002")
const CANDIDATE = joinpath(REPORTS, "candidate_final_physical.toml")
const JAC_PATH = joinpath(REPORTS, "stationarity_final", "0001_grad_jacobian.csv")
const G_PATH = joinpath(REPORTS, "stationarity_final", "0001_grad_constraints.toml")
mkpath(OUT)

# Feasible-start strictly convex QP for the local trust-box envelope.
function feasible_qp(H, c, A, b, z)
    W = Int[]
    n = length(c)
    for _ in 1:500
        if isempty(W)
            d = -(H \ (H*z + c))
            mu = Float64[]
        else
            AW = A[W, :]
            Z = nullspace(AW; rtol=1e-12)
            d = -Z * ((Z' * H * Z) \ (Z' * (H*z + c)))
            mu = pinv(AW'; rtol=1e-12) * (-(H*(z+d)+c))
        end
        if norm(d, Inf) < 2e-10
            if isempty(W) || minimum(mu) >= -1e-8
                return (; z, mu, W, status="OK")
            end
            deleteat!(W, argmin(mu))
            continue
        end
        alpha = 1.0
        blocker = 0
        for i in eachindex(b)
            i in W && continue
            den = dot(A[i, :], d)
            den <= 1e-10 && continue
            t = (b[i] - dot(A[i, :], z)) / den
            if t < alpha
                alpha = max(0.0, t)
                blocker = i
            end
        end
        z += alpha*d
        if blocker != 0
            (isempty(W) || rank(A[vcat(W, blocker), :]; rtol=1e-10) > length(W)) ||
                return (; z, mu, W, status="DEPENDENT_ACTIVE_ROWS")
            push!(W, blocker)
        end
    end
    (; z, mu=Float64[], W, status="ITERATION_LIMIT")
end

function main()
    cand = TOML.parsefile(CANDIDATE)
    g = Float64.(TOML.parsefile(G_PATH)["constraints"])
    J = Matrix(CSV.read(JAC_PATH, DataFrame))
    size(J) == (length(g), 30) || error("Unexpected tangent shape $(size(J))")

    # Continue along uniform rho increase while allowing every PLL gain to
    # retune inside a finite +/-10% log-gain box. The six event rows and modal
    # bundle come from the frozen 60-s trajectory tangent at the same design.
    jgain = J[:, 11:30]
    jrho = J[:, 1:10]
    rho0 = Float64.(cand["rho"])
    kp0 = Float64.(cand["Kp"])
    ki0 = Float64.(cand["Ki"])

    # Audit the active rightmost-pole row against fresh finite differences of
    # the full descriptor spectrum at the frozen candidate.
    ctx = PDExactDesignN.design_context(ROOT)
    sp = PDExactDesignN.spectrum(ctx, rho0, kp0, ki0)
    vals = eigen(sp.quotient' * sp.model.Ared * sp.quotient).values
    modal_ids = sort([j for j in eachindex(vals) if imag(vals[j]) >= -1e-8];
        by=j -> real(vals[j]), rev=true)[1:6]
    ds = PDExactDesignN.simple_mode_sensitivities(ctx, rho0, kp0, ki0; mode=modal_ids[1])
    j_modal_analytic = real.(vcat(ds.rho, ds.Kp .* kp0, ds.Ki .* ki0)) ./ 0.05
    alpha_at(r, p, q) = PDExactDesignN.spectrum(ctx, r, p, q).alpha
    fd_modal = zeros(30)
    hfd = 2e-5
    for j in 1:30
        rp, rm = copy(rho0), copy(rho0)
        pp, pm = copy(kp0), copy(kp0)
        qp, qm = copy(ki0), copy(ki0)
        if j <= 10
            rp[j] += hfd; rm[j] -= hfd
        elseif j <= 20
            pp[j-10] *= exp(hfd); pm[j-10] *= exp(-hfd)
        else
            qp[j-20] *= exp(hfd); qm[j-20] *= exp(-hfd)
        end
        fd_modal[j] = (alpha_at(rp,pp,qp)-alpha_at(rm,pm,qm))/(2hfd*0.05)
    end
    modal_fd_relative_error = norm(fd_modal-j_modal_analytic, Inf) /
        max(norm(fd_modal, Inf), norm(j_modal_analytic, Inf), 1e-12)
    modal_stored_relative_error = norm(J[1,:]-j_modal_analytic, Inf) /
        max(norm(j_modal_analytic, Inf), 1e-12)
    gain_radius = log(1.1)
    tau_radius = 0.01
    scales = vcat(fill(gain_radius, 20), tau_radius)
    Acore = hcat(jgain .* gain_radius, (jrho * ones(10)) .* tau_radius)
    n = 21
    A = vcat(Acore, Matrix{Float64}(I, n, n), -Matrix{Float64}(I, n, n))
    b = vcat(-g, ones(n), vcat(ones(20), 0.0))
    # Lexicographically favor more replacement, then smaller normalized gain
    # movement. This is a regularized local LP surrogate, not the nonlinear
    # optimum and not a global certificate.
    H = Matrix{Float64}(I, n, n) .* 1e-9
    c = zeros(n)
    c[end] = -tau_radius
    sol = feasible_qp(H, c, A, b, zeros(n))
    sol.status == "OK" || error("Local envelope QP failed: $(sol.status)")
    maximum(A*sol.z-b) < 1e-7 || error("Local envelope violates a linearized row")

    dk = scales[1:20] .* sol.z[1:20]
    tau = tau_radius * sol.z[end]
    delta = vcat(fill(tau, 10), dk)
    predicted = g + J*delta
    original_capacity = cand["replacement_percent"]
    # A uniform rho increment changes the dispatch-weighted replacement by the
    # same number of percentage points, irrespective of individual dispatch.
    p0 = rho0
    capacity_gain_pp = 100 * tau
    new_rho = p0 .+ tau
    all(new_rho .< 1.0) || error("Uniform direction left the interior support")

    result = Dict(
        "status" => "LOCAL_LINEAR_TANGENT_DIAGNOSTIC_ONLY",
        "base_replacement_percent" => original_capacity,
        "direction" => "uniform increase in all ten rho values",
        "gain_box_relative" => 0.10,
        "gain_box_coordinate" => "symmetric in log gain, radius log(1.1), so multiplier interval [1/1.1, 1.1]",
        "max_uniform_delta_rho" => tau,
        "capacity_gain_percentage_points" => capacity_gain_pp,
        "max_predicted_constraint_after_step" => maximum(predicted),
        "predicted_constraints" => predicted,
        "gain_log_step" => dk,
        "active_working_rows" => sol.W,
        "active_multipliers" => sol.mu,
        "rows" => length(g),
        "active_modal_gradient_finite_difference_relative_error" => modal_fd_relative_error,
        "stored_modal_row_vs_fresh_analytic_relative_error" => modal_stored_relative_error,
        "modal_gradient_fd_log_step" => hfd,
        "source_sha256" => Dict(
            "candidate" => bytes2hex(sha256(read(CANDIDATE))),
            "jacobian" => bytes2hex(sha256(read(JAC_PATH))),
            "constraints" => bytes2hex(sha256(read(G_PATH)))),
        "limitations" => [
            "first-order trajectory and modal tangent only",
            "curvature remainder not bounded",
            "finite six-event suite only",
            "does not exclude nonlinear recourse outside this local box or other rho directions",
            "does not certify global capacity optimality"
        ])
    open(joinpath(OUT, "linear_envelope.toml"), "w") do io
        TOML.print(io, result)
    end
    println("LOCAL_ENVELOPE tau=", tau, " capacity_pp=", capacity_gain_pp,
        " predicted_max_g=", maximum(predicted), " active_rows=", sol.W,
        " modal_fd_error=", modal_fd_relative_error,
        " stored_modal_error=", modal_stored_relative_error)
end

main()
