"""Frozen GFL11 equations and PowerDynamics injector transcription."""

using LinearAlgebra
using PowerDynamics
using PowerDynamics.Library
using ModelingToolkitBase
using SciCompDSL

using ModelingToolkitBase: t_nounits as t, D_nounits as Dt

const GFL11_OMEGA_B = 2π * 60.0

struct GFL11Parameters
    kp_pll::Float64
    ki_pll::Float64
    tau_p::Float64
    kp_p::Float64
    ki_p::Float64
    kp_q::Float64
    ki_q::Float64
    kp_i::Float64
    ki_i::Float64
    xf::Float64
    rf::Float64
    kp_v::Float64
    ki_v::Float64
    leak::Float64
    omega_b::Float64
end

const GFL11_DEFAULTS = GFL11Parameters(
    53.0, 1400.0, 0.03, 0.20, 8.0, 0.20, 8.0, 0.25, 6.0,
    0.15, 0.01, 2.0, 20.0, 0.05, GFL11_OMEGA_B,
)

"""Return the exact equilibrium state from the frozen initialization rule."""
function gfl11_initialize(v::Float64, a::Float64, p0::Float64, q0::Float64, w::Float64)
    terminal = v * cis(a)
    current = conj((p0 + im*q0) / w / terminal)
    local_current = current * cis(-a)
    id, iq = real(local_current), imag(local_current)
    pref, qref = v * id, -v * iq
    return [a, 0.0, pref, qref, id, -iq, id, iq,
            GFL11_DEFAULTS.rf * id, GFL11_DEFAULTS.rf * iq, qref]
end

"""Evaluate f(x,V) and system-base current/P/Q using explicit dq transforms."""
function gfl11_eval(x::AbstractVector, v::Float64, a::Float64,
                    pref::Float64, qref::Float64, vref::Float64,
                    g::Float64, w::Float64, p::GFL11Parameters=GFL11_DEFAULTS)
    θ, xpll, pf, qf, xp, xq, id, iq, xid, xiq, xv = x
    vd = v*cos(a - θ)
    vq = v*sin(a - θ)
    P = vd*id + vq*iq
    Q = vq*id - vd*iq
    err = g*(vref - v)
    qcmd = p.kp_v*err + xv
    idref = p.kp_p*(pref - pf) + xp
    iqref = -(p.kp_q*(qcmd - qf) + xq)
    ed = vd + p.kp_i*(idref - id) + xid - p.xf*iq
    eq = vq + p.kp_i*(iqref - iq) + xiq + p.xf*id
    f = [
        p.kp_pll*vq + xpll,
        p.ki_pll*vq,
        (P - pf)/p.tau_p,
        (Q - qf)/p.tau_p,
        p.ki_p*(pref - pf),
        p.ki_q*(qcmd - qf),
        (p.omega_b/p.xf)*(ed - vd - p.rf*id + p.xf*iq),
        (p.omega_b/p.xf)*(eq - vq - p.rf*iq - p.xf*id),
        p.ki_i*(idref - id),
        p.ki_i*(iqref - iq),
        p.ki_v*err - p.leak*(xv - qref),
    ]
    I = w*(id + im*iq)*cis(θ)
    return f, real(I), imag(I), w*P, w*Q
end

"""Exact PowerDynamics injector model for the same eleven states."""
@mtkmodel GFL11Injector begin
    @components begin
        terminal = Terminal()
    end
    @parameters begin
        kp_pll = 53.0
        ki_pll = 1400.0
        tau_p = 0.03
        kp_p = 0.20
        ki_p = 8.0
        kp_q = 0.20
        ki_q = 8.0
        kp_i = 0.25
        ki_i = 6.0
        xf = 0.15
        rf = 0.01
        kp_v = 2.0
        ki_v = 20.0
        leak = 0.05
        omega_b = GFL11_OMEGA_B
        g = 0.03625
        w = 1.0
        p_ref = 0.0
        q_ref = 0.0
        v_ref = 1.0
    end
    @variables begin
        theta(t), [guess=0.0]
        x_pll(t), [guess=0.0]
        p_f(t), [guess=0.0]
        q_f(t), [guess=0.0]
        x_p(t), [guess=0.0]
        x_q(t), [guess=0.0]
        i_d(t), [guess=0.0]
        i_q(t), [guess=0.0]
        x_id(t), [guess=0.0]
        x_iq(t), [guess=0.0]
        x_v(t), [guess=0.0]
        v_mag(t)
        v_d(t)
        v_q(t)
        P(t)
        Q(t)
        err(t)
        q_cmd(t)
        id_ref(t)
        iq_ref(t)
        e_d(t)
        e_q(t)
    end
    @equations begin
        v_mag ~ sqrt(terminal.u_r^2 + terminal.u_i^2)
        v_d ~ terminal.u_r*cos(theta) + terminal.u_i*sin(theta)
        v_q ~ terminal.u_i*cos(theta) - terminal.u_r*sin(theta)
        P ~ v_d*i_d + v_q*i_q
        Q ~ v_q*i_d - v_d*i_q
        err ~ g*(v_ref - v_mag)
        q_cmd ~ kp_v*err + x_v
        id_ref ~ kp_p*(p_ref - p_f) + x_p
        iq_ref ~ -(kp_q*(q_cmd - q_f) + x_q)
        e_d ~ v_d + kp_i*(id_ref - i_d) + x_id - xf*i_q
        e_q ~ v_q + kp_i*(iq_ref - i_q) + x_iq + xf*i_d
        Dt(theta) ~ kp_pll*v_q + x_pll
        Dt(x_pll) ~ ki_pll*v_q
        Dt(p_f) ~ (P - p_f)/tau_p
        Dt(q_f) ~ (Q - q_f)/tau_p
        Dt(x_p) ~ ki_p*(p_ref - p_f)
        Dt(x_q) ~ ki_q*(q_cmd - q_f)
        Dt(i_d) ~ (omega_b/xf)*(e_d - v_d - rf*i_d + xf*i_q)
        Dt(i_q) ~ (omega_b/xf)*(e_q - v_q - rf*i_q - xf*i_d)
        Dt(x_id) ~ ki_i*(id_ref - i_d)
        Dt(x_iq) ~ ki_i*(iq_ref - i_q)
        Dt(x_v) ~ ki_v*err - leak*(x_v - q_ref)
        terminal.i_r ~ w*(i_d*cos(theta) - i_q*sin(theta))
        terminal.i_i ~ w*(i_d*sin(theta) + i_q*cos(theta))
    end
end

"""Finite-difference state-space terminal-current transfer matrix."""
function gfl11_transfer(xeq, v, a, pref, qref, vref, g, w, freqs;
                        state_step=1e-7, input_step=1e-7)
    n = length(xeq)
    function fu(x, u)
        vv = hypot(u[1], u[2])
        aa = atan(u[2], u[1])
        f, ir, ii, _, _ = gfl11_eval(x, vv, aa, pref, qref, vref, g, w)
        return f, [ir, ii]
    end
    ueq = [v*cos(a), v*sin(a)]
    A = zeros(n, n); C = zeros(2, n); B = zeros(n, 2); D = zeros(2, 2)
    for j in 1:n
        h = state_step*max(1.0, abs(xeq[j]))
        xp, xm = copy(xeq), copy(xeq); xp[j] += h; xm[j] -= h
        fp, yp = fu(xp, ueq); fm, ym = fu(xm, ueq)
        A[:,j] .= (fp - fm)/(2h); C[:,j] .= (yp - ym)/(2h)
    end
    for j in 1:2
        h = input_step*max(1.0, abs(ueq[j]))
        up, um = copy(ueq), copy(ueq); up[j] += h; um[j] -= h
        fp, yp = fu(xeq, up); fm, ym = fu(xeq, um)
        B[:,j] .= (fp - fm)/(2h); D[:,j] .= (yp - ym)/(2h)
    end
    values = NamedTuple[]
    I_n = Matrix{ComplexF64}(I, n, n)
    for hz in freqs
        s = 2π*hz*im
        H = C*((s*I_n - A)\B) + D
        push!(values, (frequency_hz=hz, sigma_min=minimum(svdvals(s*I_n-A)),
                       condition=cond(s*I_n-A),
                       h11=H[1,1], h12=H[1,2], h21=H[2,1], h22=H[2,2]))
    end
    return A, B, C, D, values
end
