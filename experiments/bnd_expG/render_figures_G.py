from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
TAB = ROOT / "reports" / "experiment_G" / "tables"
OUT = ROOT / "reports" / "experiment_G" / "figures"
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"figure.dpi": 140, "savefig.dpi": 180, "font.size": 9,
                     "axes.grid": True, "grid.alpha": 0.25})

def read(name):
    p = TAB / name
    return pd.read_csv(p) if p.exists() and p.stat().st_size else None

def save(fig, name):
    fig.tight_layout()
    fig.savefig(OUT / name, bbox_inches="tight")
    plt.close(fig)

d = read("TABLE_G01_zero_singular_values.csv")
if d is not None:
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for label, g in d.groupby("gain_pattern"):
        ax.semilogy(g["index"], np.maximum(g["singular_value"], 1e-16), marker=".", label=label)
    ax.set(xlabel="Singular-value index", ylabel="Singular value of all-GFL state matrix",
           title="All-GFL zero structure")
    ax.legend()
    save(fig, "FIG_G01_zero_mode_singular_values.png")

d = read("TABLE_G01_near_zero_poles.csv")
if d is not None:
    fig, ax = plt.subplots(figsize=(6.2, 5))
    for label, g in d.groupby("gain_pattern"):
        ax.scatter(g.lambda_real, g.lambda_imag, s=22, label=label)
    ax.axvline(0, color="black", lw=.8)
    ax.axhline(0, color="black", lw=.8)
    ax.set(xlabel="Real part (s⁻¹)", ylabel="Imaginary part (rad/s)",
           title="Near-zero physical poles after gauge deflation")
    ax.legend()
    save(fig, "FIG_G02_near_zero_spectrum.png")

d = read("TABLE_G03_sg_spectral_authority.csv")
if d is not None:
    g=d.sort_values("gamma_per_MW")
    fig, ax=plt.subplots(figsize=(7.4,4.4))
    ax.bar(g.bus.astype(str), g.gamma_per_MW)
    ax.set(xlabel="Generator bus", ylabel="Collective spectral authority per MW",
           title="Retained-SG authority ranking")
    save(fig,"FIG_G03_authority_per_MW.png")

d = read("TABLE_G03_sg_spectral_authority.csv")
if d is not None and {"gamma_direct","gamma_collective"}.issubset(d.columns):
    fig, ax=plt.subplots(figsize=(7.4,4.4))
    x=np.arange(len(d)); w=.38
    ax.bar(x-w/2,d.gamma_direct,w,label="Local diagonal term")
    ax.bar(x+w/2,d.gamma_collective,w,label="Off-diagonal/network residual")
    ax.set_xticks(x,d.bus.astype(str))
    ax.set(xlabel="Generator bus",ylabel="Authority contribution",
           title="Local versus collective authority")
    ax.legend()
    save(fig,"FIG_G04_direct_vs_collective_authority.png")

d=read("TABLE_G07_beta_sweep.csv")
if d is not None:
    fig,ax=plt.subplots(figsize=(6.8,4.3))
    ax.plot(d.beta,d.small_gain_margin,marker="o")
    ax.axhline(0,color="red",lw=1)
    ax.set_xscale("symlog",linthresh=max(d.beta.max()*1e-4,1e-12))
    ax.set(xlabel="Dimensionless full-block radius β",ylabel="Small-gain margin",
           title="Certified robust radius")
    save(fig,"FIG_G05_robust_radius_vs_beta.png")

d=read("TABLE_G07_resolvent_sigma_min.csv")
if d is not None:
    fig,ax=plt.subplots(figsize=(6.8,4.3))
    ax.semilogx(np.maximum(d.omega,1e-5),d.sigma_min)
    ax.set(xlabel="Frequency ω (rad/s)",ylabel="Minimum singular value",
           title="Shifted resolvent singular-value profile")
    save(fig,"FIG_G06_resolvent_sigma_min.png")

d=read("TABLE_G12_unit_disturbance_response.csv")
if d is not None and "frequency_Hz_per_MW" in d:
    fig,ax=plt.subplots(figsize=(7.2,4.2))
    ax.plot(d.time_s,d.frequency_Hz_per_MW)
    ax.set(xlabel="Time (s)",ylabel="Frequency deviation (Hz/MW)",
           title="Linear-model COI frequency response")
    save(fig,"FIG_G07_frequency_response.png")
    fig,ax=plt.subplots(figsize=(7.2,4.2))
    ax.plot(d.time_s,d.rocof_Hz_s_per_MW)
    ax.set(xlabel="Time (s)",ylabel="RoCoF (Hz/s per MW)",
           title="Linear-model COI RoCoF response")
    save(fig,"FIG_G08_rocof_response.png")

d=read("TABLE_G11_modal_transient_residues.csv")
if d is not None and "q_abs" in d:
    g=d.sort_values("q_abs",ascending=False).head(18).sort_values("q_abs")
    fig,ax=plt.subplots(figsize=(7.2,5))
    ax.barh(g["mode"].astype(str),g.q_abs)
    ax.set(xlabel="Absolute RoCoF residue |q_j|",ylabel="Mode index",
           title="Modal RoCoF contributions")
    save(fig,"FIG_G09_modal_rocof_contributions.png")

d=read("TABLE_G06_coordinate_searches.csv")
if d is not None:
    fig,ax=plt.subplots(figsize=(6.8,4.3))
    ax.scatter(d.retained_sg_MW,d.robustified_abscissa,c="tab:blue")
    for _,r in d.iterrows(): ax.annotate(r["order"],(r.retained_sg_MW,r.robustified_abscissa),fontsize=7)
    ax.axhline(-.05,color="red",ls="--",label="Required robust boundary")
    ax.set(xlabel="Retained SG (MW)",ylabel="Robustified spectral abscissa (s⁻¹)",
           title="Exact-spectrum coordinate candidates")
    ax.legend()
    save(fig,"FIG_G10_robust_margin_vs_retained_SG.png")

d=read("TABLE_G12_design_sensitivity_beta_rocof.csv")
if d is not None:
    b=d[d.sweep=="beta"].sort_values("parameter_value")
    fig,ax=plt.subplots(figsize=(6.8,4.3))
    ax.plot(b.parameter_value,b.max_GFL_MW,marker="o")
    ax.set_xscale("log")
    ax.set(xlabel="Dimensionless robust radius β",ylabel="Converted GFL power (MW)",
           title="Maximum replacement versus robust radius")
    save(fig,"FIG_G11_max_GFL_vs_beta.png")
    r=d[d.sweep=="rocof"].sort_values("parameter_value")
    fig,ax=plt.subplots(figsize=(6.8,4.3))
    ax.plot(r.parameter_value,r.max_GFL_MW,marker="o")
    ax.set(xlabel="RoCoF limit (Hz/s)",ylabel="Converted GFL power (MW)",
           title="Maximum replacement versus RoCoF limit")
    save(fig,"FIG_G12_max_GFL_vs_rocof.png")

d=read("TABLE_G13_final_per_generator_design.csv")
if d is not None:
    fig,ax=plt.subplots(figsize=(7.2,4.2))
    ax.bar(d.bus.astype(str),d.rho_converted)
    ax.set(xlabel="Generator bus",ylabel="Converted fraction ρ",ylim=(0,1),
           title="Final synchronous/GFL split")
    save(fig,"FIG_G13_rho_by_bus.png")
    fig,ax=plt.subplots(figsize=(7.2,4.4))
    ax.plot(d.bus,d.Kp,"o-",label="Kp")
    ax.plot(d.bus,d.Ki,"s-",label="Ki")
    ax.set(xlabel="Generator bus",ylabel="PLL gain",title="Independent PLL gains by bus")
    ax.legend()
    save(fig,"FIG_G14_Kp_Ki_by_bus.png")

d=read("TABLE_G15_eigenvalue_comparison.csv")
if d is not None:
    fig,ax=plt.subplots(figsize=(6.8,4.8))
    ax.scatter(d.critical_real,d.critical_imag,s=35,label="Critical pole")
    for _,r in d.iterrows(): ax.annotate(r["case"],(r.critical_real,r.critical_imag),fontsize=7)
    ax.axvline(-.05,color="red",ls="--",label="Nominal requirement")
    ax.axhline(0,color="black",lw=.7)
    ax.set(xlabel="Real part (s⁻¹)",ylabel="Imaginary part (rad/s)",
           title="Nominal and robustified critical poles")
    ax.legend()
    save(fig,"FIG_G15_nominal_robust_eigenvalue_map.png")

d=read("TABLE_G14_blinded_ExpE_comparison.csv")
if d is not None:
    vals={r.metric:r.value for _,r in d.iterrows()}
    fig,ax=plt.subplots(figsize=(6.8,4.3))
    names=["ExpE provisional","ExpG frozen"]
    yy=[]
    for k in ("ExpE provisional retained SG MW","ExpG retained SG MW"):
        try: yy.append(float(vals[k]))
        except Exception: yy.append(np.nan)
    ax.bar(names,yy,color=["0.55","tab:blue"])
    ax.set(ylabel="Retained SG (MW)",title="Blinded ExpE/ExpG retained-capacity comparison")
    save(fig,"FIG_G16_ExpE_vs_ExpG_summary.png")

d = read("TABLE_G16_powerdynamics_spectral_validation.csv")
if d is not None and "alpha_PD_s_inv" in d:
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    x = np.arange(len(d))
    ax.bar(x, d.alpha_PD_s_inv, color=["0.55", "tab:blue"], label="PowerDynamics")
    analytic = d.alpha_ExpG_s_inv.to_numpy(dtype=float)
    keep = np.isfinite(analytic)
    ax.scatter(x[keep], analytic[keep], marker="D", color="black", zorder=4,
               label="ExpG analytical reduced model")
    ax.axhline(-.05, color="red", ls="--", label="Required strict boundary")
    ax.set_xticks(x, d["case"].str.replace("_", " "))
    ax.set(ylabel="Spectral abscissa α (s⁻¹)",
           title="Post-freeze PowerDynamics spectral check")
    ax.legend()
    save(fig, "FIG_G17_powerdynamics_spectral_validation.png")

d = read("TABLE_G17_powerdynamics_tds_trajectories.csv")
if d is not None:
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    for pulse, g in d.groupby("pulse_fraction"):
        ax.plot(g.time_s, g.retained_sg_coi_frequency_deviation_Hz * 1e6,
                label=f"{pulse:g} setpoint fraction")
    ax.axvspan(1.0, 1.1, color="0.8", alpha=.45, label="Load pulse window")
    ax.set(xlabel="Time (s)", ylabel="Retained-SG COI frequency deviation (μHz)",
           title="Nonlinear PowerDynamics response to two load pulses")
    ax.legend()
    save(fig, "FIG_G18_powerdynamics_tds_frequency.png")

    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    for pulse, g in d.groupby("pulse_fraction"):
        ax.plot(g.time_s, g.max_voltage_deviation_pu * 1e6,
                label=f"{pulse:g} setpoint fraction")
    ax.axvspan(1.0, 1.1, color="0.8", alpha=.45, label="Load pulse window")
    ax.set(xlabel="Time (s)", ylabel="Maximum bus-voltage deviation (μpu)",
           title="Nonlinear voltage response to two load pulses")
    ax.legend()
    save(fig, "FIG_G19_powerdynamics_tds_voltage.png")

print("Rendered ExpG figures to", OUT)
