import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

TABLES, FIGURES, MATRICES = sys.argv[1:4]
os.makedirs(FIGURES, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.22})


def table(name):
    path = os.path.join(TABLES, name)
    return pd.read_csv(path) if os.path.exists(path) else pd.DataFrame()


def save(name, title, draw):
    fig, ax = plt.subplots(figsize=(7.4, 4.5), constrained_layout=True)
    ax.set_title(title)
    try:
        draw(fig, ax)
    except Exception as exc:
        ax.text(.5, .5, f"Unavailable in this run\n{type(exc).__name__}: {exc}",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
    fig.savefig(os.path.join(FIGURES, name), dpi=180)
    plt.close(fig)


def no_data(ax, text="No converged or finite data"):
    ax.text(.5, .5, text, ha="center", va="center", transform=ax.transAxes)


z = table("TABLE_C03_zero_frequency_scaling.csv")
save("FIG_C01_zero_singularity_states.png", "C33 condensed block: smallest singular values", lambda fig, ax: (
    ax.semilogy(np.arange(1, len(z) + 1), np.sort(z.Tcc_sigma_min.values)[::-1], ".-"),
    ax.set_xlabel("zero-frequency sample index"), ax.set_ylabel(r"$\sigma_{min}(sI-A_{cc})$")))
save("FIG_C02_zero_frequency_scaling.png", "Empirical low-frequency self-energy scaling", lambda fig, ax: (
    [ax.loglog(g.s_abs, g.Pi_norm2, "o-", label=d) for d, g in z.groupby("direction")],
    ax.set_xlabel("|s| (s⁻¹)"), ax.set_ylabel(r"$\|\Pi_q(s)\|_2$"), ax.legend()))

b = table("TABLE_C04_backbone_diagnostics.csv")
save("FIG_C03_graph_spectrum.png", "Primary and secondary synchronizing-backbone spectra", lambda fig, ax: (
     [ax.plot(np.arange(len(g)) + 1, g.backbone_eigenvalue, "o-", label=f"{case}, {bb}")
     for (case, bb), g in b.groupby(["case", "backbone"]) if case.startswith("C33")],
    ax.axhline(0, color="black", lw=.8), ax.set_xlabel("mode index"),
    ax.set_yscale("symlog", linthresh=10),
    ax.set_ylabel("synchronizing-backbone eigenvalue (s⁻², symlog scale)"), ax.legend()))
save("FIG_C04_gauge_mode.png", "C33 primary backbone spectrum with retained gauge mode", lambda fig, ax: (
    ax.scatter(np.arange(len(b[(b.case.str.startswith("C33")) & (b.backbone == "primary")])),
               b[(b.case.str.startswith("C33")) & (b.backbone == "primary")].backbone_eigenvalue,
               c=np.where(np.abs(b[(b.case.str.startswith("C33")) & (b.backbone == "primary")].backbone_eigenvalue)<1e-8, "crimson", "steelblue")),
    ax.axhline(0, color="black", lw=.8), ax.set_xlabel("mode index"), ax.set_ylabel("eigenvalue (s⁻²)")))

f = table("TABLE_C06_Psi_frequency_metrics.csv")
save("FIG_C05_psi_component_norms.png", "C33 graph-coordinate self-energy source split", lambda fig, ax: (
    [ax.loglog(g.frequency_hz, g[col], label=label) for col, label in
     [("D0_term_norm", r"$s\hat D_0$"), ("DeltaL_norm", r"$\widehat{\Delta L}$"), ("pi_norm", r"$\hat\Pi_q$"), ("psi_norm", r"$\hat\Psi$" )]
     for g in [f[f.case.str.startswith("C33")]]],
    ax.set_xlabel("frequency (Hz)"), ax.set_ylabel("Frobenius norm"), ax.legend()))

def heatmap_psi(fig, ax):
    p = os.path.join(MATRICES, "C33_Psihat_sample.csv")
    if not os.path.exists(p): return no_data(ax)
    d = pd.read_csv(p).iloc[:, 1:].values
    im = ax.imshow(np.abs(d), origin="lower", aspect="auto", cmap="magma")
    fig.colorbar(im, ax=ax, label=r"$|\hat\Psi_{ij}|$")
    ax.set_xlabel("graph mode j"); ax.set_ylabel("graph mode i")
save("FIG_C06_psihat_heatmap.png", "C33 |Psi-hat| at s = 0.3 + 1.1j s⁻¹", heatmap_psi)

save("FIG_C07_diagonal_modal_terms.png", "C33 modal coupling magnitude across frequency", lambda fig, ax: (
    ax.loglog(f[f.case.str.startswith("C33")].frequency_hz,
              f[f.case.str.startswith("C33")].psi_norm, label="total"),
    ax.loglog(f[f.case.str.startswith("C33")].frequency_hz,
              f[f.case.str.startswith("C33")].offdiagonal_norm, label="off-diagonal"),
    ax.set_xlabel("frequency (Hz)"), ax.set_ylabel("Frobenius norm"), ax.legend()))
save("FIG_C08_DG_eff_spectrum.png", "C33 generalized dissipative operator spectrum", lambda fig, ax: (
    ax.semilogx(f[f.case.str.startswith("C33")].frequency_hz,
                f[f.case.str.startswith("C33")].DG_eff_min, label="minimum eigenvalue"),
    ax.semilogx(f[f.case.str.startswith("C33")].frequency_hz,
                f[f.case.str.startswith("C33")].DG_eff_max, label="maximum eigenvalue"),
    ax.axhline(0, color="black", lw=.8), ax.set_xlabel("frequency (Hz)"),
    ax.set_ylabel(r"eigenvalue of $D_G^{eff}$ (s⁻¹)"), ax.legend()))
save("FIG_C09_offdiagonal_ratio.png", "C33 off-diagonal fraction of Psi-hat", lambda fig, ax: (
    ax.semilogx(f[f.case.str.startswith("C33")].frequency_hz,
                f[f.case.str.startswith("C33")].offdiag_ratio),
    ax.set_xlabel("frequency (Hz)"), ax.set_ylabel("off-diagonal Frobenius ratio")))
save("FIG_C10_commutator_vs_frequency.png", "C33 graph/dynamic basis mismatch", lambda fig, ax: (
    ax.loglog(f[f.case.str.startswith("C33")].frequency_hz,
              f[f.case.str.startswith("C33")].chi_comm),
    ax.set_xlabel("frequency (Hz)"), ax.set_ylabel(r"$\chi_{comm,\Psi}$")))

poles = table("TABLE_C07_physical_pole_graph_mapping.csv")
save("FIG_C11_physical_pole_graph_projection.png", "C33 physical poles: graph-projection concentration", lambda fig, ax: (
    ax.scatter(poles[poles.case.str.startswith("C33")].frequency_hz,
               poles[poles.case.str.startswith("C33")].graph_concentration,
               c=poles[poles.case.str.startswith("C33")].dominant_graph_mode, cmap="viridis"),
    ax.set_xlabel("frequency (Hz)"), ax.set_ylabel("largest graph projection fraction")))

g = table("TABLE_C08_Gamma_exact.csv")
save("FIG_C12_exact_Gamma_vs_frequency.png", "C33 exact intermodal self-energy at analyzed poles", lambda fig, ax: (
    ax.scatter(g[g.case.str.startswith("C33")].evaluation_s_imag.abs() / (2*np.pi),
               g[g.case.str.startswith("C33")].gamma_abs,
               c=g[g.case.str.startswith("C33")].graph_mode_k, cmap="plasma"),
    ax.set_xlabel("frequency (Hz)"), ax.set_ylabel(r"$|\Gamma_k|$ (s⁻²)")))

path = table("TABLE_C11_exact_pathways.csv")
def pathway_heatmap(fig, ax):
    d = path[path.case.str.startswith("C33")]
    if d.empty: return no_data(ax)
    k = d.mode_k.iloc[0]
    d = d[d.mode_k == k]
    pivot = d.pivot_table(index="mode_l", columns="mode_m", values="contribution_abs", aggfunc="sum", fill_value=0)
    im=ax.imshow(pivot.values, origin="lower", aspect="auto", cmap="magma")
    fig.colorbar(im, ax=ax, label="exact pathway magnitude (s⁻²)")
    ax.set_xlabel("complement mode m"); ax.set_ylabel("complement mode l")
save("FIG_C13_exact_pathway_heatmap.png", "C33 exact pathway magnitude matrix", pathway_heatmap)

susc = table("TABLE_C12_pairwise_susceptibility.csv")
save("FIG_C14_intermodal_susceptibility_ranking.png", "C33 approximate pairwise susceptibility", lambda fig, ax: (
    ax.bar(susc[susc.case.str.startswith("C33")].mode_l.astype(str),
           susc[susc.case.str.startswith("C33")].susceptibility_abs),
    ax.set_xlabel("complementary graph mode"), ax.set_ylabel(r"$|\mathcal{S}_{k\leftarrow\ell}|$ (s⁻²)")))
save("FIG_C15_coupling_vs_susceptibility_vs_exact_pathway.png", "C33 coupling, detuning-adjusted susceptibility, exact pathway", lambda fig, ax: (
    ax.scatter(susc[susc.case.str.startswith("C33")].coupling_product_abs,
               susc[susc.case.str.startswith("C33")].exact_diagonal_path_abs, label="coupling product", marker="o"),
    ax.scatter(susc[susc.case.str.startswith("C33")].susceptibility_abs,
               susc[susc.case.str.startswith("C33")].exact_diagonal_path_abs, label="susceptibility", marker="x"),
    ax.set_xscale("log"), ax.set_yscale("log"), ax.set_xlabel("proxy magnitude (s⁻⁴ or scaled)"),
    ax.set_ylabel("exact diagonal pathway magnitude (s⁻²)"), ax.legend()))
save("FIG_C16_real_resonance_map.png", "C33 dynamic detuning versus exact pathway contribution", lambda fig, ax: (
    ax.scatter(susc[susc.case.str.startswith("C33")].dynamic_detuning_abs,
               susc[susc.case.str.startswith("C33")].exact_diagonal_path_abs,
               c=susc[susc.case.str.startswith("C33")].susceptibility_abs, cmap="viridis"),
    ax.set_xscale("log"), ax.set_yscale("log"), ax.set_xlabel(r"$|t_{\ell}^{0}(s_k^0)|$ (s⁻²)"),
    ax.set_ylabel("exact diagonal pathway (s⁻²)")))

pred = table("TABLE_C10_pole_predictor.csv")
def poles_plot(fig, ax):
    d=pred[pred.case.str.startswith("C33")]
    if d.empty: return no_data(ax)
    ax.scatter(d.s0_real,d.s0_imag/(2*np.pi),label="uncoupled root",marker="o")
    ax.scatter(d.pred_real,d.pred_imag/(2*np.pi),label="predicted",marker="x")
    ax.scatter(d.exact_real,d.exact_imag/(2*np.pi),label="PowerDynamics",marker="+")
    ax.set_xlabel("real part (s⁻¹)"); ax.set_ylabel("imaginary part / 2π (Hz)"); ax.legend()
save("FIG_C17_exact_vs_predicted_poles.png", "C33 diagonal roots, first-order prediction, actual poles", poles_plot)
save("FIG_C18_predictor_error_vs_eta.png", "C33 pole predictor error versus relative shift", lambda fig, ax: (
    ax.loglog(pred[pred.case.str.startswith("C33")].eta_shift,
              pred[pred.case.str.startswith("C33")].rel_error, "o"),
    ax.set_xlabel(r"$|\Delta s|/|s_0|$"), ax.set_ylabel("relative pole error")))
save("FIG_C19_predictor_error_vs_Trr_condition.png", "C33 predictor error versus complement conditioning", lambda fig, ax: (
    ax.loglog(1/np.maximum(pred[pred.case.str.startswith("C33")].Trr_sigma_min,1e-300),
              pred[pred.case.str.startswith("C33")].rel_error, "o"),
    ax.set_xlabel(r"$1/\sigma_{min}(T_{rr})$"), ax.set_ylabel("relative pole error")))

comp=table("TABLE_C13_baseline_vs_replacement.csv")
save("FIG_C20_baseline_vs_bus33_mechanism.png", "All-SG versus bus-33 GFL matched pole shifts", lambda fig, ax: (
    ax.scatter(comp[comp.replacement_bus==33].replacement_graph_mode,
               comp[comp.replacement_bus==33].delta_real, label="real shift"),
    ax.set_xlabel("matched graph mode"), ax.set_ylabel("replacement minus baseline real part (s⁻¹)"), ax.legend()))
cross=table("TABLE_C15_cross_bus_summary.csv")
save("FIG_C21_cross_bus_summary.png", "Cross-case graph and intermodal identity status", lambda fig, ax: (
    ax.bar(cross.case, cross.number_analyzed_modes, color=np.where(cross.graph_identity_pass, "steelblue", "firebrick")),
    ax.set_xlabel("case"), ax.set_ylabel("analyzed non-conjugate poles")))

abl=table("TABLE_C17_component_ablation.csv")
save("FIG_C22_component_ablation.png", "Operator-only Psi component ablation (not a realizable plant)", lambda fig, ax: (
    [ax.loglog(g.frequency_hz,g.gamma_abs,"o-",label=op) for op,g in
     abl[abl.case.str.startswith("C33")].groupby("operator")],
    ax.set_xlabel("frequency (Hz)"),ax.set_ylabel("diagnostic |Gamma| (s⁻²)"),ax.legend()))

sens=table("TABLE_C14_backbone_sensitivity.csv")
save("FIG_C23_backbone_sensitivity.png", "C33 primary versus secondary backbone sensitivities", lambda fig, ax: (
    ax.scatter(sens[sens.case.str.startswith("C33")].gamma_primary,
               sens[sens.case.str.startswith("C33")].gamma_secondary,
               c=sens[sens.case.str.startswith("C33")].subspace_overlap,cmap="viridis"),
    ax.set_xscale("log"),ax.set_yscale("log"),ax.set_xlabel("primary |Gamma| (s⁻²)"),
    ax.set_ylabel("secondary |Gamma| (s⁻²)")))
