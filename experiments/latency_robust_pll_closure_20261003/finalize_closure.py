"""Build only evidence-backed closure tables and figures from executed runs."""
from __future__ import annotations

import csv
import json
import math
import re
import tomllib
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent
OLD=HERE.parent/"latency_robust_pll_codesign_20261003"

def rows(path:Path):
    with path.open(newline="") as f:return list(csv.DictReader(f))

def write(name:str,records:list[dict]):
    if not records:return
    with (HERE/name).open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)

def design_data(cid:str):
    spec=HERE/"evaluations"/cid/"SPECTRAL.csv"
    event=HERE/"event_validations"/cid/"Q0_RESULT.toml"
    if not spec.exists() or not event.exists():return None
    s=rows(spec)[0];e=tomllib.loads(event.read_text())
    roots=HERE/"evaluations"/cid/"ROOTS.csv"
    rootrows=rows(roots) if roots.exists() else []
    critical=min((float(r["local_crossing_ms"]) for r in rootrows),default=math.nan)
    coverage_file=HERE/"evaluations"/cid/"ROOT_COVERAGE.csv"
    if not coverage_file.exists():coverage_file=HERE/"evaluations"/cid/"ROOT_COVERAGE_HIGH.csv"
    coverage=rows(coverage_file)[0]["status"] if coverage_file.exists() else "MISSING"
    extra=HERE/"evaluations"/cid/"ADDITIONAL_COUNTS.csv"
    additional=rows(extra) if extra.exists() else []
    boundary=(s["complete_contour_bracket_pass"].lower()=="true" or
              (any(r["status"]=="SAFE" and float(r["tau_ms"])<critical for r in additional)
               and any(r["status"]=="UNSAFE" and float(r["tau_ms"])>critical for r in additional)))
    validated=bool(e["all_five_events_pass"]) and math.isfinite(critical) and boundary and coverage=="COMPLETE_NUMERICAL_ROOT_COVERAGE_AT_UNSAFE_CONTOUR"
    return dict(design_id=cid,tau_crit_ms=critical,
        alpha_zero_delay_s_inv=float(s["alpha_zero_delay_s_inv"]),
        actuator_slack=min(float(x["min_SG_actuator_fraction_slack"]) for x in e["events"]),
        all_five_events_pass=bool(e["all_five_events_pass"]),
        max_frequency_deviation_hz=max(float(x["F_peak_Hz"]) for x in e["events"]),
        max_rocof_hz_s=max(float(x["RoCoF_peak_Hz_s"]) for x in e["events"]),
        root_coverage=coverage,boundary_contour_pass=boundary,
        status="VALIDATED_SPECTRUM_AND_FIVE_ZERO_DELAY_EVENTS" if validated else "NOT_VALIDATED")

def main():
    parent=rows(HERE/"F00_PARENT_REPRODUCTION.csv")
    previous=parent[-1]
    candidate_ids=["pending_41ms","pending_next_lp","full20_next_lp"]
    candidate_ids += sorted(p.stem for p in (HERE/"designs").glob("full20_step*.toml"))
    candidate_ids += sorted(p.stem for p in (HERE/"designs").glob("reserve*.toml"))
    candidates=[design_data(cid) for cid in candidate_ids]
    candidates=[x for x in candidates if x]
    write("TABLE_F2_FINAL_OPTIMIZATION_TRACE.csv",candidates)
    feasible=[x for x in candidates if x["status"]=="VALIDATED_SPECTRUM_AND_FIVE_ZERO_DELAY_EVENTS"]
    best=max(feasible,key=lambda x:x["tau_crit_ms"]) if feasible else None
    if not best:raise RuntimeError("no closure design passed five events")
    cid=best["design_id"]
    d=tomllib.loads((HERE/"designs"/(cid+".toml")).read_text())
    z=tomllib.loads((HERE/"designs/Z_zero_delay_tuned.toml").read_text())
    gains=[dict(design_id=cid,bus=30+i,rho=d["rho"][i],Kp=d["Kp"][i],Ki=d["Ki"][i],
                Kp_over_Z=d["Kp"][i]/z["Kp"][i],Ki_over_Z=d["Ki"][i]/z["Ki"][i]) for i in range(10)]
    write("TABLE_F4_FINAL_OPTIMAL_GAINS.csv",gains)
    audit=[]
    for sid,kind,value,base_design in (("uncertainty_Kp33_minus1pct","gain_Kp33_factor",0.99,"pending_41ms"),
                           ("uncertainty_event_plus1MW","event_bus16_plus_MW",101.0,"pending_41ms"),
                           ("uncertainty_final_Kp33_minus1pct","gain_Kp33_factor",0.99,"full20_step6_multi"),
                           ("uncertainty_final_event_plus1MW","event_bus16_plus_MW",101.0,"full20_step6_multi"),
                           ("uncertainty_best_event_plus1MW","event_bus16_plus_MW",101.0,"full20_step11_maxtrust"),
                           ("uncertainty_reserve_event_plus1MW","event_bus16_plus_MW",101.0,"reserve_near_best_20e6"),
                           ("uncertainty_tiny_event_plus1MW","event_bus16_plus_MW",101.0,"full20_step12_tiny"),
                           ("uncertainty_tiny_Kp33_minus1pct","gain_Kp33_factor",0.99,"full20_step12_tiny"),
                           ("uncertainty_step13_event_plus1MW","event_bus16_plus_MW",101.0,"full20_step13_tiny"),
                           ("uncertainty_step13_event_plus0p01MW","event_bus16_plus_MW",100.01,"full20_step13_tiny"),
                           ("uncertainty_step13_event_plus0p02MW","event_bus16_plus_MW",100.02,"full20_step13_tiny")):
        file=HERE/"event_screens"/sid/"Q0_EVENT_METRICS.csv"
        if file.exists():
            r=rows(file)[0]
            audit.append(dict(case_id=sid,kind=kind,value=value,base_design_id=base_design,complete=r["complete"],
                              retcode=r["retcode"],actuator_slack=r["min_SG_actuator_fraction_slack"],
                              pass_screen=r["pass"],scope="ONE_EVENT_ONLY"))
    for sid,value in (("load_plus1pct",1.01),("load_minus1pct",0.99)):
        file=HERE/"load_uncertainty"/sid/"Q0_EVENT_METRICS.csv"
        if file.exists():
            r=rows(file)[0]
            audit.append(dict(case_id=sid,kind="all_ZIP_load_factor",value=value,base_design_id="pending_41ms",complete=r["complete"],
                              retcode=r["retcode"],actuator_slack=r["min_SG_actuator_fraction_slack"],
                              pass_screen=r["pass"],scope="ONE_EVENT_ONLY_REEQUILIBRATED"))
    write("TABLE_F9_UNCERTAINTY_AUDIT.csv",audit)
    kkt=HERE/f"KKT_LP_{cid}.csv"
    if kkt.exists():write("TABLE_F10_KKT_ACTIVE_SET.csv",rows(kkt))

    plt.rcParams.update({"figure.dpi":150,"savefig.dpi":220,"font.size":10,
                          "axes.spines.top":False,"axes.spines.right":False})
    navy="#14213d";teal="#087e8b";orange="#e77c40"
    short={"previous_best":"old best","pending_41ms":"41 ms candidate","pending_next_lp":"sparse LP",
           "full20_next_lp":"full-20 step 1","full20_step2_lp":"full-20 step 2",
           "full20_step3_lp":"full-20 step 3","full20_step4_lp":"full-20 step 4",
           "full20_step5_small":"step 5 single-mode","full20_step5_wide":"step 5 wide",
           "full20_step5_multi":"step 5 minimax",
           "reserve20e6_lp":"early reserve","reserve_final_20e6":"first reserve",
           "reserve_near_best_20e6":"reserve attempt"}
    def label(design_id):
        if design_id in short:return short[design_id]
        m=re.match(r"full20_step(\d+)_(.*)",design_id)
        return f"step {m.group(1)} {m.group(2)}" if m else design_id
    trajectory=sorted((x for x in candidates
                       if x["status"]=="VALIDATED_SPECTRUM_AND_FIVE_ZERO_DELAY_EVENTS"),
                      key=lambda x:x["tau_crit_ms"])
    milestones={"pending_41ms","full20_step2_lp","full20_step5_multi",
                "full20_step8_wide","full20_step11_maxtrust","reserve_near_best_20e6",cid}
    plot_trajectory=[x for x in trajectory if x["design_id"] in milestones]
    labels=["Z","N","old best"]+[label(x["design_id"]) for x in plot_trajectory]
    vals=[float(parent[0]["tau_crit_ms"]),float(parent[1]["tau_crit_ms"]),
          float(parent[2]["tau_crit_ms"])]+[x["tau_crit_ms"] for x in plot_trajectory]
    plt.figure(figsize=(9,4.5));plt.plot(range(len(vals)),vals,"o-",color=teal,lw=2)
    for i,v in enumerate(vals):plt.annotate(f"{v:.3f}",(i,v),xytext=(0,8),textcoords="offset points",ha="center")
    plt.xticks(range(len(labels)),labels,rotation=20,ha="right");plt.ylabel("Critical uniform delay (ms)")
    plt.title("Validated milestones at fixed 88.455% GFL replacement")
    plt.tight_layout();plt.savefig(HERE/"FIG_F1_LATENCY_MARGIN_OPTIMIZATION.png");plt.close()

    modal=[]
    previous_roots=rows(OLD/"evaluations"/json.loads((OLD/"RESULT_SUMMARY.json").read_text())["best_fully_validated_found_id"]/"ROOTS.csv")
    for r in previous_roots:modal.append((0,float(r["local_crossing_ms"]),float(r["frequency_hz"])))
    for i,x in enumerate(plot_trajectory,1):
        path=HERE/"evaluations"/x["design_id"]/"ROOTS_EXPANDED.csv"
        if not path.exists():path=HERE/"evaluations"/x["design_id"]/"ROOTS.csv"
        for r in rows(path):modal.append((i,float(r["local_crossing_ms"]),float(r["frequency_hz"])))
    fig,(ax0,ax1)=plt.subplots(2,1,figsize=(9,6.5),sharex=True,
                              gridspec_kw={"height_ratios":[1,1.25]})
    positions=range(1+len(plot_trajectory))
    floors={i:min(x[1] for x in modal if x[0]==i) for i in positions}
    ax0.plot(list(positions),[floors[i] for i in positions],"o-",color=teal,lw=2)
    ax0.set_ylabel("First crossing (ms)")
    ax0.set_title("Delayed-mode competition along validated designs")
    scatter=ax1.scatter([x[0] for x in modal],
                        [1000*(x[1]-floors[x[0]]) for x in modal],
                        c=[x[2] for x in modal],cmap="viridis",s=65,zorder=3)
    ax1.set_yscale("symlog",linthresh=1)
    ax1.set_ylabel("Gap to first crossing (µs)")
    ax1.grid(axis="y",alpha=.2)
    ax1.set_xticks(list(positions),["old best"]+[label(x["design_id"]) for x in plot_trajectory],rotation=35,ha="right")
    fig.colorbar(scatter,ax=[ax0,ax1],label="Mode frequency (Hz)",shrink=.8)
    fig.savefig(HERE/"FIG_F2_MODAL_FAMILY_EQUALIZATION.png",bbox_inches="tight");plt.close(fig)

    allpoints=[dict(design_id=r["label"],tau_crit_ms=float(r["tau_crit_ms"]),
                    alpha_zero_delay_s_inv=float(r["alpha_zero_delay_s_inv"]),
                    actuator_slack=float(r["minimum_actuator_slack"]),
                    status="VALIDATED_SPECTRUM_AND_FIVE_ZERO_DELAY_EVENTS") for r in parent]+candidates
    unvalidated=[x for x in allpoints if x["status"]!="VALIDATED_SPECTRUM_AND_FIVE_ZERO_DELAY_EVENTS"]
    accepted=[x for x in allpoints if x["status"]=="VALIDATED_SPECTRUM_AND_FIVE_ZERO_DELAY_EVENTS"]
    featured={"Z","N","previous_best",cid}
    xz=-float(parent[0]["alpha_zero_delay_s_inv"])
    xnom=lambda x:1e6*(-x["alpha_zero_delay_s_inv"]-xz)
    plt.figure(figsize=(7.5,4.5))
    plt.scatter([xnom(x) for x in unvalidated],
        [x["tau_crit_ms"] for x in unvalidated],color="#a2aab3",s=35,label="not five-event validated")
    plt.scatter([xnom(x) for x in accepted],
        [x["tau_crit_ms"] for x in accepted],color=teal,s=55,label="five-event validated")
    for x in allpoints:
        if x["design_id"] in featured:
            plt.annotate(label(x["design_id"]),(xnom(x),x["tau_crit_ms"]),fontsize=8)
    plt.xlabel("Zero-delay margin change vs Z (10⁻⁶ s⁻¹)");plt.ylabel("Critical delay (ms)")
    plt.title("Observed designs only — optimized Pareto frontier not established")
    plt.legend();plt.tight_layout();plt.savefig(HERE/"FIG_F3_NOMINAL_LATENCY_PARETO.png");plt.close()

    fig,axes=plt.subplots(1,2,figsize=(10,4.5),sharey=True)
    for ax in axes:
        ax.scatter([x["actuator_slack"]-0.002 for x in unvalidated],
            [x["tau_crit_ms"] for x in unvalidated],color="#a2aab3",s=35,label="not five-event validated")
        ax.scatter([x["actuator_slack"]-0.002 for x in accepted],
            [x["tau_crit_ms"] for x in accepted],color=orange,s=55,label="five-event validated")
        ax.axvline(0,color=navy,ls="--",lw=1)
        ax.set_xlabel("Actuator reserve above frozen limit")
    axes[0].set_ylabel("Critical delay (ms)")
    axes[0].set_title("All observed designs")
    axes[1].set_xlim(-2e-6,4e-5)
    axes[1].set_title("Near frozen actuator limit")
    axes[1].ticklabel_format(axis="x",style="sci",scilimits=(-3,3))
    axes[0].legend(loc="lower right",fontsize=8)
    fig.suptitle("Observed actuator–latency tradeoff")
    fig.tight_layout();fig.savefig(HERE/"FIG_F4_ACTUATOR_LATENCY_FRONTIER.png");plt.close(fig)

    trace=HERE/f"TIME_DOMAIN_V2_{cid}.csv"
    if trace.exists():
        rr=rows(trace);plt.figure(figsize=(7.5,4.5))
        for factor in (0.9,0.98,1.02):
            take=[r for r in rr if abs(float(r["factor"])-factor)<1e-9]
            plt.semilogy([float(r["t_s"]) for r in take],[float(r["relative_amplitude"]) for r in take],label=f"{factor:.0%} τcrit")
        plt.xlabel("Time (s)");plt.ylabel("Relative complex modal amplitude")
        plt.title("Exact linear DDE method of steps (V2; not nonlinear validation)")
        plt.legend();plt.tight_layout();plt.savefig(HERE/"FIG_F5_DELAYED_TIME_DOMAIN_VALIDATION.png");plt.close()
    buses=np.arange(30,40);plt.figure(figsize=(8.5,4.5))
    plt.plot(buses,[x["Kp_over_Z"] for x in gains],"o-",label="Kp / Kp(Z)")
    plt.plot(buses,[x["Ki_over_Z"] for x in gains],"s-",label="Ki / Ki(Z)")
    plt.axhline(1,color=navy,ls="--",lw=1);plt.xticks(buses);plt.xlabel("Generator bus")
    plt.ylabel("Gain ratio");plt.title("Final validated gain retuning at unchanged replacement")
    plt.legend();plt.tight_layout();plt.savefig(HERE/"FIG_F6_GAIN_CHANGES_BY_BUS.png");plt.close()
    (HERE/"BEST_VALIDATED_DESIGN_ID.txt").write_text(cid+"\n")
    print("best",cid,best["tau_crit_ms"],"figures rendered")

if __name__=="__main__":main()
