"""Freeze the targeted latency experiment before new evaluations."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, subprocess, sys, platform, tomllib

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MEGA = ROOT / "experiments/analytical_delay_codesign_mega_20261002"
paths = {
    "design_A": MEGA / "seed_uniform_875.toml",
    "design_B": MEGA / "M1_ZERO_DELAY_DESIGN.toml",
    "trace_count_source": MEGA / "m3_a_trace_integral.jl",
    "root_refinement_source": MEGA / "m3_c_refine_roots.jl",
    "delay_characteristic_source": ROOT / "experiments/delay_dressed_replacement_frontier_20261002/DelayCharacteristic.jl",
    "action_space_source": MEGA / "q56_validate_action_space.jl",
}
source = {name: {"path":p.relative_to(ROOT).as_posix(),"sha256":hashlib.sha256(p.read_bytes()).hexdigest()}
          for name,p in paths.items()}
protocol = {
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "git_head": subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
    "git_dirty_at_start": bool(subprocess.check_output(["git","status","--porcelain"],cwd=ROOT,text=True).strip()),
    "python": platform.python_version(),
    "status": "FROZEN_BEFORE_TARGETED_EVALUATION",
    "objective": "maximize uniform PLL measurement delay tolerated by exact linear DDE at fixed nodal rho, optimizing only in frozen gain box if materiality gate passes",
    "spectral_margin_s_inv": -0.05,
    "uniform_delay_transition_grid_ms": [20,25,30,32,34,36,38,40],
    "tau_critical_bracket_target_ms": 0.1,
    "fixed_gain_replacement_percent_levels": [87.5,87.75,88.0,88.25,88.45514078456176],
    "fixed_gain_rho_path": "linear interpolation in generator-MW-weighted replacement between stored seed A and stored best-zero-delay B rho vectors; use exact solved nodal vector at each requested percent",
    "fixed_gain_Kp_Ki": "uniform stored seed A gains, held fixed at all five replacement levels",
    "materiality_gate_ms": {"weak_below":2.0,"strong_at_least":5.0},
    "allow_positive_delay_nonlinear_safety_claim": False,
    "do_not_commit_or_push": True,
    "source": source,
}
(HERE / "FROZEN_PROTOCOL.json").write_text(json.dumps(protocol,indent=2)+"\n",encoding="utf-8")
print("FROZEN_PROTOCOL", HERE / "FROZEN_PROTOCOL.json")
