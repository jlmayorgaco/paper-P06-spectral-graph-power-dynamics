# Model and execution provenance

Baseline Git HEAD: dfbfae6f5b36898d276780ace80fcba67bcc0a24.
The working tree was dirty before this work; INPUT_MANIFEST.json records it.
Unrelated edits were neither staged nor changed by this campaign.

Physical model: the canonical fixed-support IEEE39 SG/GFL sharing model. Its
204-state differential model follows exact elimination of the algebraic voltage
variables in a regular chart, not truncation of dynamic states. The234-coordinate
affine descriptor retains204states,20voltage coordinates and10PLL detector
coordinates. Ten active generator sites are buses30–39. Initialized generator
dispatch totals5402.761089978847MW; never substitute net bus injections.

rho semantics: SG retained dispatch and rating scale by1-rho; GFL by rho;
the archived equilibrium sharing/current law is used unchanged. The exported
parametric Jacobian is verified against the canonical nonlinear ReducedDAE
implementation in Julia by ForwardDiff. This is an independent differentiation
path through the same physical model, not a second independent simulator.

Delay channel: only the actual PLL detector error e=-sin(theta)*v_r+cos(theta)*v_i.
The PLL frequency filter tf=1/(300*2pi) remains modeled. The gain channels are
Kp*e_delayed/tf and Ki*e_delayed. No delayed power/current/setpoint substitution.
All executed cases have tau_i=.04s; tau is not a decision variable. Spectral
characteristics and nonlinear method-of-steps histories use the pure delay.

Inputs: experiments/graph_gsp_codesign_20261003/model/*.csv, baseline.toml and
the pre-existing root catalog. Exact source paths and SHA256 are recorded in
INPUT_MANIFEST.json and DEPENDENCY_MANIFEST.json. The archived Julia sources,
Project/Manifest and frozen operating data are included in DEPENDENCIES.zip.

Security contract: all finite physical roots satisfy alpha<=-.05/s after
deflating rotational gauge; numerical full-contour check is distinguished from
certified local root boxes. Nonlinear5case contract uses nominal constant-Z
load changes (8,-100),(16,+100),(16,-100),(29,+100),(29,-100)MW,60s horizon,
Rodas5P method of steps, dtmax.01s,tolerance1e-9. All39bus phase measurements
use.5s windows. Limits:F<=.5Hz,RoCoF<=.5Hz/s,V[.9,1.1],SGslack>=.002.
No converter current limiter is included. Julia versions are in JULIA_VERSIONS.csv.

Rigorous arithmetic: python-flint0.8.0,160bits, exact binary64 source inputs.
Its resolvent series is an enclosure of numerical evaluation, not a delay
approximation. Physical/model/CSV-export uncertainty is not enclosed.
The full-spectrum Julia contour is floating point, not an interval certificate.

Reproduction status: eight endpoint root certificates, four common counts,
six independent differentiation/parity cases, six numerical full spectra,
and five nonlinear events. Source/figure/table relationships are reproducible;
native editor compiler is unavailable, installed MiKTeX compiled successfully.
