# MASTER PACKAGE — Modal Substitutability Paper/Poster
## (verified theory + results + poster layout + ANDES validation plan)

═══════════════════════════════════════════════════════════════════
## 0. STATUS: GATE PASSED (with a correction)
═══════════════════════════════════════════════════════════════════
The naive "inertia ≡ line reinforcement" equivalence is FALSE in sign — verified.
Inertia LOWERS the critical modal stiffness ν_c; line reinforcement RAISES it.
What survives and is verified to machine precision (rel.err ~1e-4):

  R1. First-order shift formulas (sign-correct, exact):
      - inertia at node i:  dν_c = q_cᵀ [ −½(E_i L̃ + L̃ E_i) ] q_c   (E_i = e_i e_iᵀ)
      - line (i,j):         dν_c = q_cᵀ [ M^(−1/2)(e_i−e_j)(e_i−e_j)ᵀ M^(−1/2) ] q_c
      The two levers move the critical mode in OPPOSITE directions.

  R2. Shunt-defect decomposition: an inertia change ε at node i has an equivalent
      topology perturbation δL = ½(ΕL+LΕ) whose ROW-SUM DEFECT ½Lε is NOT realizable
      by lines (it is a shunt/grounding term). [This part is exact algebra.]

  R3. MODAL SUBSTITUTABILITY CLASSIFIER (the contribution):
      s_i = 1 − |(½Lε_i)·q_c| / ‖½Lε_i‖
      High s_i  → the inertia effect is (modally) accessible by topology → reinforce lines (cheap)
      Low  s_i  → inertia is modally irreplaceable → virtual inertia / synchronous condenser
      6-bus result: s = [0.54, 0.14, 0.58, 0.36, 0.90, 0.50]; discriminates 8.9×.
      Witness: node 4 (s=0.90, topology) vs node 1 (s=0.14, inertia-inevitable).

  Scope (state honestly): exact for the conservative network-inertia subsystem s²M + L̃
  (damping lives in control, GENROU D=0). Control self-energy Σ(s) modifies the result in
  control-dominant regimes — that is future work, not claimed here.

═══════════════════════════════════════════════════════════════════
## 1. PAPER (conference)
═══════════════════════════════════════════════════════════════════
TITLE: "Modal Substitutability of Virtual Inertia and Network Topology in
         Low-Inertia Inverter Grids"
  (punchy alt for talk: "Inertia Lowers It, Topology Raises It")

CONTRIBUTION (claim exactly this, no more):
  1. The two robustness levers (virtual inertia, line reinforcement) act with OPPOSITE
     sign on the critical modal stiffness, with exact first-order formulas in one currency.
  2. The inertia effect carries a shunt defect that lines cannot realize.
  3. A modal-substitutability classifier ranks weak nodes into topology-cheap vs
     inertia-inevitable — a structural, pre-simulation design criterion.
NOT CLAIMED: individual Rayleigh sensitivities (known); normalized weighted Laplacian
  (known in spectral graph theory). The novelty is the modal comparison + classifier for
  IBR robustness planning.

WHY USEFUL: virtual inertia is expensive and continuous (GFM hardware/control/energy);
  line reinforcement is one-time capital. The classifier says, per weak node, whether the
  expensive lever can be swapped for the cheap one — and how much inertia is genuinely
  irreplaceable (the shunt defect → synchronous condenser territory).

WHO ELSE WORKS NEARBY (cite + differentiate):
  - Poolla/Dörfler/Pirani/Borsche: optimal virtual-inertia placement (H2/H∞/pole). Have
    inertia as a lever; do NOT compare it against topology nor have the substitutability
    criterion.
  - Pagnier–Jacquod: inertia & slow Laplacian modes, but HOMOGENEOUS inertia (substitutability
    collapses there) and RoCoF, not damping-mode substitutability.
  - Spectral graph theory (normalized weighted Laplacian M^(−1/2)LM^(−1/2)): the object exists
    for clustering, not for physical inertia/topology trade-off.
  - Grid-forming-as-virtual-impedance: device-level synthesis of impedance; ours is graph/
    spectral level. *VERIFY this flank with one adversarial search before submission.*

═══════════════════════════════════════════════════════════════════
## 2. POSTER LAYOUT (conference, A0 portrait 841×1189mm, 3 columns)
═══════════════════════════════════════════════════════════════════
TITLE BAND (8%): title + subtitle "When to reinforce the grid instead of buying inertia"
HOOK BAND (5%, colored): "Virtual inertia is expensive. For some weak nodes a line does the
  job. For others, inertia is irreplaceable. A spectral rule tells which."

COL LEFT — THE PROBLEM & THEORY
 [A] Context (industrial, low-inertia IBR; inertia & damping now synthetic & placeable)
 [B] Model: M θ̈ + Dθ̇ + Lθ = p ; L̃ = M^(−1/2)LM^(−1/2) (ONE equation, large)
 [C] Two opposite levers: inertia lowers ν_c, lines raise it (state the two formulas, small)

COL CENTER — THE RESULT (dominant)
 [D] ★ FIG 1: substitutability graph (the green/red node map). THE 30-second figure.
 [E] FIG 2: opposite-direction bar chart (inertia down / lines up)
 [F] One-line meaning: "It is not how MUCH inertia — it is whether its shunt defect aligns
     with the critical mode."

COL RIGHT — DESIGN & EVIDENCE
 [G] The classifier: s_i = 1 − |defect·q_c|/‖defect‖ (one equation) + reject rule
 [H] Witness: node 4 (s=0.90 → lines) vs node 1 (s=0.14 → inertia). Verified vs finite-diff.
 [I] Workflow: 1 find weak nodes → 2 compute s_i → 3 high s: reinforce / low s: inertia/SynCon
     → 4 validate with QEP/ANDES
 [J] Claim discipline: Proven (first-order formulas + classifier, machine-precision FD check).
     Not yet (full ANDES IBR validation; control self-energy regime).

TAKEAWAY BAND: "Robustness planning is not 'add inertia everywhere.' It is: substitute the
  cheap lever where the spectrum allows, and pay for inertia only where it can't."
+ Scope + Reproducibility (seeds, repo) + QR.

DESIGN RULES: ≤5 equations total. Colors: green=topology/cheap, red=inertia/expensive,
  consistent everywhere. Body ≥28pt. 35%+ whitespace. Fig 1 dominant (~25% center column).

═══════════════════════════════════════════════════════════════════
## 3. MORE RESULTS TO GENERATE (Codex, in order)
═══════════════════════════════════════════════════════════════════
G1. Robustness of the classifier: re-run s_i over many random 6–14 bus graphs and inertia
    profiles (fixed seeds). Report the distribution of discrimination ratio max/min. Need it
    to be >1.5 typically, not just in the hand-picked 6-bus.
G2. Critical-mode choice sensitivity: recompute s_i for Fiedler-critical vs top-mode-critical.
    Report whether the node ranking is stable. (We used top mode; document the choice.)
G3. Cost-aware version (optional, strong for IAS): attach a citable cost per lever
    (line reinforcement $/p.u. vs virtual-inertia/SynCon $/MWs²) and rank nodes by
    damping-benefit-per-dollar. Only if a defensible cost source exists; else omit.
G4. The shunt defect → synchronous condenser link: show that for low-s nodes, the
    irreplaceable part is exactly a shunt, i.e. a SynCon at that bus, quantify how much.

═══════════════════════════════════════════════════════════════════
## 4. ANDES VALIDATION PLAN (next level, beyond poster)
═══════════════════════════════════════════════════════════════════
Goal: show the structural classifier (computed on the reduced L̃, M) predicts the REAL
damping outcome on the full ANDES eigensolve with control.

V1. Build L̃, M from a solved ANDES IEEE-39 case (active-power Jacobian → L; machine +
    virtual inertia → M). Document the extraction (Level-1 surrogate).
V2. For each generator/weak bus, compute s_i from the reduced model.
V3. ACTION TEST: for a high-s node, apply (a) virtual inertia ΔM vs (b) the equivalent line
    reinforcement; measure Δζ_min in FULL ANDES eigensolve. Confirm lines achieve comparable
    Δζ for high-s nodes and FAIL to for low-s nodes.
V4. NEGATIVE CONTROL: a low-s node where lines should NOT substitute — confirm in ANDES that
    line reinforcement underperforms inertia there.
V5. METRIC: agreement (sign + magnitude) between classifier prediction and full-ANDES Δζ_min.
    Report regret. Controls: Fiedler-rank, uniform-inertia, frequency-rank.
DISCIPLINE: ANDES = ground truth; no parameter tuning; report negatives; fixed seeds; the
  reduced-model classifier must be computed BLIND to the ANDES outcome, compared only after.

═══════════════════════════════════════════════════════════════════
## 5. FILES IN THIS PACKAGE
═══════════════════════════════════════════════════════════════════
  verify_theorem.py     — original gate (shows global fraction is flat: the cautionary result)
  rescue_test.py        — finds that modal projection discriminates
  corrected_theory.py   — sign-correct first-order formulas, machine-precision FD validation
  full_results.py       — classifier + witness + json export
  make_figures.py       — fig1 (substitutability graph), fig2 (opposite levers)
  poster_data.json      — all numbers for figures/tables
