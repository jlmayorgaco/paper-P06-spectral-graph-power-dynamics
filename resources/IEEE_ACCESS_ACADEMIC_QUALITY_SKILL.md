---
name: ieee-access-academic-quality
version: 1.0
purpose: >-
  Produce, audit, and revise high-quality IEEE Access engineering research articles
  with the standards of a skeptical technical reviewer, a reproducibility reviewer,
  and an expert reader. Optimized for power systems, PMU/synchrophasor research,
  physics-informed learning, Bayesian/dynamic estimation, and simulation-heavy work,
  but applicable to related engineering papers.
---

# IEEE Access Academic Quality Skill

## 1. Mission

Produce a paper that is not merely correct or publishable, but memorable, technically defensible, reproducible, visually intelligible, and genuinely useful to the research community.

The target is a paper for which a demanding reader can answer, after one reading:

1. **What exact unresolved problem was addressed?**
2. **Why were existing methods insufficient for that problem?**
3. **What is the one central technical idea?**
4. **What evidence shows that the idea works, and where does it fail?**
5. **What can another researcher now do that they could not do before?**

A great paper is not a collection of methods, figures, equations, and benchmarks. It is a **compressed scientific argument**:

> **problem -> gap -> insight -> mechanism -> test -> evidence -> boundary -> consequence**

Every section, equation, figure, table, experiment, and citation must serve that argument.

---

# 2. Evidence Basis for This Skill

This skill is informed by:

- Current IEEE Access author and submission guidance.
- IEEE Access reviewer criteria and reviewer best practices.
- IEEE Author Center guidance on article structure, equations, figures, references, and reproducibility.
- IEEE Access reproducibility criteria: documentation, completeness, and exercisability.
- A close structural review of representative IEEE Access papers in PMU/power-system research, including:
  - *Transfer Learning for Event Detection From PMU Measurements With Scarce Labels* (IEEE Access, 2021).
  - *An Open-Access Repository of Synchrophasor Data Quality Examples: Curation and Example Applications* (IEEE Access, 2025).
  - *Model Validation and Calibration of Large Power Systems Using Hierarchical Topological Cuts* (IEEE Access, 2026).
  - Related IEEE Access event-classification, model-validation, state-estimation, and PMU-data papers from approximately the last six years.

Patterns observed in strong papers include:

- a narrow scientific question rather than an encyclopedic scope;
- a related-work section that leads directly to the contribution;
- explicit comparison to prior art rather than citation catalogues;
- a method described in enough detail to reproduce;
- experiments organized around hypotheses/questions rather than arbitrary experiment numbers;
- multiple regimes or stress conditions rather than one favorable benchmark;
- meaningful baselines and ablations;
- figures that carry scientific reasoning, not decoration;
- quantitative results linked directly to claims;
- explicit limitations and deployment boundaries;
- synthetic evidence supplemented by independent/realistic evidence when feasible;
- conclusions that are narrower than the ambition of the introduction, not broader.

Official IEEE Access reviewers explicitly assess contribution to knowledge, technical soundness, comprehensiveness, sufficiency of references, quality of benchmarking and validation, reproducibility, and whether conclusions are supported by the data. This skill treats those as hard design constraints rather than end-stage checks.

---

# 3. Non-Negotiable Scientific Principles

## 3.1 One paper, one central scientific claim

A paper may have 3–4 contributions, but they must support one scientific thesis.

Bad scope:

> DAE + Bayes + GNN + Transformer + federated learning + continual learning + edge hardware + open-set detection + cyber security.

Good scope:

> Sparse PMU event diagnosis is reformulated as physical-intervention inference; a hybrid physics/learning method is tested for source generalization and diagnosability under changing sensor availability.

If two contributions could become independent papers without weakening either story, strongly consider splitting them.

## 3.2 Novelty is not a list of ingredients

Never claim novelty from combining known terms alone.

Do not write:

> “We combine physics, Bayesian inference, machine learning, and edge computing.”

Instead identify a precise unresolved contract:

> “Existing source-labeled models cannot localize assets absent from training; the proposed formulation evaluates those assets as physical hypotheses and tests whether learned model discrepancy improves the physical likelihood without sacrificing source generalization.”

The contribution is the **new capability, theorem, formulation, or evidence**, not the presence of fashionable components.

## 3.3 Every claim needs an evidence class

Maintain a claim–evidence register. Every important claim must be one of:

- **Theoretical**: follows from stated assumptions and proof.
- **Empirical-primary**: supported by the main frozen test.
- **Empirical-secondary**: supported by ablation/stress/transfer evidence.
- **Historical**: supported by prior work or legacy experiments.
- **Contextual**: supported by literature.
- **Hypothesis/future work**: not yet established.

Never upgrade one class into another in prose.

## 3.4 Negative results are scientific assets

If a method fails under unseen sources, domain shift, sensor loss, noise, long-stream scanning, or calibration, preserve that result.

A trustworthy paper can say:

> “The electrical descriptor improves line-transfer ranking but does not improve shifted ANDES performance.”

This is stronger science than hiding the negative condition and claiming general superiority.

## 3.5 Do not confuse task definitions

Keep distinct:

- event detection;
- event-family classification;
- physical-source localization;
- measurement-integrity diagnosis;
- state estimation;
- onset estimation;
- alarm-stream precision;
- event recall;
- uncertainty calibration;
- operational severity.

A method that succeeds at one does not automatically solve the others.

---

# 4. Before Writing: Required Research Map

Do not start by polishing prose. First construct the research map.

## 4.1 Formulate the paper in five sentences

Write exactly five internal sentences:

1. **Problem** — what real scientific/engineering problem remains unresolved?
2. **Gap** — why do the strongest existing approaches fail to solve this exact version?
3. **Insight** — what new idea changes the problem?
4. **Test** — what experiment could falsify that idea?
5. **Consequence** — if successful, what becomes possible?

If these five sentences are vague, the paper is not ready to write.

## 4.2 Build a literature capability matrix

For each directly relevant work, extract at minimum:

| Field | Required extraction |
|---|---|
| Problem | exact task solved |
| System/data | field, HIL, simulator, benchmark, number of PMUs/sensors |
| Event classes | what physical/integrity phenomena are included |
| Physics | model-based, physics-guided features, constraints, none |
| Learning | model family and what it learns |
| Inference | deterministic score, probability, posterior, residual, optimizer |
| Temporal | causal, offline, retrospective, windowed |
| Source generalization | known sources only or held-out physical sources |
| Sensor configuration | fixed, missing-data robust, dynamic join/leave |
| Uncertainty | calibrated or uncalibrated; abstention or forced decision |
| Transfer | operating-point, grid, simulator, utility, topology |
| Baselines | what was actually compared |
| Metrics | exact denominators and units |
| Reproducibility | code/data/seeds/manifests |
| Limitation | strongest limitation stated or inferred |
| Relevance to this paper | baseline, component, contrast, or excluded alternative |

Do not cite a paper as a baseline unless its experiment can be reproduced under a comparable contract.

## 4.3 Literature search strategy

For a modern engineering paper:

- Search the most recent ~6 years aggressively.
- Retain older seminal papers where they define foundational methods, standards, or theory.
- Search IEEE Xplore and Web of Science/Scopus where available.
- Use Google Scholar for forward/backward citation tracing.
- Use arXiv/preprints only for genuinely emerging results and identify them as such.
- Search direct competitors by method and by problem, not only by keywords from the proposed title.
- Search for negative/limitation papers, not only successful methods.
- Search adjacent fields for mature concepts when useful (sequential testing, Bayesian evidence, information geometry, continual learning, distributed inference), but connect them rigorously to the power-system problem.

### Saturation rule

The review is not complete until new search queries stop revealing **new methodological families or materially different evaluation contracts**.

### Citation verification rule

Before a reference enters the paper:

1. verify that the work exists;
2. verify title, authors, venue, year, DOI when available;
3. inspect enough of the full work to confirm the claimed contribution;
4. ensure the cited sentence is actually supported;
5. check that the work is not retracted.

Never cite from title/abstract alone when the claim depends on methodology or limitations.

---

# 5. What Strong IEEE Access Papers Commonly Do

The following are **house rules derived from observed high-quality papers**, not IEEE mandates.

## 5.1 Page allocation for an 18–20 page Research Article

A strong technical paper in this domain typically benefits from approximately:

| Section | Internal target |
|---|---:|
| Abstract/front matter | 0.4–0.6 page |
| Introduction | 1.2–1.6 pages |
| Related Work / Gap | 1.3–2.0 pages |
| Problem formulation / theory | 2–3 pages |
| Proposed method | 2–3 pages |
| Experimental methodology | 1.3–2 pages |
| Results | **4–5 pages** |
| Discussion / limitations | 0.8–1.2 pages |
| Conclusion | 0.3–0.5 page |
| References + biographies | 2–3 pages |

The Results section should normally be one of the largest sections.

A paper with 5 pages of conceptual method and 1 page of evidence is usually under-validated.

## 5.2 Reference density

IEEE does not impose a target reference count. Use only relevant references.

For a broad 18–20 page PMU / physics-informed / Bayesian / ML paper, a useful internal target is often **~45–70 verified references**, depending on scope.

The bibliography should include:

- seminal foundations;
- direct modern competitors;
- recent state of the art;
- methods used as baselines;
- datasets/simulators/standards;
- uncertainty/statistical methodology when central;
- relevant negative or robustness work.

Do not pad references. A shorter accurate bibliography is better than a long generic one.

## 5.3 Contribution count

Prefer **2–4 contributions**.

Each contribution should satisfy this grammar:

> **We [technical action] that [new capability], and we demonstrate/test it by [specific evidence].**

Weak contribution:

> “A comprehensive literature review is presented.”

Strong contribution:

> “We formulate PMU source diagnosis as source-agnostic physical-intervention inference and evaluate it on assets excluded from all source-labeled training.”

---

# 6. Section-by-Section Quality Standard

## 6.1 Title

IEEE guidance: specific, concise, descriptive; avoid “new” and “novel.”

A strong title contains:

- the core method or concept;
- the core task;
- the application domain.

Avoid listing every technology.

Bad:

> “A Novel Comprehensive Robust AI-Based Framework for Smart Grid Monitoring.”

Better:

> “Physics-Informed Bayesian Event Inference and Localization From Sparse PMU Measurements.”

### Title test

A researcher should know from the title whether the paper is relevant without reading the abstract.

---

## 6.2 Abstract

Official IEEE guidance: one paragraph, self-contained, up to 250 words; no references, footnotes, or numbered equations.

Internal target: **~180–220 words**.

Recommended structure:

1. 1–2 sentences: problem and why it matters.
2. 1–2 sentences: specific unresolved gap.
3. 2–3 sentences: proposed method/mechanism.
4. 3–4 sentences: the most important quantitative results.
5. 1 sentence: implication and boundary.

The final abstract must contain actual evidence. Do not submit a final abstract that says “planned evaluation,” “will be tested,” or “expected to.”

### Abstract memory test

After reading only the abstract, an expert should remember:

- one problem;
- one idea;
- two or three important numbers;
- one implication.

---

## 6.3 Introduction

The introduction is an argument, not background lecture notes.

Recommended seven-paragraph logic:

1. **Operational/scientific problem.**
2. **What strong existing methods already solve.**
3. **Precise remaining gap.**
4. **Why the gap matters and why obvious solutions fail.**
5. **Core technical insight.**
6. **2–4 contributions.**
7. **Headline evidence / scope boundary.**

Avoid spending half a page explaining what a PMU is to an IEEE power-systems audience.

Avoid generic openings such as:

> “Power systems are the backbone of modern society.”

### Introduction quality test

At the end of the introduction, the skeptical reader should be able to state:

> “I understand exactly what these authors claim is missing, and I know what experiment would prove them wrong.”

---

## 6.4 Related Work / Literature Review

Do not write a citation catalogue.

Bad pattern:

> A did X. B did Y. C did Z.

Preferred pattern:

> This family solves X under assumptions A and B. It is strong when C is available, but published evidence does not establish D. This motivates the specific contract tested here.

Organize by **methodological families and unresolved capabilities**, e.g.:

- PMU event detection/classification/localization;
- model-based and dynamic estimation;
- physics-informed / graph learning;
- sparse observability / placement / identifiability;
- Bayesian / hybrid event-state inference;
- open-set, domain shift, uncertainty, and data quality;
- edge/federated/continual inference.

Finish with a compact gap statement or capability matrix.

### Related-work rule

For every major family answer:

1. What does it solve?
2. What assumptions does it require?
3. What kind of generalization has actually been demonstrated?
4. What failure mode remains relevant to this paper?
5. Is it a baseline, component, or excluded alternative here?

### Fairness rule

Do not invent weaknesses in prior work to make the proposed method look better.

If a competitor solves a different problem, say so rather than claiming superiority.

---

## 6.5 Problem Formulation

A strong formulation fixes the scientific contract before the algorithm begins.

Define explicitly:

- system model;
- observed variables;
- hidden variables;
- event/intervention variables;
- measurement/integrity variables;
- hypothesis/candidate space;
- causal data boundary;
- outputs;
- what counts as correct;
- what is outside scope.

Readers should not discover in Results that “localization” actually means a different denominator or candidate set.

---

## 6.6 Theory and Equations

Use equations because they remove ambiguity, not because mathematical density looks sophisticated.

### Equation test

Every numbered equation must do at least one of:

- define the plant/measurement model;
- define a new quantity;
- state an optimization/inference objective;
- define a decision rule;
- support a theorem/proposition;
- define an evaluation measure unique enough to matter.

Delete decorative equations that merely restate textbook material.

### Theorem/proposition quality pattern

For every formal result provide:

1. assumptions;
2. precise statement;
3. proof or justified derivation;
4. interpretation;
5. limitation;
6. experiment that tests whether the approximation is useful in the nonlinear/noisy setting.

Prefer **one strong proposition and one meaningful corollary** over many trivial lemmas.

Do not claim novelty for standard KL, Fisher, Mahalanobis, Bayes, or information-additivity facts.

### Theory-to-evidence bridge

A theoretical quantity must predict or explain something empirical.

Example:

> pairwise information -> posterior ambiguity -> candidate-set size -> time-to-resolution.

If the theory has no measurable consequence, question why it belongs in the main paper.

---

## 6.7 Proposed Method

The method section must be executable in the reader’s head.

Specify:

- inputs;
- state retained between time steps;
- initialization;
- preprocessing;
- hypothesis/candidate generation;
- learned components;
- physical components;
- objective/likelihood/loss;
- thresholds and where they are selected;
- inference update;
- stopping/abstention rule;
- output;
- computational complexity;
- causal/noncausal boundary.

If the algorithm is complex, include one concise algorithm block or flow diagram.

### Physics-informed ML rule

Clearly state **where physics enters** and **what ML learns**.

Distinguish:

- physics-aware feature engineering;
- physics-guided constraints;
- physics-informed architecture;
- learned model discrepancy;
- physics-based likelihood plus learned correction.

Do not use “physics-informed” as a decorative adjective.

### Bayesian rule

Do not call scores “probabilities” or “posteriors” unless a probabilistic model and calibration/normalization justify it.

---

## 6.8 Experimental Methodology

The experimental protocol must be designed to falsify claims, not maximize accuracy.

### Required structure

For each campaign define:

- research question;
- independent experimental unit;
- data source/simulator;
- train/development/test split;
- leakage prevention;
- source/event exclusion rule;
- operating conditions;
- noise/integrity conditions;
- seeds/repetitions;
- baselines;
- metrics;
- uncertainty/statistical analysis;
- threshold selection;
- frozen test protocol.

### Hard anti-leakage rule

Time-series windows or samples from the same physical trajectory/event parent must not leak across train/test merely because rows were randomly shuffled.

For source-generalization claims, the **physical source asset itself** must be excluded from source-labeled training/tuning, not just a replica.

### Development/test rule

Tune thresholds and model choices on development data.

Freeze them.

Evaluate the final test once.

Do not select the best test row retrospectively.

### Independent-unit rule

Confidence intervals and significance tests must use the independent unit (trajectory, event, parent record, operating condition), not millions of correlated PMU samples.

### Multiple regimes

A strong simulation-heavy paper should usually include more than one of:

- nominal condition;
- operating-point shift;
- parameter shift;
- unseen physical source;
- sensor loss;
- noise / corruption;
- cross-simulator transfer;
- long normal exposure;
- independent field/HIL evidence when available.

---

# 7. Baselines, Ablations, and Comparisons

## 7.1 Baselines must answer scientific alternatives

Include baselines representing different paradigms when relevant:

- legacy/feature-based ML;
- strong contemporary data-driven model;
- physics-only/model-based method;
- hybrid physics/learning method;
- dynamic estimator;
- proposed method.

Do not include a weak baseline just to win.

## 7.2 Explain why each baseline is included or excluded

For every major literature method, record:

- directly comparable? yes/no;
- implementation available? yes/no;
- data contract compatible? yes/no;
- included baseline? why/why not.

This prevents unfair claims such as comparing a bus-localization method against a paper that only classifies event family.

## 7.3 Ablation standard

Ablations should map to contribution claims.

Example hybrid stack:

- ML only;
- physics only;
- physics + legacy engineered features;
- physics + learned discrepancy;
- + Bayesian uncertainty;
- + diagnosability-aware abstention.

Each ablation asks **what scientific ingredient contributes**, not merely whether deleting a random layer changes accuracy.

---

# 8. Metrics and Statistics

## 8.1 Declare denominator and unit

Every metric must identify its population.

Bad:

> “Localization accuracy = 91%.”

Good:

> “Exact physical-source Top-1 = 91% over 131 held-out event parents.”

## 8.2 Do not let one metric hide another task

Report separately when relevant:

- detection precision/recall/F1;
- macro-F1 event family;
- Top-1/Top-3 exact source;
- joint event + source correctness;
- candidate-set coverage/size;
- Brier/NLL/ECE for calibrated probabilistic outputs;
- false alarms per hour;
- event recall;
- joint alarm precision;
- detection delay;
- source-resolution delay;
- runtime/memory/communication.

## 8.3 Confidence intervals

Use confidence intervals for important differences, ideally at the scenario/parent level.

Prefer paired comparisons when methods are evaluated on the same cases.

Do not use p-values as decoration. Report effect size and uncertainty.

## 8.4 Accuracy is not operational reliability

High event recovery can coexist with alarm flooding.

Any continuous-monitoring paper must report alarm burden on normal exposure and clarify how duplicate/repeated detections are scored.

---

# 9. Results Section

Results should be organized around **scientific questions**, not experiment file names.

Bad:

- Experiment 1
- Experiment 2
- Experiment 3

Better:

- Can unseen physical sources be localized?
- Does the theoretical diagnosability quantity predict empirical ambiguity?
- What changes when informative PMUs are lost or added?
- Does learned discrepancy improve transfer without damaging source generalization?
- How does the method behave in continuous scanning?

### Result paragraph pattern

Each result paragraph should normally contain:

1. question;
2. observation/number;
3. comparison;
4. interpretation;
5. limitation or exception when relevant.

Example:

> “Under complete source holdout, the hybrid model retains X% Top-1 compared with Y% for the source-labeled baseline. The improvement is concentrated in fault and outage hypotheses, while load localization remains ambiguous. This supports source-generalizable physical ranking but does not establish transfer to unseen network topology.”

### No future tense in final Results

The submission-ready paper must not contain:

> “Campaign A will report…”

or

> “Results will be evaluated…”

---

# 10. Discussion and Limitations

The discussion is not a second results section and not a marketing section.

It must answer:

1. What did the evidence establish?
2. Why did it happen?
3. Which competing explanation remains possible?
4. Where does the method fail?
5. Which assumptions matter most?
6. What would have to be true for operational deployment?
7. What should the next experiment test?

A strong limitation increases credibility when it is precise.

Bad:

> “Future work will improve robustness.”

Good:

> “The cross-simulator campaign covers line interventions under two operating blocks; it does not establish multi-family transfer or robustness to independently calibrated PMU current semantics.”

---

# 11. Conclusion

The conclusion must answer the paper’s question using only evidence already shown.

Recommended structure:

1. one sentence: problem;
2. one sentence: central mechanism;
3. two or three sentences: strongest results;
4. one sentence: main limitation;
5. one sentence: implication/next step.

No new references, experiments, theories, or claims.

Do not repeat the abstract verbatim.

---

# 12. Figures and Tables

## 12.1 Every figure must earn space

A figure must do one of:

- explain the problem;
- explain the method;
- reveal a pattern;
- compare alternatives;
- validate a theory;
- show a failure mode;
- quantify a tradeoff.

Delete decorative diagrams.

## 12.2 Conceptual vs empirical figures

Use TikZ/vector diagrams for:

- system architecture;
- inference flow;
- generative models;
- theoretical mechanisms;
- sensor join/leave logic.

Generate empirical plots programmatically from frozen results.

Never hand-draw numerical results.

## 12.3 IEEE graphic quality

Prefer vector PDF/EPS/PS for line art when possible.

For raster:

- color/grayscale >300 dpi;
- black/white line art >600 dpi.

Design at approximately:

- 3.5 in one-column width;
- 7.16 in two-column width.

Figures must remain interpretable in grayscale and for readers with color-vision deficiency. Use marker shape, line style, and brightness in addition to color.

## 12.4 Visual hierarchy

Use consistent:

- fonts;
- line widths;
- axis styles;
- label capitalization;
- abbreviations;
- symbol meanings;
- panel labeling.

Avoid tiny text, giant legends, rainbow palettes, 3D bars, shadows, gradients, and PowerPoint-style boxes.

## 12.5 Figure density

For an 18–20 page technical paper, approximately **7–10 high-information numbered figures** is often a productive range, frequently using multi-panel composites. This is a house guideline, not an IEEE rule.

Prefer one excellent composite over four trivial plots.

## 12.6 Tables

Use tables for exact values and experimental contracts.

A strong paper often needs only ~4–7 main tables:

- literature/gap matrix;
- model/event definitions;
- dataset/experimental contract;
- main benchmark;
- ablation;
- robustness;
- runtime/resources.

Push exhaustive hyperparameters and giant tables to supplementary material.

---

# 13. Writing Style: Technical, Human, Restrained

## 13.1 Default voice

Write clear technical English with restrained confidence.

Prefer:

> “Electrical descriptors improve fault localization under the tested transfer condition.”

over:

> “The proposed innovative framework demonstrates remarkable robustness and superior performance.”

## 13.2 Prohibited AI-style habits

Avoid or aggressively reduce:

- “In today’s rapidly evolving…”
- “plays a pivotal/crucial role” when the statement is generic;
- “leverages the power of…”
- “seamlessly”;
- “groundbreaking”;
- “revolutionary”;
- “comprehensive framework” as self-praise;
- “remarkable” without quantitative justification;
- repetitive “Furthermore / Moreover / Additionally” paragraph starts;
- repeated three-item rhetorical lists;
- excessive em dashes;
- a summary sentence after every paragraph;
- restating the section title as the first sentence;
- “It can be clearly seen…” instead of describing the result;
- fake quotations or imagined reviewer language;
- inflated novelty wording.

## 13.3 Paragraph discipline

A paragraph should have one technical function.

Typical paragraph:

- topic claim;
- supporting mechanism/evidence;
- consequence or transition.

Avoid paragraphs that merely collect unrelated facts.

## 13.4 Sentence discipline

Prefer subject–verb–object and concrete nouns.

Keep long sentences only when logical dependencies require them.

Do not create artificial sentence fragments for dramatic effect.

## 13.5 Terminology discipline

One concept, one term.

Do not switch casually among:

- anomaly/event/disturbance/contingency;
- source/location/origin;
- confidence/probability/score;
- sensor/PMU/node;
- online/real-time/causal.

Define distinctions and preserve them.

---

# 14. Claim Discipline Dictionary

Use the following words only when the corresponding evidence exists.

## “Physics-informed”
State exactly how physics constrains representation, objective, architecture, likelihood, or learning.

Physics-aware engineered features alone are better described as **physics-guided features** unless the learning process itself is constrained/informed by physics.

## “Bayesian”
Requires a stated probabilistic model/prior/likelihood/posterior or a justified approximation.

Classifier votes are not automatically Bayesian probabilities.

## “Causal”
A decision at time k may use only data available through k.

Retrospective smoothing is not causal.

## “Real-time”
Requires data availability and measured end-to-end runtime/latency consistent with the application.

Fast model scoring alone is insufficient.

## “Distributed”
Computation or inference must genuinely be partitioned across nodes or modeled as such, with fusion assumptions and communication defined.

## “Robust”
Must specify robust to what: noise, topology, parameter shift, missing sensors, adversarial corruption, etc.

## “Generalizes”
State the axis of generalization: new replica, operating point, source, simulator, topology, utility, hardware, event family.

## “Unseen source / zero-shot source”
The physical source asset must be absent from source-labeled training and tuning.

## “State of the art”
Use only after a fair common-contract benchmark against relevant contemporary methods.

## “First”
Use only after an explicit literature search targeted to every clause of the novelty statement; prefer “to the best of our knowledge” when appropriate.

---

# 15. Reproducibility as a Design Constraint

IEEE Access explicitly values reproducible research and evaluates code artifacts on documentation, completeness, and exercisability.

Design the repository so that every paper result has provenance:

> **manifest -> raw/frozen outputs -> aggregation -> table/figure -> manuscript**

## Minimum reproducibility package

- exact environment/dependency specification;
- simulator/model versions;
- network model hash/version;
- seed policy;
- train/dev/test manifests;
- data-generation configuration;
- baseline configurations;
- metric implementation;
- scripts to regenerate tables/figures;
- expected outputs;
- runtime/hardware notes;
- README instructions.

## No manual transcription rule

If a numerical value can be generated, do not type it manually into LaTeX.

Generate macros/tables programmatically from frozen result artifacts.

## Failed-solve accounting

Report both successful-solve denominators and planned-experiment denominators where simulation failures can bias the evaluation.

---

# 16. Three Internal Review Modes

Before considering a manuscript submission-ready, run all three modes independently.

## Mode A — Reviewer 1: Skeptical technical expert

Assume the reviewer is a power-systems/control/inference expert who believes the novelty is overstated until proven otherwise.

Ask:

1. Is the central mathematical/physical formulation correct?
2. Are all assumptions explicit?
3. Is any theorem trivial, circular, or mislabeled as novel?
4. Does the proposed method actually solve the stated problem?
5. Are causal/Bayesian/distributed claims technically justified?
6. Could results be explained by leakage or source memorization?
7. Are baselines strong and fair?
8. Does an ablation isolate the contribution?
9. Are statistical units independent?
10. Is there a simpler established method that would solve the same problem?
11. Are negative cases investigated?
12. Does the theory predict an empirical effect?

### Reviewer 1 veto conditions

Any of the following blocks submission until fixed:

- central theorem/model error;
- train/test leakage;
- source-holdout claim with source identity leakage;
- future information used in a causal claim;
- classifier score mislabeled as posterior;
- missing direct baseline when implementable;
- primary conclusion unsupported by a frozen test;
- statistically invalid independence assumptions hidden in evaluation.

---

## Mode B — Reviewer 2: Experimental/reproducibility expert

Assume the reviewer accepts the idea but distrusts the evidence.

Ask:

1. Can the experiment be reproduced from the paper/repository?
2. Are data-generation and simulator semantics documented?
3. Are all splits and thresholds frozen before test?
4. Are failure cases counted?
5. Are denominators and units explicit?
6. Are reported differences larger than variation across seeds/scenarios?
7. Does robustness cover plausible field imperfections rather than only Gaussian noise?
8. Is domain shift genuine or merely a new random seed?
9. Are figures traceable to artifacts?
10. Are runtime and hardware measurements meaningful?
11. Are field/known/retrospective/synthetic data labels honest?
12. Are code/data availability claims accurate?

### Reviewer 2 veto conditions

- irreproducible primary numbers;
- missing split manifests;
- cherry-picked examples presented as general evidence;
- no normal exposure for alarm claims;
- no uncertainty across relevant experimental units;
- mixing incompatible campaigns into one aggregate;
- undocumented preprocessing differences between methods.

---

## Mode C — The expert reader looking for “the paper of the decade”

This reader is not impressed by complexity. They want an idea worth remembering.

### 30-second test

After title + abstract + Fig. 1:

- Is the problem unmistakable?
- Is the new idea visible?
- Is there a concrete result worth continuing for?

### 3-minute test

After introduction + architecture + main-results figure:

- Can the reader explain the contribution to a colleague?
- Is one figure memorable enough to teach from?
- Is the main result stronger than “our model is 1% more accurate”?

### Next-day recall test

The paper should leave behind:

- **one sentence**;
- **one figure**;
- **one quantitative result**;
- **one limitation**.

### Citation test

Write the sentence a future paper would use when citing this work.

If that sentence is generic (“Mayorga et al. used AI for PMUs”), the contribution is too weak.

A better citation target is:

> “Mayorga et al. showed that PMU source diagnosability predicts posterior ambiguity under changing sensor coverage and demonstrated source localization on assets excluded from source-labeled training.”

### Engineering-decision test

Ask:

> “Would this result change how a researcher/operator designs, evaluates, or trusts a PMU diagnosis system?”

If no, deepen the contribution rather than decorating the paper.

---

# 17. Quantitative Internal Scoring Rubric

Score the manuscript out of 100 before submission.

| Dimension | Points |
|---|---:|
| Scientific question and significance | 10 |
| Novelty and literature positioning | 12 |
| Mathematical/physical correctness | 12 |
| Method clarity and reproducibility | 10 |
| Experimental design and leakage control | 14 |
| Baseline fairness and ablation quality | 10 |
| Statistical evidence and uncertainty | 8 |
| Results depth and failure analysis | 8 |
| Writing/narrative clarity | 6 |
| Figures/tables and visual reasoning | 5 |
| Reproducibility artifact quality | 3 |
| Claim integrity and limitations | 2 |
| **Total** | **100** |

Interpretation:

- **95–100**: exceptional; still run external human review.
- **90–94**: strong submission candidate.
- **85–89**: promising but substantive revision recommended.
- **80–84**: reviewer rejection risk remains high.
- **<80**: do not submit.

A Reviewer-1 or Reviewer-2 veto overrides the numerical score.

---

# 18. Fatal Rejection Patterns

Treat these as pre-submission blockers:

- the novelty is only a combination of buzzwords;
- related work ignores obvious direct competitors;
- results use the same data for selection and evaluation;
- source generalization is claimed from new replicas of known sources;
- the paper calls retrospective inference “real-time”;
- the paper calls tree/softmax scores “Bayesian probabilities” without justification;
- only accuracy is reported for a continuous alarm problem;
- a simulation-only paper implies field validation;
- one successful example carries the central claim;
- massive feature engineering obscures what causes performance;
- method details are insufficient to reproduce;
- equations are mathematically ornamental;
- figures have tiny text or unreadable legends;
- conclusion is stronger than Results;
- limitations are generic or absent;
- bibliography is thin for the breadth of claims;
- citations are unverified or unrelated;
- final paper still contains template residue, TODOs, “planned evaluation,” or placeholder results.

---

# 19. Quality Workflow for Creating or Revising a Paper

When this skill is invoked, perform the following sequence.

## Stage 1 — Audit

Inspect:

- manuscript;
- repository structure;
- bibliography;
- figures/tables;
- experiment manifests/results;
- relevant prior manuscripts;
- target journal rules.

Produce:

1. one-paragraph scientific diagnosis;
2. novelty/gap map;
3. claim–evidence matrix;
4. missing evidence list;
5. section/page budget;
6. top Reviewer-1 risks;
7. top Reviewer-2 risks.

Do not edit prose yet if the scientific contract is unstable.

## Stage 2 — Literature build

Construct the capability matrix.

Identify:

- 5–10 direct competitors;
- 5–10 foundational/adjacent methods;
- recent robust/data-quality/transfer works;
- direct baselines that must be implemented;
- methods that cannot be fairly compared and why.

Update the gap claim only after this stage.

## Stage 3 — Paper architecture

Freeze:

- working title;
- central question;
- 2–4 contributions;
- section outline;
- figure roadmap;
- table roadmap;
- experimental campaigns;
- main result expected *type* (never fabricate expected value).

## Stage 4 — Theory/method

Check every equation, symbol, assumption, theorem, and algorithm.

Build a notation table internally if the notation is dense.

Test simplified/synthetic cases before running large campaigns.

## Stage 5 — Frozen experiments

Preregister/freeze where practical:

- manifests;
- seeds;
- selection metric;
- thresholds;
- evaluation code;
- baseline configs.

Then run test campaigns.

## Stage 6 — Evidence-first Results

Generate all tables/plots programmatically.

Write Results from the artifacts, not from memory.

Write Discussion only after Results stabilize.

Write abstract last.

## Stage 7 — Adversarial review

Run Reviewer 1, Reviewer 2, and Reader modes.

For each major concern:

- concern;
- why it matters;
- evidence needed;
- concrete fix;
- status.

Do not “answer” a reviewer concern with prose if it actually requires a new experiment.

## Stage 8 — Final integrity pass

Verify:

- title/abstract match actual work;
- no future-tense results;
- claims match evidence class;
- references verified;
- figures readable at final size;
- all denominators explicit;
- no stale template metadata;
- biographies/affiliations correct;
- supplementary material referenced correctly;
- code/data statements accurate;
- AI assistance disclosure handled according to current IEEE policy when applicable.

---

# 20. P07-Specific Application Rules

For the P07 PMU paper, apply these additional rules unless the authors explicitly change scope.

## 20.1 Preserve legacy SGSMA work as baseline and scientific evidence

The hackathon/legacy method should be represented compactly as a **physics-guided hierarchical ML baseline**, grouped into feature families such as:

- local robust/statistical features;
- multiscale temporal/dynamic features;
- physical V/I/frequency/ROCOF and power-proxy features;
- spatial/network/electrical-distance features;
- integrity/missing/bad-data features.

Do not devote journal pages to listing tens of thousands of feature names.

Use the legacy system to answer:

> How far can supervised physics-guided feature engineering go under known-source training?

## 20.2 Reuse legacy features intelligently

A stronger hybrid use is to compute engineered features on **physics-model residuals**:

> residual = measured PMU response - candidate physical prediction.

This allows legacy expertise to become a model-discrepancy representation rather than a source-label memorizer.

## 20.3 Required method paradigms for comparison

Whenever feasible compare:

1. data-driven/legacy ML;
2. physics-only hypothesis inference;
3. hybrid physics + learned discrepancy;
4. full Bayesian/diagnosability-aware proposal.

## 20.4 Required event composition

Keep physical interventions distinct from measurement integrity.

Physical families may include:

- fault;
- line outage;
- generation change;
- load change.

Measurement states may include:

- nominal;
- missing/dropout;
- corrupted/bad data;
- timing/quality issue when modeled.

Do not represent “missing + physical” as one fundamental physical class when the framework can represent both states simultaneously.

## 20.5 Required hard experiments

The strongest version of P07 should prioritize:

- **held-out physical sources**;
- diagnosability vs empirical ambiguity;
- informative PMU loss/join;
- operating/parameter shift;
- cross-simulator transfer;
- long continuous normal exposure;
- explicit false alarms/hour and resolution delay;
- physics-only vs ML-only vs hybrid ablation;
- knowledge-transfer campaign only if sufficiently implemented.

## 20.6 Current manuscript completion gate

The current blueprint is not submission-ready while:

- proposed-method Results remain unexecuted;
- figures/tables are absent or placeholders;
- bibliography is too thin for the scope;
- contributions are still phrased as framework commitments rather than empirically supported findings.

The final manuscript should transform the current conceptual structure into an evidence-driven article rather than merely expanding its prose.

---

# 21. Final Standard

Before submission, be able to state all of the following without qualification:

- **The question is important.**
- **The exact gap is supported by a current literature review.**
- **The contribution is not merely an assembly of known tools.**
- **The method is mathematically and physically coherent.**
- **The experiments could falsify the claims.**
- **The strongest competitors are treated fairly.**
- **The main results survive at least one meaningful shift/stress condition.**
- **The failure modes are visible, quantified, and discussed.**
- **The figures make the argument easier to understand.**
- **The paper is reproducible enough for an independent researcher to rebuild its key results.**
- **Every strong adjective can be replaced by a number, mechanism, or citation.**
- **A reader remembers the paper for an idea, not for a model acronym.**

The goal is not to make the manuscript look impressive.

The goal is to make the scientific contribution difficult to dismiss.
