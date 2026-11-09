# ECM parameter identification with an agentic workflow

A small practice project. The goal is to identify the parameters of an
equivalent circuit model (ECM) of a Li-ion cell from a public drive-cycle
dataset, using least squares (LS), recursive least squares (RLS) and nonlinear
least squares (NLS), and to check that the results are physically sound and
reproducible. An AI coding agent did the coding. I guided it and reviewed the
results.

## 1. Task

- **Data:** `data/ECM1.mat` (public dataset). Source: _TODO: add link/reference_.
  A 2000 s drive cycle at 10 Hz: current −251 to +132 A, voltage 3.06 to
  4.07 V (single cell), plus a 41-point OCV(SOC) table.
- **Model:** Thevenin ECM with 0, 1 or 2 RC branches:

  ```
  V = OCV(SOC) + R0·I + Σ V_i,     dV_i/dt = −V_i/(R_i·C_i) + I/C_i
  SOC = SOC0 + ∫ I dt / (3600·Q)          (I < 0 = discharge)
  ```

- **Unknowns:** R0, R_i, τ_i = R_i·C_i, and also the capacity Q and the initial
  SOC0 (not given with the data).
- **Fit/validation split:** fit on t < 1400 s, validate on t ≥ 1400 s.

## 2. Methods

| Method | Type | Idea | Strength | Weakness |
|---|---|---|---|---|
| **NLS** (`nls.py`) | batch, output error | simulate the model, minimise the voltage error with `scipy.optimize.least_squares` | no noise bias; also fits Q and SOC0; gives standard errors | needs start values and bounds; slower |
| **LS-ARX** (`ls.py`) | batch, equation error | discrete 1RC as `y[k] = a·y[k−1] + b0·I[k] + b1·I[k−1]` | linear, very fast, classic | **biased by voltage noise** (τ much too small) |
| **LS-separable** (`ls.py`) | batch | for a fixed τ the model is linear in R0, R_i → plain LS; scan τ on a grid and refine | unbiased; the τ-scan shows if τ is identifiable | needs Q and SOC0 from elsewhere |
| **RLS-ARX** (`rls.py`) | online | recursive form of LS-ARX with forgetting factor λ | runs online, tracks changes | same noise bias as LS-ARX |
| **RLS fixed-τ** (`rls.py`) | online | τ fixed (from separable LS), RLS tracks R0, R1 and an OCV offset | unbiased, tracks changes | τ must be known |
| **RLS R0 + offset** (`rls.py`) | online | tracks R0 and the offset V − OCV − R0·I | shows drift between model and data | no RC dynamics |
| **PSO** (`pso.py`) | batch, global search | particle swarm on the voltage MSE (log-scale R and τ), optional LS polish | no start values needed; finds global minimum | slow; no standard errors |
| **Difference LS** (`ls.py`) | batch | ΔV = R0·ΔI, OCV cancels out | R0 without Q, SOC0 or OCV | ignores RC dynamics |
| **LS/RLS-ARX2** (`ls.py`, `rls.py`) | batch / online | 2nd-order ARX for 2RC, exact pole mapping | exact without noise | breaks down with noise |
| **SOC-segmented LS** (`segmented.py`) | batch | R0 (and R1) per SOC bin → look-up table | parameters vs SOC | needs enough data per bin |

The last four methods come from the MATLAB reference scripts (section 6).

**OCV model:** two options are compared automatically, linear interpolation
of the table and a 10th-order polynomial of the table (as in the MATLAB
scripts). The pipeline picks the OCV model and the number of RC branches
together, by BIC.

LS and RLS work on `y = V − OCV(SOC)`, so they need Q and SOC0. These are
taken from the recommended NLS model ("SOC set-up").

## 3. Agent workflow

1. **Inspect first, no code.** The agent loaded the .mat file and reported
   signals, units, sampling, sign convention (negative = discharge), the noise
   level (~10 mV) and a first R0 estimate from current steps (2.5 mΩ).
2. **Rules file.** The facts and rules went into `CLAUDE.md`: units, sign
   convention, model equations, plausible parameter ranges, fit/validation
   split, and how the agent must work (plan first, small steps, tests).
3. **Explore before building.** Quick experiments showed two things that
   shaped the design: (a) classic ARX-LS is strongly biased by the voltage
   noise, and (b) on this dataset the RC branch hardly changes the fit.
4. **Build in small modules** (`data`, `ocv`, `model`, `ls`, `rls`, `nls`,
   `checks`, `synthetic`, `plots`, `pipeline`), each with tests.
5. **Synthetic check first.** Every method is tested on synthetic data with
   *known* parameters (same current profile and OCV table, 10 mV noise,
   fixed seed). Only methods that recover the truth there are trusted on the
   real data.
6. **Review.** I reviewed the plan, the code changes and the plots. See
   `notes/agent_log.md` for what was changed and the open questions.
7. **Reference code review.** The agent read the MATLAB course scripts, ported
   the useful methods to Python, and listed their bugs with line numbers
   (`notes/matlab_review.md`). This review found the OCV-model issue that
   explained the validation error (section 5b).

A reusable skill for the repeated task "run a fit and check that it is
physically reasonable" is in `.claude/skills/ecm-fit-check/SKILL.md`.

## 4. Validation

| Check | How | Where |
|---|---|---|
| Data file unchanged | SHA-256 checksum | `tests/test_data.py` |
| Units | voltage range of a single cell in V; SOC as fraction 0–1 | `check_data`, tests |
| Sign convention | net charge < 0, voltage falls, dV/dI at steps > 0; a flipped sign is detected | `check_data`, tests |
| Model equations | zero current → OCV; Coulomb counting units; exact RC step response; discharge lowers V and SOC | `tests/test_model.py` |
| Methods work | parameter recovery on synthetic data, with and without noise | `tests/test_identify.py` |
| Known weakness | ARX-LS bias is written down as a test | `tests/test_identify.py` |
| Plausibility | all parameters > 0, inside ranges, not on a bound | `checks.plausibility` |
| Identifiability | RC branch voltage vs noise; flat τ-scan; NLS standard errors and correlations | `checks.rc_identifiability`, `tau_scan.png` |
| Residuals | RMSE vs noise level, histogram, autocorrelation, slow drift | `checks.fit_quality`, `residuals.png` |
| Generalisation | separate validation window (t ≥ 1400 s) | all methods |
| Model choice | BIC over (OCV model × RC branches), only among models with plausible parameters | `params.json` |
| SOC sign | SOC must fall on discharge (catches a flipped Coulomb-counting sign) | `check_soc_direction` |
| Optimizer | PSO: several seeds, convergence plot, result not on a bound | `tests/test_reference.py` |
| Reproducibility | pipeline runs twice → identical `params.json`; fixed seed; pinned versions | `tests/test_pipeline.py` |

## 5. Results

### 5a. Synthetic data (known truth: R0 = 2.5 mΩ, R1 = 1.5 mΩ, τ1 = 30 s, Q = 54 Ah)

| Method | R0 error | R1 error | τ1 error | fit / val RMSE [mV] |
|---|---|---|---|---|
| NLS 1RC | −0.07 % | +1.0 % | +1.3 % | 10.04 / 10.05 |
| LS-separable 1RC | −0.07 % | +1.0 % | +1.3 % | 10.04 / 10.05 |
| RLS fixed-τ 1RC | −0.02 % | +0.9 % | +1.3 % | 10.04 / 10.05 |
| LS-ARX (1 s avg) | +0.9 % | −22 % | −31 % | 11.92 / 10.52 |
| RLS-ARX (1 s avg) | +0.8 % | −23 % | −31 % | 11.96 / 10.54 |
| LS-ARX (10 Hz) | −2.5 % | −74 % | −98 % | 27.49 / 16.11 |

NLS also finds Q = 54.08 Ah (+0.15 %) and SOC0 = 0.5000. Model selection
correctly picks the 1RC model. **Conclusion:** NLS, separable LS and fixed-τ
RLS are reliable. The ARX methods are biased by noise, as expected.

### 5b. Real data (`ECM1.mat`)

| Quantity | Value |
|---|---|
| **R0** | **2.501 mΩ** (NLS std ±0.002 mΩ; all methods 2.46–2.51 mΩ) |
| RC branch | **not supported by the data** (see below) |
| OCV model | **10th-order polynomial of the table** (chosen by BIC) |
| Q | **59.9 Ah** (stable: 59.9–60.1 Ah over different windows) |
| SOC0 | 0.500 |
| Voltage noise (model-free estimate) | 9.98 mV |
| Fit RMSE / validation RMSE | 9.85 / 10.02 mV (both = noise level, white residuals) |
| Checks | all OK for the recommended model |

**Why no RC branch?**

- NLS 1RC and 2RC push R1 and τ1 to their range limits (R1 → 0.1 mΩ,
  τ1 → 200 s), and the BIC is better for the R0-only model. PSO finds the same.
- In the separable-LS τ-scan the fit error changes by less than 0.001 mV over
  τ = 1 to 200 s, so the data cannot tell the time constant.
- The RC branch voltage (about 1.4 mV RMS) is far below the 10 mV noise.
- On synthetic data the same code finds the RC branch, so this is a property
  of the data, not a code bug.

**How the validation problem was solved (a good example of review work)**

- First version (linear OCV interpolation): validation RMSE 15.6 mV,
  a slow drift up to +20 mV at low SOC, and Q depended on the fit window
  (51–56 Ah). The checks flagged this as "something is not modelled".
- The MATLAB scripts use a 10th-order polynomial OCV, Q = 60 Ah and
  SOC0 = 0.5. With the polynomial OCV the fit finds Q = 59.9 Ah and
  SOC0 = 0.500 by itself, the drift disappears, and validation RMSE = noise.
- So the data was most likely generated with the polynomial OCV. The
  polynomial differs from the table by up to 15 mV, which is exactly the
  size of the drift.
- Lesson: the OCV model is part of the model and must be validated like the
  parameters. On real cell data the measured table may be the better choice,
  so the pipeline compares both and does not assume either one.

**Final model for this dataset:** R0-only ECM, R0 = 2.50 mΩ, polynomial
OCV, Q ≈ 60 Ah, SOC0 = 0.50.

All numbers: `results/real/summary.md` and `results/real/params.json`.

## 6. MATLAB reference scripts → Python

The course scripts `Exercise3.m`, `Exercise4.m` and `PSO_battery_parameter.m`
were used as references. What was taken over, and how it was checked:

| MATLAB | Python | Result on ECM1 / synthetic |
|---|---|---|
| OCV `polyfit(…,10)` | `ocv.PolyOCV`, OCV-model selection | explains ECM1 (val RMSE 15.6 → 10.0 mV) |
| LS for R0 (Method 1) | `separable_ls(n_rc=0)` | 2.501 mΩ |
| Difference method ΔV = R0·ΔI (Method 2) | `difference_ls` | 2.510 mΩ, needs no SOC/OCV |
| RLS for R0, λ and P0 studies | `r0_rls_sensitivity` | λ = 0.97 → R0 std 0.18 mΩ; λ = 1 → 0.001 mΩ |
| 1RC ARX LS / RLS (Ex. 4) | `arx_ls`, `arx_rls` | biased by noise (τ 31–98 % too low) |
| 2RC ARX RLS (Ex. 3 end) | `arx2_ls`, `arx2_rls` (exact pole mapping) | exact without noise; non-physical with 10 mV noise |
| SOC-segmented LSQ (Ex. 4) | `segmented_ls`, `simulate_varying` | synthetic R0(SOC) within 1.5 %; ECM1 R0 flat (2.50 ± 0.01) |
| PSO for R0, R1, C1 | `pso.pso_fit` (log scale, τ, seeds, polish) | synthetic within 3.5 % (1.3 % with polish); ECM1 = NLS |
| PSO with the MATLAB settings | `pso.matlab_style_pso` | swarm gets stuck on the lower bounds (RMSE 25–46 mV) |
| CSV export for ANN | `export.export_timeseries_csv` | `timeseries_ECM1.csv`, columns documented |

Bugs and risky assumptions found in the MATLAB code (wrong signs, wrong R1
formula, wrong RC input gain, SOC assumed linear 0→1, broken lags in
segmented LSQ, clamps that hide failures, PSO bounds, …) are listed with
line numbers in **`notes/matlab_review.md`**.

Full results: `results/reference_methods/summary.md`.

## 7. How to run

```bash
pip install -r requirements.txt
pytest                                           # 62 tests, about 25 s
PYTHONPATH=src python -m ecm.pipeline            # real + synthetic, about 8 s
PYTHONPATH=src python -m ecm.reference_methods   # MATLAB methods + PSO, about 15 s
# pipeline options: --only real | --only synthetic, --ocv linear poly10,
#                   --no-plots, --data <file>, --out <dir>
```

Or `pip install -e .` and then `python -m ecm.pipeline`.

## 8. Project layout

```
data/ECM1.mat            public data (read-only, checksum in CLAUDE.md)
src/ecm/
  config.py              ranges, split time, RLS settings, synthetic truths
  data.py                loading, data checks, noise estimate, SOC sign check
  ocv.py                 OCV(SOC): table interpolation and polynomial
  model.py               ECM simulation (ZOH), also with time-varying parameters
  ls.py                  LS-ARX, ARX2, separable LS, difference LS
  rls.py                 RLS-ARX, ARX2-RLS, fixed-tau RLS, R0+offset RLS, lambda/P0 study
  nls.py                 output-error NLS (+ Q, SOC0), BIC
  pso.py                 particle swarm optimisation (+ MATLAB set-up replica)
  segmented.py           SOC-segmented LS and look-up tables
  checks.py              metrics, residual and plausibility checks
  synthetic.py           synthetic data with known parameters
  export.py              CSV export
  plots.py               figures
  pipeline.py            main pipeline, writes results/real and results/synthetic
  reference_methods.py   MATLAB methods study, writes results/reference_methods
tests/                   pytest tests
results/                 real/, synthetic/, reference_methods/
notes/agent_log.md       agent log and open questions
notes/matlab_review.md   review of the MATLAB reference scripts
.claude/skills/          reusable skill: fit + plausibility check
CLAUDE.md                rules for the agent
```

## 9. Open questions (for the team)

See `notes/agent_log.md`. Short version:

1. How do you normally get the capacity Q and the OCV curve: separate tests,
   or fitted? Here the OCV model choice changed Q from 51 to 60 Ah.
2. Do you use OCV tables (interpolated) or fitted polynomials/splines in
   your models, and how do you validate them?
3. Do you identify RC branches from pulse tests (HPPC) rather than drive
   cycles? This drive cycle with 10 mV noise does not excite them enough.
4. Which forgetting factor and covariance handling do you use for online RLS?
5. When would you use a global optimizer (PSO, GA) instead of NLS with
   multiple starts?

## 10. Next steps (stretch goals)

- EKF-based SOC estimation using the identified model.
- Repeat on data at other temperatures or aging states and compare R0.
- Instrumental-variable (IV) method to remove the ARX noise bias.
- Run the ported methods on `ECM2.mat` (not available here) and compare
  with the true values stored in that file.
