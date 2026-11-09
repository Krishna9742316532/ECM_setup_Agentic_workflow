# Review of the MATLAB reference scripts

Files: `Exercise3.m` (ECM1), `Exercise4.m` (ECM2), `PSO_battery_parameter.m` (ECM2).
They were used as **reference only**. Every method was re-implemented in
Python and checked with the project's validation rules. Results:
`results/reference_methods/summary.md`.

Limits of this review: I did not run the MATLAB code (no MATLAB here), and
`ECM2.mat` was not available. Problems marked **(tested)** were reproduced
in Python on ECM1 or on synthetic data. The others are found by reading the
code.

Severity: **High** = gives wrong numbers, **Medium** = gives fragile or
misleading results, **Low** = cosmetic or a code smell.

## 1. Most important finding (positive)

The scripts use `Q = 60 Ah`, `SOC0 = 0.5` and a 10th-order polynomial OCV.
With exactly this OCV model, the Python NLS fit finds **Q = 59.9 Ah,
SOC0 = 0.500, R0 = 2.501 mΩ**. The validation error drops from 15.6 mV to
10.0 mV (the noise level), and Q no longer depends on the time window
(59.9–60.1 Ah). **(tested)**

So ECM1 was most likely generated with the polynomial OCV, not with the
table. The project now compares both OCV models and selects one
automatically (BIC). On a real cell the measured table may be the better
choice, so the comparison has to be repeated for every dataset. The script's
assumption should not simply be copied.

## 2. Exercise3.m

| Line(s) | Severity | Problem | Python handling |
|---|---|---|---|
| 5 | Medium | `polyfit(…,10)` on 41 points: no check of the fit. It misses the table by up to 14.7 mV (more than the noise); high orders without centring/scaling can also be badly conditioned in MATLAB. **(tested)** | `PolyOCV` uses scaled `Polynomial.fit`, reports the table error and monotonicity |
| 9–10 | Medium | Q = 60 Ah and SOC0 = 0.5 are assumed, not estimated or checked | NLS fits Q and SOC0; a capacity-window sensitivity check |
| 20–24 | Low | Loop computes `v0 = Voltage(k) - vOC'` (whole vector each time) and `R_not` is overwritten; the result is never used | removed |
| 49, 97–98 | **High** | `V_predicted = vOC(k) - R0*Current`: uses **one scalar** `vOC(k)` (the last value of the loop variable) for all times, and the **minus sign** is wrong for this data (V = OCV + R0·I, negative = discharge) | `simulate()` uses OCV(SOC(t)) and the checked sign |
| 41–46 | Low | Method 2 (difference) gives two estimates, θ(1) and θ(2). Good idea: OCV cancels out. The two should agree, but this is not checked | `difference_ls` returns both and the 1-column LS; R0 = 2.510 mΩ, no SOC/OCV needed **(tested)** |
| 144, 205, 288, 338 | Medium | `Time = 1:length(Y)` overwrites the time vector: plots and anything later use samples, not seconds | time is never overwritten |
| 143, 204 | Medium | λ = 0.97–0.99 at 10 Hz means a memory of only 3–10 s. R0 is very noisy (std 0.04–0.18 mΩ vs 0.001 mΩ for λ = 1) **(tested)** | λ study in section D; default λ = 0.9998 (~500 s) |
| 287 | Medium | P0 = 1 … 1e4 all behave the same, because with currents of ~100 A any P0 above ~1e-3 already means "no trust in the start value". The study cannot show an effect **(tested)** | P0 study uses 1e-8 … 1e4; small P0 converges much more slowly |
| 413, 568 | **High** | `z(k+1) = z(k) - (Ts/Q)*Current(k)`: **flipped sign**. In ECM1, negative current = discharge, so SOC rises during discharge (0.50 → 0.62). Fit RMSE goes up to 160 mV **(tested)** | `check_soc_direction()` catches this; sign checked on the raw data |
| 415, 569 | Medium | SOC clamped to [0, 1]: this hides a wrong sign or a wrong Q instead of showing it | no clamp; out-of-range SOC is reported |
| 432–449 | **High** | Regressor `[Uoc(k-1), Uoc(k-2), I(k), I(k-1), I(k-2)]` has no output lags, so it is not an ARX model of the RC circuit. The mapping `R0 = -k4/(k3+k4)` is marked "example mapping" and is not derived | `arx2_ls` / `arx2_rls` with the exact zero-order-hold pole mapping |
| 518–524 | Medium | A second variant with a constant term but again no output lags; no mapping to physical parameters | not ported |
| 586–632 | Medium | ARX2 with `[-y(k-1), -y(k-2), …]` and a mapping from a paper (continuous-time, approximate). `Delta < 0` (complex poles) is possible and the guard is commented out → complex τ | exact mapping; complex or out-of-range poles give NaN + FAIL |
| 582 | Low | Comment says "Large initial covariance" for `P = 1*eye(5)` | P0 is a parameter |
| 643 | Low | One-step prediction with the final θ is shown as "model output". A one-step prediction always looks good; a free-run simulation is the real test | all methods are evaluated by free-run simulation on a separate validation window |
| — | Medium | No fit/validation split, no residual analysis, no parameter plausibility check | all present in the pipeline |

## 3. Exercise4.m (ECM2, not available)

| Line(s) | Severity | Problem | Python handling |
|---|---|---|---|
| 24 | Medium | `inv(Phi'*Phi)*Phi'*Y` (normal equations): numerically worse than `Phi\Y` | `lstsq` / `nnls` |
| 22–34 | Medium | ARX-LS on noisy data is **biased** (τ much too small). With 10 mV noise on ECM1-type data τ is 31–98 % too low **(tested)** | separable LS / fixed-τ RLS / NLS as unbiased alternatives |
| 38, 141 | Low | `C = (-1/(R1*log(α)))/10`: correct only because Ts = 0.1 s. It should be `-Ts/(R1*log(α))` | ZOH formulas with Ts |
| 283, 451 | **High** | `v1(k+1) = α·v1(k) + (1-α)·I(k)·Ts/C`: wrong input gain. The ZOH form is `R1·(1-α)·I(k)`. `(1-α)·Ts/C` is smaller by the factor Ts/τ | single `rc_voltage()` used everywhere |
| 292–303 | Medium | "Measured Voltage" in the legend and RMSE are computed against `V_pred` (an earlier model), not against the measurement | metrics always against measured V |
| 363 | **High** | `vOC = polyval(Pocv, linspace(0, 1, N))`: assumes SOC runs linearly from 0 to 100 % during the test. Not physical | Coulomb counting |
| 368, 433 | **High** | Regressor uses `Voltage(k-1)` instead of `Y(k-1)` (OCV not removed), so α also has to model the OCV; the mapping to R1 and C is then invalid | regressors always on y = V − OCV |
| 633 | **High** | `Y = vOC - Voltage` (sign flipped) with the formulas for `Y = V - vOC` → R0 negative | one sign convention, checked by tests |
| 669 | **High** | `R1 = (θ3 − R0·(1−α))/(1−α)`: wrong formula; correct is `(θ3 + R0·α)/(1−α)` | `arx_to_physical`, tested by round-trip |
| 682–684, 983–985 | Medium | Clamping R0, R1 to [1e-6, 0.1] Ω and **C to [1, 5000] F**. For a 50–60 Ah cell C1 is typically 1e4–1e5 F, so the clamp cuts real values. Clamps also hide failures | no clamps; out-of-range values are reported as FAIL |
| 618, 741, 897, 1050 | Medium | Initial SOC changes between sections (1.0, 0.7) | one SOC set-up, fitted |
| 791–807, 948–964 | **High** | Segmented LSQ takes `idx = find(SOC in bin)` and builds `Y(k-1)` from the previous **index in the bin**. When SOC moves back and forth over a bin edge (regen), that sample is not the previous one in time → wrong regression | `segmented_ls`: no lags inside bins; RC state from the full profile |
| 545, 724, 880, 1036 | Medium | RLS / segmented results are called "ground truth for ANN training". They are biased estimates, not ground truth | CSV header says "estimates, NOT ground truth" |
| 1245, 1272 | Low | `Ro_internal_resistance`, `sim_SOC` are not defined in this script | `export_timeseries_csv` with documented columns |

## 4. PSO_battery_parameter.m (ECM2, not available)

| Line(s) | Severity | Problem | Python handling |
|---|---|---|---|
| 75–102 then 104… | **High** | The local function `objective` is in the middle of the script, but MATLAB needs local functions at the **end** of a script file. As written, the file does not run (the commented copy at the end suggests this was noticed) | separate module `pso.py` |
| 18–19 | **High** | Bounds R0, R1 ∈ [1 mΩ, 0.5 Ω] on a **linear** scale: almost all particles start near 0.25 Ω, 100× too large. C1 ≤ 5000 F limits τ = R1·C1 to a few seconds | log-scale search, τ instead of C, bounds from `config.RANGES` |
| 13–15, 33 | Medium | 50 particles × 15 iterations, one seed, no convergence check | 60 iterations, several seeds, convergence plot |
| — | **High** | With this set-up, my Python replica of the swarm **collapses onto the lower bounds** (R0 = R1 = 1 mΩ, C1 = 10 F) on both ECM1 and synthetic data. RMSE is 25–46 mV vs 10 mV noise; 100 iterations or other seeds do not help **(tested in Python; please confirm in MATLAB)** | `matlab_style_pso` reproduces it; a test documents it |
| 94 | Low | RC update uses `Current(k)` instead of `Current(k-1)`: a 1-sample (0.1 s) shift vs the other scripts | one ZOH model everywhere |
| 79, 7 | Medium | SOC0 = 0.7 and Q = 60 Ah fixed inside the objective | Q and SOC0 can be fitted |
| 112–152 | Medium | Lower-bound sweep: if the optimum moves with the bound, the parameter sits on the bound (not identifiable or wrong bounds). The sweep shows this but does not flag it | `at_bound` + plausibility FAIL/WARN |
| 147, 160 | Low | `results(i,j) = gbest_obj % Store RMSE` stores the **MSE**; the heatmap title says RMSE | RMSE in mV everywhere |
| 112–254 | Low | Section 3 duplicates section 2 | one function |

**Our PSO** (log scale, τ, 60 iterations) recovers the synthetic truth within
about 3.5 % without polish and within about 1.3 % with a short LS polish. On
ECM1 it agrees with NLS (R0 = 2.50 mΩ, RC branch on the bounds = not
identifiable). Different seeds give the same result after polish.

## 5. What is good in the scripts

- The difference method (ΔV = R0·ΔI) is a neat OCV-free R0 estimate. It is
  kept as a cross-check.
- The λ and P0 studies are the right questions; they just needed a wider P0
  range and a noise measure.
- The ZOH form `α·v1 + R1·(1−α)·I` in Exercise 4 (lines 48, 149) is correct.
- SOC-dependent parameters are a sensible next step for cells where R0
  changes with SOC (not the case in ECM1: R0 is flat at 2.50 ± 0.01 mΩ).
- PSO as a global search is useful when good start values are unknown.
