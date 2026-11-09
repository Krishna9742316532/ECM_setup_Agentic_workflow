# Agent log

Short notes on each step: what I asked, what the agent proposed or did, what
I checked or changed, and open questions for the team.

## Step 1 – Inspect the data and set up the project (2026-11-09)
- **Asked:** check the data, then set up the folder, rules file (CLAUDE.md) and
  pinned versions.
- **Agent did:** inspected `ECM1.mat` before writing code: 20005 samples at
  10 Hz, current −251 to +132 A, voltage 3.06–4.07 V, 41-point OCV table.
  Found: negative current = discharge; noise about 10–12 mV; R0 ≈ 2.5 mΩ from
  current steps; no long rests, so a global fit is needed; Q and SOC0 unknown
  (rough guess 55 Ah, 0.50).
- **Review:** facts checked against the plots. Data confirmed as public.

## Step 2 – Explore methods before building (2026-11-09)
- **Asked:** use LS and RLS (and other methods) in Python, build the full project.
- **Agent did quick experiments and found:**
  1. Classic ARX-LS is strongly biased by the 10 mV voltage noise (on synthetic
     data τ = 30 s came out as 0.6 s). Block averaging helps but does not fix it.
  2. On the real data, adding an RC branch hardly reduces the error; NLS
     pushes τ to its upper bound and trades it against Q.
  3. Q and SOC0 are correlated and only weakly fixed by one short profile.
- **Decision:** keep ARX-LS/RLS (classic, and a good lesson), and add unbiased
  versions: separable LS (τ grid + linear LS) and fixed-τ RLS. Use NLS as the
  reference. Test everything on synthetic data with known parameters first.

## Step 3 – Code and tests (2026-11-09)
- Modules: `data`, `ocv`, `model`, `ls`, `rls`, `nls`, `checks`, `synthetic`,
  `plots`, `pipeline`. 40 tests, all passing.
- **Issues found and fixed during review:**
  - `tests/__init__.py` broke the `conftest` import → removed, pytest paths
    set in `pyproject.toml`.
  - Noise-free RLS-ARX test failed (R1 off by 2 %): a finite initial
    covariance P0 acts like a prior when the pole is close to 1. Made P0 a
    parameter; the test uses a large P0.
  - NLS 2RC: after sorting the branches by τ, the standard errors were still
    in the old order → fixed with a rename map.
  - Model selection first allowed parameters *near* a bound → rule tightened:
    only models with all parameters clearly inside their ranges.
  - τ-scan plot used an axis offset that hid the flatness → now plots the
    RMSE rise above the best value and prints a "flat curve" note.
  - RLS on the real data gave only "FAIL" rows (no RC) → added LS R0-only and
    RLS R0 + offset, which are the right tools when no RC is present.
- **My review:** _TODO: write what you checked yourself (code, plots, numbers)._

## Step 4 – Results (2026-11-09)
- Synthetic: NLS, separable LS and fixed-τ RLS recover R0, R1, τ1 within
  about 1.5 %; ARX methods underestimate τ by 31–98 % (noise bias).
- Real: R0 = 2.50 mΩ (very robust). RC branch not identifiable. Validation
  RMSE 15.6 mV vs noise 10 mV, caused by a slow drift at low SOC; Q depends
  on the fit window (51–56 Ah). *(Explained and fixed in step 5: OCV model.)*

## Step 5 – MATLAB reference scripts (2026-11-09)
- **Asked:** use `Exercise3.m` and `Exercise4.m` as references and set the
  methods up in the Python project; later also `PSO_battery_parameter.m`.
- **Plan (sent before coding):** port the useful methods as new modules, test
  them on synthetic data with known values, and review the MATLAB code.
- **Agent did:** `PolyOCV`, `difference_ls`, `arx2_ls`/`arx2_rls` (exact pole
  mapping), `r0_rls_sensitivity` (λ, P0), `segmented_ls` + `simulate_varying`,
  CSV export, `pso.py` (own PSO + a replica of the MATLAB set-up), and the
  runner `ecm.reference_methods`. 22 new tests (62 in total).
- **Key finding:** with the 10th-order polynomial OCV used in the scripts,
  NLS finds Q = 59.9 Ah and SOC0 = 0.500 by itself and the validation error
  drops to the noise level. The pipeline now selects the OCV model together
  with the number of RC branches (BIC). The earlier "open question" about Q
  and the low-SOC drift is explained.
- **Other findings:** flipped SOC sign in the 2RC part of Exercise3 (SOC rises
  on discharge); wrong R1 formula and sign flip in Exercise4; wrong RC input
  gain `(1-α)·Ts/C`; segmented LSQ breaks the time order of lags; the MATLAB
  PSO set-up (linear scale, C1 ≤ 5000 F) collapses onto the lower bounds in
  my Python replica. Full list: `notes/matlab_review.md`.
- **Issues fixed during the work:** the P0 study first showed no effect →
  P0 range extended to 1e-8 (scale is set by 1/I²); PSO legend repeated per
  seed → one entry; markdown table rows for skipped bins fixed.
- **Not verified:** I could not run MATLAB or use `ECM2.mat`. The PSO collapse
  is reproduced in Python only.
- **My review:** _TODO_

## Open questions for the team
1. How are the capacity Q and the OCV curve normally obtained: separate
   tests, or fitted? Here the OCV model choice changed Q from 51 to 60 Ah.
2. Table interpolation or fitted polynomial/spline for OCV in your models?
   How do you validate the OCV model? Is OCV hysteresis expected?
3. Are RC parameters usually identified from HPPC pulse tests rather than
   drive cycles? This profile with 10 mV noise does not excite them enough.
4. What forgetting factor and covariance handling do you use for online RLS
   (for example covariance reset, or limits on P during rest)?
5. Should the fit/validation split be by time window, or by a separate
   profile (another drive cycle or temperature)?
6. How much should the agent do alone versus stop for approval at each step?
7. When do you use a global optimizer (PSO, GA) instead of NLS with several
   start points?
8. Were the MATLAB course scripts' bugs (see `matlab_review.md`) known, or
   should they be fixed in the course material?
