# ECM parameter identification – project rules for the AI agent

## Goal
Identify the parameters of an equivalent circuit model (ECM) of a Li-ion cell
from a public drive-cycle dataset, and validate that the result is physically
sound and reproducible. Start with a 1RC (Thevenin) model; try 2RC only if the
1RC residuals show a clear pattern.

## Data: `data/ECM1.mat` (public dataset, do not modify)
SHA-256: c0c2b75d87b5e133bcb6821e860f65a61255b96db4421f008278815563e37097

| Variable  | Shape     | Unit | Notes |
|-----------|-----------|------|-------|
| `Time`    | (20005,)  | s    | 0 to 2000.4 s, uniform dt = 0.1 s |
| `Current` | (20005,)  | A    | -251 to +132 A, dynamic drive cycle |
| `Voltage` | (20005,)  | V    | 3.06 to 4.07 V, single cell |
| `OCV`     | (41, 2)   | -, V | col 0 = SOC (1.0 down to 0.0, step 0.025), col 1 = OCV in V |

Known facts (checked on the raw data):
- **Sign convention: negative current = discharge.** Net charge is about -7.2 Ah.
- Voltage noise is about 12 mV (std) even at low current. An RMSE near 12 mV is
  the practical floor; do not chase lower values (that means overfitting).
- No long rest periods, so R1/C1 cannot be read from relaxation curves.
  Use a global fit over the time series.
- Instant-step estimate: R0 ≈ 2.5 mOhm.

Unknown (must be estimated or fitted; state as assumptions in README):
- Cell capacity Q. Rough estimate from OCV: about 55 Ah.
- Initial SOC. Rough estimate from OCV: about 0.50 (ends near 0.37).
- Temperature (assume constant).

## Model conventions
- Use SI units inside all code: A, V, Ohm, F, s, As. Convert to Ah or mOhm only
  for printing and plots.
- SOC(t) = SOC0 + (1 / (3600 * Q_Ah)) * ∫ I dt   (I negative on discharge → SOC falls)
- 1RC: V = OCV(SOC) + I*R0 + V1,  dV1/dt = -V1/(R1*C1) + I/C1
  (with I negative on discharge, V drops below OCV, as expected)
- Discretize the RC branch with the exact zero-order-hold form:
  V1[k+1] = exp(-dt/tau) * V1[k] + R1 * (1 - exp(-dt/tau)) * I[k]
- OCV lookup: linear interpolation on the table; SOC must be sorted ascending
  for `np.interp`.

## Plausible parameter ranges (use as optimizer bounds and as test checks)
| Parameter | Range |
|-----------|-------|
| R0        | 0.5 – 10 mOhm |
| R1        | 0.1 – 10 mOhm |
| tau1 = R1*C1 | 1 – 200 s |
| Q         | 30 – 100 Ah (only if fitted) |
| SOC0      | 0.40 – 0.60 (only if fitted) |
All parameters must be positive. A fitted value sitting on a bound is a warning
sign: report it, do not hide it.

## Fit and validation split
- Fit on t < 1400 s, validate on t ≥ 1400 s. Never tune on the validation part.
- Report RMSE and max absolute error for both parts, in mV.
- Save plots to `results/`: measured vs model voltage, residuals, OCV curve.

## Identification methods (see README section 2)
- NLS output error (`src/ecm/nls.py`) is the reference; it also fits Q and SOC0.
- PSO (`src/ecm/pso.py`) is a global search: search R and tau on a log scale,
  use several seeds, and polish with LS. Check that no parameter is on a bound.
- LS-ARX / RLS-ARX are biased by voltage noise. Never report their tau or R1
  as final without comparing with NLS or separable LS.
- Separable LS and fixed-tau RLS are unbiased; they need Q and SOC0 from NLS.
- Every method must pass the synthetic recovery tests before its real-data
  result is trusted.

## Known results on ECM1.mat (update if the data or code changes)
- R0 = 2.50 mOhm, very robust across all methods.
- OCV model: a 10th-order polynomial of the table explains the data
  (Q = 59.9 Ah stable over windows, SOC0 = 0.500, val RMSE = noise).
  Linear table interpolation leaves a drift (val RMSE 15.6 mV, Q 51-56 Ah).
  Always compare OCV models; never assume one.
- The RC branch is NOT identifiable from this data (RC voltage << noise,
  flat tau-scan, NLS/PSO push R1/tau1 to the bounds). Recommended: R0-only.

## MATLAB reference scripts
- Exercise3.m, Exercise4.m, PSO_battery_parameter.m are references only.
  Do not copy their assumptions (Q = 60 Ah, SOC0, polyfit, bounds, clamps)
  without checking them. Known bugs: `notes/matlab_review.md`.

## Commands
- Tests: `pytest` (about 25 s)
- Full run: `PYTHONPATH=src python -m ecm.pipeline` -> `results/real`, `results/synthetic`
- MATLAB methods + PSO study: `PYTHONPATH=src python -m ecm.reference_methods` -> `results/reference_methods`
- Repeated task "fit and check plausibility": `.claude/skills/ecm-fit-check/SKILL.md`

## How the agent must work
1. **Plan before code.** For any new step, first write a short plan and wait for
   approval.
2. **Small steps.** One module or function per change. Show the diff.
3. **Tests with every change.** Put them in `tests/`, run with `pytest`.
4. **Reproducible.** No randomness without a fixed seed. Use the pinned versions
   in `requirements.txt`. Do not add new dependencies without asking.
5. **Do not edit the data file.** Do not change these rules without asking.
6. Log important decisions and open questions in `notes/agent_log.md`.

## Project layout
```
data/      raw data (read-only)
src/ecm/   model, data loading, LS / RLS / NLS, checks, plots, pipeline
tests/     pytest tests
results/   generated plots and parameter files (real/ and synthetic/)
notes/     agent log and open questions
.claude/skills/ecm-fit-check/   reusable fit + plausibility skill
```
