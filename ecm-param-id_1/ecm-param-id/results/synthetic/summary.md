# Results: synthetic

## Data

- 20005 samples, dt = 0.1 s, duration 2000.4 s
- Current -251.428 to 132.311 A (negative = discharge), net charge -7.22 Ah
- Voltage 2.992 to 3.989 V
- Estimated voltage noise: **9.85 mV** (model-free estimate)
- R0 from current steps: 2.485 mOhm
- Data checks: all passed

## Model selection (NLS, fitted on t < split)

| OCV model | Model | params | Q [Ah] | SOC0 | fit RMSE [mV] | BIC | params on bound | supported |
|---|---|---|---|---|---|---|---|---|
| linear | NLS R0-only | 3 | 49.93 | 0.4900 | 21.850 | -107031.0 | - | yes |
| linear | NLS 1RC | 5 | 54.08 | 0.5000 | 10.042 | -128781.0 | - | yes |
| linear | NLS 2RC | 7 | 54.12 | 0.5000 | 10.041 | -128763.0 | - | no |
| poly10 | NLS R0-only | 3 | 55.65 | 0.4876 | 21.292 | -107755.0 | - | yes |
| poly10 | NLS 1RC | 5 | 64.04 | 0.4952 | 10.059 | -128732.0 | - | yes |
| poly10 | NLS 2RC | 7 | 64.13 | 0.4953 | 10.058 | -128716.0 | R1 | no |

Recommended: **NLS 1RC** with OCV model **linear** (lowest BIC among (OCV model, RC branches) whose parameters are all inside their ranges and not on a bound).

## SOC set-up used by LS and RLS

Q = 54.08 Ah, SOC0 = 0.5000 (from NLS 1RC, OCV linear); SOC range 0.367 to 0.500.

## All methods

| Method | R0 [mOhm] | R1 [mOhm] | tau1 [s] | R2 [mOhm] | tau2 [s] | Q [Ah] | SOC0 | fit RMSE [mV] | val RMSE [mV] | check |
|---|---|---|---|---|---|---|---|---|---|---|
| NLS R0-only | 2.658 | - | - | - | - | 49.93 | 0.4900 | 21.85 | 26.60 | WARN |
| NLS 1RC | 2.498 | 1.515 | 30.39 | - | - | 54.08 | 0.5000 | 10.04 | 10.05 | OK |
| NLS 2RC | 2.497 | 0.100 | 18.24 | 1.422 | 31.69 | 54.12 | 0.5000 | 10.04 | 10.05 | WARN |
| LS R0-only | 2.803 | - | - | - | - | 54.08 | 0.5000 | 28.72 | 16.64 | WARN |
| LS-ARX (10 Hz) | 2.437 | 0.397 | 0.75 | - | - | 54.08 | 0.5000 | 27.49 | 16.11 | FAIL |
| LS-ARX (1 s avg) | 2.522 | 1.166 | 20.56 | - | - | 54.08 | 0.5000 | 11.92 | 10.52 | WARN |
| LS-separable 1RC | 2.498 | 1.515 | 30.39 | - | - | 54.08 | 0.5000 | 10.04 | 10.05 | OK |
| LS-separable 2RC | 2.497 | 0.224 | 22.33 | 1.294 | 32.18 | 54.08 | 0.5000 | 10.04 | 10.05 | WARN |
| RLS R0 + offset | 2.633 | - | - | - | - | 54.08 | 0.5000 | 29.69 | 17.12 | WARN |
| RLS-ARX (1 s avg) | 2.521 | 1.159 | 20.70 | - | - | 54.08 | 0.5000 | 11.96 | 10.54 | WARN |
| RLS-fixed-tau 1RC | 2.499 | 1.513 | 30.39 | - | - | 54.08 | 0.5000 | 10.04 | 10.05 | OK |

LS and RLS methods use the Q and SOC0 from the SOC set-up above.

## Known truth (synthetic data)

R0_mOhm = 2.5, R1_mOhm = 1.5, tau1_s = 30, C1_F = 20000, Q_Ah = 54, soc0 = 0.5

| Method | R0 err [%] | R1 err [%] | tau1 err [%] | Q err [%] |
|---|---|---|---|---|
| NLS R0-only | +6.33 | - | - | -7.53 |
| NLS 1RC | -0.07 | +1.0 | +1.3 | +0.15 |
| NLS 2RC | -0.11 | -5.2 | +5.6 | +0.22 |
| LS R0-only | +12.11 | - | - | - |
| LS-ARX (10 Hz) | -2.52 | -73.5 | -97.5 | - |
| LS-ARX (1 s avg) | +0.87 | -22.2 | -31.5 | - |
| LS-separable 1RC | -0.07 | +1.0 | +1.3 | - |
| LS-separable 2RC | -0.11 | -13.7 | +7.3 | - |
| RLS R0 + offset | +5.30 | - | - | - |
| RLS-ARX (1 s avg) | +0.84 | -22.7 | -31.0 | - |
| RLS-fixed-tau 1RC | -0.02 | +0.9 | +1.3 | - |

## Findings per method

**NLS R0-only** (WARN)
- [OK] all parameters positive and inside plausible ranges
- [WARN] validation RMSE 26.60 mV is clearly above the noise (9.8 mV) - something is not modelled
- [WARN] slow residual drift 18.39 mV (white noise would give ~0.98 mV) - check SOC/capacity or slow dynamics

**NLS 1RC** (OK)
- [OK] all parameters positive and inside plausible ranges
- [OK] RC branch 1: RMS voltage 30.09 mV (noise 9.8 mV)
- [OK] errors are at the noise level and residuals look white

**NLS 2RC** (WARN)
- [WARN] R1 = 0.0001 sits on/near a range limit [0.0001, 0.01]
- [WARN] RC branch 1: RMS voltage 2.23 mV is below the noise (9.8 mV) - weakly identifiable
- [OK] RC branch 2: RMS voltage 28.00 mV (noise 9.8 mV)
- [OK] errors are at the noise level and residuals look white

**LS R0-only** (WARN)
- [OK] all parameters positive and inside plausible ranges
- [WARN] validation RMSE 16.64 mV is clearly above the noise (9.8 mV) - something is not modelled
- [WARN] slow residual drift 24.96 mV (white noise would give ~0.98 mV) - check SOC/capacity or slow dynamics

**LS-ARX (10 Hz)** (FAIL)
- [FAIL] tau1 = 0.7539 outside plausible range [1, 200]
- [OK] RC branch 1: RMS voltage 16.26 mV (noise 9.8 mV)
- [WARN] validation RMSE 16.11 mV is clearly above the noise (9.8 mV) - something is not modelled
- [WARN] slow residual drift 24.22 mV (white noise would give ~0.98 mV) - check SOC/capacity or slow dynamics

**LS-ARX (1 s avg)** (WARN)
- [OK] all parameters positive and inside plausible ranges
- [OK] RC branch 1: RMS voltage 25.27 mV (noise 9.8 mV)
- [WARN] slow residual drift 6.43 mV (white noise would give ~0.98 mV) - check SOC/capacity or slow dynamics

**LS-separable 1RC** (OK)
- [OK] all parameters positive and inside plausible ranges
- [OK] RC branch 1: RMS voltage 30.09 mV (noise 9.8 mV)
- [OK] errors are at the noise level and residuals look white

**LS-separable 2RC** (WARN)
- [OK] all parameters positive and inside plausible ranges
- [WARN] RC branch 1: RMS voltage 4.75 mV is below the noise (9.8 mV) - weakly identifiable
- [OK] RC branch 2: RMS voltage 25.40 mV (noise 9.8 mV)
- [OK] errors are at the noise level and residuals look white

**RLS R0 + offset** (WARN)
- [OK] all parameters positive and inside plausible ranges
- [WARN] validation RMSE 17.12 mV is clearly above the noise (9.8 mV) - something is not modelled
- [WARN] slow residual drift 27.42 mV (white noise would give ~0.98 mV) - check SOC/capacity or slow dynamics

**RLS-ARX (1 s avg)** (WARN)
- [OK] all parameters positive and inside plausible ranges
- [OK] RC branch 1: RMS voltage 25.07 mV (noise 9.8 mV)
- [WARN] slow residual drift 6.51 mV (white noise would give ~0.98 mV) - check SOC/capacity or slow dynamics

**RLS-fixed-tau 1RC** (OK)
- [OK] all parameters positive and inside plausible ranges
- [OK] RC branch 1: RMS voltage 30.06 mV (noise 9.8 mV)
- [OK] errors are at the noise level and residuals look white

## Capacity sensitivity (information only)

Same recommended model fitted on different time windows.

| Window [s] | Q [Ah] | Q std [Ah] | SOC0 | R0 [mOhm] | RMSE [mV] |
|---|---|---|---|---|---|
| 0-1400 | 54.08 | 0.11 | 0.5000 | 2.498 | 10.04 |
| 0-2000.4 | 54.16 | 0.06 | 0.5000 | 2.499 | 10.04 |
| 600-2000.4 | 54.11 | 0.15 | 0.5001 | 2.499 | 10.05 |

Versions: {'python': '3.13.18', 'numpy': '2.5.3', 'scipy': '1.18.1'}
