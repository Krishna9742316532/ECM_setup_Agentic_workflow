# Results: ECM1.mat

## Data

- 20005 samples, dt = 0.1 s, duration 2000.4 s
- Current -251.428 to 132.311 A (negative = discharge), net charge -7.22 Ah
- Voltage 3.059 to 4.067 V
- Estimated voltage noise: **9.98 mV** (model-free estimate)
- R0 from current steps: 2.501 mOhm
- Data checks: all passed

- SHA-256: `c0c2b75d87b5e133bcb6821e860f65a61255b96db4421f008278815563e37097`

## Model selection (NLS, fitted on t < split)

| OCV model | Model | params | Q [Ah] | SOC0 | fit RMSE [mV] | BIC | params on bound | supported |
|---|---|---|---|---|---|---|---|---|
| linear | NLS R0-only | 3 | 51.03 | 0.5053 | 9.899 | -129201.0 | - | yes |
| linear | NLS 1RC | 5 | 51.58 | 0.5056 | 9.909 | -129154.0 | R1, tau1 | no |
| linear | NLS 2RC | 7 | 52.16 | 0.5060 | 9.946 | -129031.0 | R1, tau1, R2, tau2 | no |
| poly10 | NLS R0-only | 3 | 59.92 | 0.5001 | 9.854 | -129327.0 | - | yes |
| poly10 | NLS 1RC | 5 | 60.72 | 0.5003 | 9.865 | -129278.0 | R1, tau1 | no |
| poly10 | NLS 2RC | 7 | 61.54 | 0.5006 | 9.891 | -129185.0 | R1, tau1, R2, tau2 | no |

Recommended: **NLS R0-only** with OCV model **poly10** (lowest BIC among (OCV model, RC branches) whose parameters are all inside their ranges and not on a bound).

## SOC set-up used by LS and RLS

Q = 59.92 Ah, SOC0 = 0.5001 (from NLS R0-only, OCV poly10); SOC range 0.380 to 0.500.

## All methods

| Method | R0 [mOhm] | R1 [mOhm] | tau1 [s] | R2 [mOhm] | tau2 [s] | Q [Ah] | SOC0 | fit RMSE [mV] | val RMSE [mV] | check |
|---|---|---|---|---|---|---|---|---|---|---|
| NLS R0-only | 2.501 | - | - | - | - | 59.92 | 0.5001 | 9.85 | 10.02 | OK |
| NLS 1RC | 2.500 | 0.100 | 200.00 | - | - | 60.72 | 0.5003 | 9.86 | 10.07 | WARN |
| NLS 2RC | 2.499 | 0.100 | 200.00 | 0.100 | 200.00 | 61.54 | 0.5006 | 9.89 | 10.23 | WARN |
| LS R0-only | 2.501 | - | - | - | - | 59.92 | 0.5001 | 9.85 | 10.02 | OK |
| LS-ARX (10 Hz) | 2.457 | 0.045 | 0.03 | - | - | 59.92 | 0.5001 | 9.85 | 10.02 | FAIL |
| LS-ARX (1 s avg) | 2.499 | - | - | - | - | 59.92 | 0.5001 | - | - | FAIL |
| LS-separable 1RC | 2.498 | 0.004 | 1.00 | - | - | 59.92 | 0.5001 | 9.85 | 10.02 | FAIL |
| LS-separable 2RC | 2.498 | 0.003 | 1.00 | 0.001 | 8.96 | 59.92 | 0.5001 | 9.85 | 10.02 | FAIL |
| RLS R0 + offset | 2.499 | - | - | - | - | 59.92 | 0.5001 | 9.85 | 10.02 | OK |
| RLS-ARX (1 s avg) | 2.500 | 0.000 | 0.36 | - | - | 59.92 | 0.5001 | 9.85 | 10.02 | FAIL |
| RLS-fixed-tau 1RC | 2.498 | 0.002 | 1.00 | - | - | 59.92 | 0.5001 | 9.85 | 10.02 | FAIL |

LS and RLS methods use the Q and SOC0 from the SOC set-up above.

## Findings per method

**NLS R0-only** (OK)
- [OK] all parameters positive and inside plausible ranges
- [OK] errors are at the noise level and residuals look white

**NLS 1RC** (WARN)
- [WARN] R1 = 0.0001 sits on/near a range limit [0.0001, 0.01]
- [WARN] tau1 = 200 sits on/near a range limit [1, 200]
- [WARN] RC branch 1: RMS voltage 1.44 mV is below the noise (10.0 mV) - weakly identifiable
- [OK] errors are at the noise level and residuals look white

**NLS 2RC** (WARN)
- [WARN] R1 = 0.0001 sits on/near a range limit [0.0001, 0.01]
- [WARN] tau1 = 200 sits on/near a range limit [1, 200]
- [WARN] R2 = 0.0001 sits on/near a range limit [0.0001, 0.01]
- [WARN] tau2 = 200 sits on/near a range limit [1, 200]
- [WARN] RC time constants are not ordered or are equal
- [WARN] RC branch 1: RMS voltage 1.44 mV is below the noise (10.0 mV) - weakly identifiable
- [WARN] RC branch 2: RMS voltage 1.44 mV is below the noise (10.0 mV) - weakly identifiable
- [OK] errors are at the noise level and residuals look white

**LS R0-only** (OK)
- [OK] all parameters positive and inside plausible ranges
- [OK] errors are at the noise level and residuals look white

**LS-ARX (10 Hz)** (FAIL)
- [FAIL] R1 = 4.456e-05 outside plausible range [0.0001, 0.01]
- [FAIL] tau1 = 0.02604 outside plausible range [1, 200]
- [WARN] RC branch 1: RMS voltage 1.97 mV is below the noise (10.0 mV) - weakly identifiable
- [OK] errors are at the noise level and residuals look white

**LS-ARX (1 s avg)** (FAIL)
- [FAIL] R1 is not a valid number (e.g. ARX pole outside 0..1)
- [FAIL] tau1 is not a valid number (e.g. ARX pole outside 0..1)

**LS-separable 1RC** (FAIL)
- [FAIL] R1 = 3.609e-06 outside plausible range [0.0001, 0.01]
- [WARN] tau1 = 1 sits on/near a range limit [1, 200]
- [WARN] RC branch 1: RMS voltage 0.14 mV is below the noise (10.0 mV) - weakly identifiable
- [OK] errors are at the noise level and residuals look white

**LS-separable 2RC** (FAIL)
- [FAIL] R1 = 3.102e-06 outside plausible range [0.0001, 0.01]
- [WARN] tau1 = 1 sits on/near a range limit [1, 200]
- [FAIL] R2 = 7.207e-07 outside plausible range [0.0001, 0.01]
- [WARN] RC branch 1: RMS voltage 0.12 mV is below the noise (10.0 mV) - weakly identifiable
- [WARN] RC branch 2: RMS voltage 0.02 mV is below the noise (10.0 mV) - weakly identifiable
- [OK] errors are at the noise level and residuals look white

**RLS R0 + offset** (OK)
- [OK] all parameters positive and inside plausible ranges
- [OK] errors are at the noise level and residuals look white

**RLS-ARX (1 s avg)** (FAIL)
- [FAIL] R1 = 8.852e-08 outside plausible range [0.0001, 0.01]
- [FAIL] tau1 = 0.357 outside plausible range [1, 200]
- [WARN] RC branch 1: RMS voltage 0.00 mV is below the noise (10.0 mV) - weakly identifiable
- [OK] errors are at the noise level and residuals look white

**RLS-fixed-tau 1RC** (FAIL)
- [FAIL] R1 = 1.791e-06 outside plausible range [0.0001, 0.01]
- [WARN] tau1 = 1 sits on/near a range limit [1, 200]
- [WARN] RC branch 1: RMS voltage 0.07 mV is below the noise (10.0 mV) - weakly identifiable
- [OK] errors are at the noise level and residuals look white

## Capacity sensitivity (information only)

Same recommended model fitted on different time windows.

| Window [s] | Q [Ah] | Q std [Ah] | SOC0 | R0 [mOhm] | RMSE [mV] |
|---|---|---|---|---|---|
| 0-1400 | 59.92 | 0.11 | 0.5001 | 2.501 | 9.85 |
| 0-2000.4 | 59.97 | 0.08 | 0.5000 | 2.501 | 9.91 |
| 600-2000.4 | 60.05 | 0.14 | 0.4999 | 2.499 | 9.87 |

Versions: {'python': '3.13.18', 'numpy': '2.5.3', 'scipy': '1.18.1'}
