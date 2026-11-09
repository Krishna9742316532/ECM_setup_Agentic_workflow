# MATLAB reference methods - Python results (ECM1)

Methods from `Exercise3.m` / `Exercise4.m`, re-implemented and checked. Review of the MATLAB code: `notes/matlab_review.md`.

## A. OCV model: interpolation vs polynomial (MATLAB polyfit order 10)

- Polynomial error at the table points: max 14.7 mV, RMS 5.3 mV
- In this test's SOC window: max 14.7 mV (noise is about 10 mV)
- Monotonic: True

| OCV model | R0 [mOhm] | Q [Ah] | SOC0 | fit RMSE [mV] | val RMSE [mV] |
|---|---|---|---|---|---|
| linear interpolation | 2.501 | 51.03 | 0.5053 | 9.90 | 15.57 |
| polynomial order 10 | 2.501 | 59.92 | 0.5001 | 9.85 | 10.02 |

OCV model chosen by the main pipeline: **poly10**.

## B. SOC assumptions: MATLAB scripts vs this project

| Set-up | SOC start -> end | SOC direction | LS R0 [mOhm] | fit RMSE [mV] | val RMSE [mV] |
|---|---|---|---|---|---|
| MATLAB Ex3 (poly10, Q=60 Ah, SOC0=0.5) | 0.500 -> 0.380 | OK | 2.501 | 9.85 | 10.02 |
| MATLAB Ex3 2RC part (same, flipped sign z-Ts/Q*I) | 0.500 -> 0.620 | SOC rises during discharge - Coulomb counting sign is flipped | 3.648 | 160.23 | 314.88 |
| this project before (linear OCV, Q and SOC0 fitted) | 0.505 -> 0.364 | OK | 2.501 | 9.90 | 15.57 |
| this project now (poly10 OCV, Q=59.9 Ah, SOC0=0.500, fitted) | 0.500 -> 0.380 | OK | 2.501 | 9.85 | 10.02 |

## C. R0 from the simple methods (fit window)

| Method | R0 [mOhm] | needs Q / SOC0 / OCV |
|---|---|---|
| Method 1: LS on V - OCV(SOC) | 2.5009 | yes |
| Method 2: difference LS dV = R0*dI (1 column) | 2.5102 | no |
| Method 2: MATLAB 2-column form, theta(1) | 2.5107 | no |
| Method 2: MATLAB 2-column form, theta(2) | 2.5097 | no |
| RLS (lambda = 1, P0 = 1) on V - OCV(SOC) | 2.5009 | yes |
| NLS R0-only, poly10 OCV (reference, also fits Q and SOC0) | 2.5009 | no |

## D. RLS for R0: forgetting factor and initial covariance

| lambda | memory [s] | R0 at split [mOhm] | std of R0(t) [mOhm] |
|---|---|---|---|
| 0.97 | 3.3 | 2.5047 | 0.1765 |
| 0.98 | 5.0 | 2.4917 | 0.1039 |
| 0.99 | 10.0 | 2.4925 | 0.0374 |
| 0.999 | 100.0 | 2.4958 | 0.0051 |
| 0.9998 | 500.0 | 2.4994 | 0.0017 |
| 1 | inf | 2.5009 | 0.0011 |

P0 study (lambda = 0.999, start value 10 mOhm, reference R0 = 2.5009 mOhm):

| P0 | R0 at split [mOhm] | time to first reach +-5 % of reference [s] |
|---|---|---|
| 1e-08 | 2.4965 | 753.8 |
| 1e-07 | 2.4958 | 649.9 |
| 1e-06 | 2.4958 | 501.1 |
| 1e-05 | 2.4958 | 92.9 |
| 1 | 2.4958 | 11.0 |
| 10000 | 2.4958 | 0.5 |

## E. ARX methods (1RC and 2RC)

Synthetic 2RC truth: R0_mOhm = 2.5, R1_mOhm = 1, tau1_s = 5, C1_F = 5000, R2_mOhm = 1.5, tau2_s = 60, C2_F = 40000, Q_Ah = 54, soc0 = 0.5

| Data | Method | R0 [mOhm] | R1 [mOhm] | tau1 [s] | R2 [mOhm] | tau2 [s] | physical |
|---|---|---|---|---|---|---|---|
| ECM1 (real) | LS-ARX (10 Hz) | 2.457 | 0.045 | 0.03 | - | - | yes |
| ECM1 (real) | LS-ARX (1 s avg) | 2.499 | - | - | - | - | NO |
| ECM1 (real) | RLS-ARX (1 s avg) | 2.500 | 0.000 | 0.36 | - | - | yes |
| ECM1 (real) | LS-ARX2 (10 Hz) | 2.455 | - | - | - | - | NO |
| ECM1 (real) | LS-ARX2 (1 s avg) | 2.499 | - | - | - | - | NO |
| ECM1 (real) | RLS-ARX2 (1 s avg) | 2.500 | - | - | - | - | NO |
| synthetic 2RC, noise 0 mV | LS-ARX2 (10 Hz) | 2.500 | 1.000 | 5.00 | 1.500 | 60.00 | yes |
| synthetic 2RC, noise 0 mV | LS-ARX2 (1 s avg) | 2.608 | 0.253 | 2.00 | 1.813 | 15.18 | yes |
| synthetic 2RC, noise 10 mV | LS-ARX2 (10 Hz) | 2.518 | - | - | - | - | NO |
| synthetic 2RC, noise 10 mV | LS-ARX2 (1 s avg) | 2.585 | - | - | - | - | NO |
| synthetic 2RC, noise 10 mV | NLS 2RC | 2.495 | 0.997 | 4.93 | 1.527 | 60.47 | yes |

## F. SOC-segmented LS

Real data (R0-only, 2.5 % SOC bins, fit window):

| SOC bin | samples | R0 [mOhm] | RMSE [mV] |
|---|---|---|---|
| 0.350-0.375 | 0 | - | skipped: too few samples |
| 0.375-0.400 | 0 | - | skipped: too few samples |
| 0.400-0.425 | 3463 | 2.4970 | 9.70 |
| 0.425-0.450 | 2938 | 2.4993 | 9.82 |
| 0.450-0.475 | 1762 | 2.5043 | 9.85 |
| 0.475-0.500 | 5798 | 2.5027 | 9.93 |
| 0.500-0.525 | 39 | - | skipped: too few samples |

| Model | fit RMSE [mV] | val RMSE [mV] |
|---|---|---|
| constant R0 | 9.85 | 10.02 |
| R0(SOC) table | 9.85 | 10.02 |

Synthetic data with R0(SOC) (tau fixed at 32.1 s from separable LS):

| SOC bin | R0 [mOhm] | true R0 at bin centre | error [%] | R1 [mOhm] |
|---|---|---|---|---|
| 0.375-0.400 | 3.0189 | 3.0625 | -1.42 | 1.477 |
| 0.400-0.425 | 2.9093 | 2.9375 | -0.96 | 1.566 |
| 0.425-0.450 | 2.8319 | 2.8125 | +0.69 | 1.509 |
| 0.450-0.475 | 2.6774 | 2.6875 | -0.38 | 1.540 |
| 0.475-0.500 | 2.5717 | 2.5625 | +0.36 | 1.530 |

| Model | fit RMSE [mV] | val RMSE [mV] |
|---|---|---|
| constant R0, R1 | 11.97 | 12.40 |
| R0(SOC), R1(SOC) table | 10.18 | 10.29 |

## G. Particle swarm optimisation (PSO)

Synthetic truth: R0_mOhm = 2.5, R1_mOhm = 1.5, tau1_s = 30, C1_F = 20000, Q_Ah = 54, soc0 = 0.5

MATLAB bounds: {'lb [R0, R1, C1]': [0.001, 0.001, 10.0], 'ub [R0, R1, C1]': [0.5, 0.5, 5000.0]}

| Data | Method | seed | R0 [mOhm] | R1 [mOhm] | tau1 [s] | C1 [F] | Q [Ah] | RMSE [mV] | on bound | error vs truth (R0 / R1 / tau1) |
|---|---|---|---|---|---|---|---|---|---|---|
| synthetic 1RC | PSO 1RC | 1 | 2.497 | 1.520 | 30.00 | 19734 | 54.01 | 10.04 | - | -0.1 % / +1.3 % / -0.0 % |
| synthetic 1RC | PSO 1RC | 2 | 2.506 | 1.494 | 31.04 | 20776 | 54.67 | 10.07 | - | +0.2 % / -0.4 % / +3.5 % |
| synthetic 1RC | PSO 1RC | 3 | 2.508 | 1.489 | 30.71 | 20626 | 54.13 | 10.06 | - | +0.3 % / -0.7 % / +2.4 % |
| synthetic 1RC | PSO 1RC + LS polish | 1 | 2.498 | 1.515 | 30.39 | 20061 | 54.08 | 10.04 | - | -0.1 % / +1.0 % / +1.3 % |
| ECM1 real (poly10 OCV) | PSO R0-only + LS polish | 1 | 2.501 | - | - | - | 59.92 | 9.85 | - | - |
| ECM1 real (poly10 OCV) | PSO R0-only + LS polish | 2 | 2.501 | - | - | - | 59.92 | 9.85 | - | - |
| ECM1 real (poly10 OCV) | PSO R0-only + LS polish | 3 | 2.501 | - | - | - | 59.92 | 9.85 | - | - |
| ECM1 real (poly10 OCV) | PSO 1RC + LS polish | 1 | 2.500 | 0.100 | 200.00 | 2000000 | 60.72 | 9.86 | lnR1, lntau1 | - |
| ECM1 real (poly10 OCV) | PSO 1RC + LS polish | 2 | 2.500 | 0.100 | 200.00 | 2000000 | 60.72 | 9.86 | lnR1, lntau1 | - |
| ECM1 real (poly10 OCV) | PSO 1RC + LS polish | 3 | 2.500 | 0.100 | 200.00 | 2000000 | 60.72 | 9.86 | lnR1, lntau1 | - |
| ECM1 real (MATLAB set-up: poly10, Q=60, SOC0=0.5) | PSO MATLAB set-up (15 iter) | 1 | 1.000 | 1.000 | 0.01 | 10 | 60.00 | 24.67 | R0, R1, C1 | - |
| synthetic 1RC (MATLAB set-up, true Q and SOC0) | PSO MATLAB set-up (15 iter) | 1 | 1.000 | 1.000 | 0.01 | 10 | 54.00 | 45.67 | R0, R1, C1 | -60.0 % / -33.3 % / -100.0 % |

## H. CSV export

`timeseries_ECM1.csv` (20005 rows). The header explains every column. Identified parameters are estimates, not ground truth.
