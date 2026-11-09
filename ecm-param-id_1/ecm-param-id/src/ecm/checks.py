"""Validation checks: error metrics, residual analysis and plausibility.

Every check returns plain data (numbers or a list of findings), so it can be
used in tests, in the pipeline and in reports.
"""

from dataclasses import dataclass

import numpy as np

from .model import ECMParams, rc_voltage


@dataclass
class Finding:
    level: str     # "FAIL", "WARN" or "OK"
    message: str

    def __str__(self):
        return f"[{self.level}] {self.message}"


def metrics(V, Vm, mask):
    e = (Vm - V)[mask]
    return {
        "rmse_mV": float(np.sqrt(np.mean(e ** 2)) * 1e3),
        "max_abs_mV": float(np.max(np.abs(e)) * 1e3),
        "mean_mV": float(np.mean(e) * 1e3),
    }


def acf(x, lags=(1, 10, 100)):
    x = np.asarray(x) - np.mean(x)
    d = float(x @ x)
    return {int(k): float(x[:-k] @ x[k:] / d) for k in lags}


def slow_drift_mV(res, window=100):
    """RMS of the moving-average residual [mV].

    For white noise this is about noise_std/sqrt(window). A much larger value
    means a slow, unmodelled effect (for example an SOC/capacity error).
    """
    k = np.ones(window) / window
    return float(np.sqrt(np.mean(np.convolve(res, k, mode="valid") ** 2)) * 1e3)


def _near_edge(x, lo, hi, tol=0.02):
    """True if x is within `tol` (relative, log scale) of a bound."""
    span = np.log(hi) - np.log(lo)
    return (np.log(x) - np.log(lo)) < tol * span or (np.log(hi) - np.log(x)) < tol * span


def plausibility(p: ECMParams, ranges, check_soc=True):
    """Check signs, ranges and bound hits of the parameters."""
    out = []
    vals = [("R0", p.R0, ranges["R0"])]
    for i, (r, t) in enumerate(zip(p.R, p.tau), start=1):
        vals += [(f"R{i}", r, ranges["R"]), (f"tau{i}", t, ranges["tau"])]
    if check_soc:
        vals += [("Q_Ah", p.Q_Ah, ranges["Q_Ah"]), ("soc0", p.soc0, ranges["soc0"])]

    for name, x, (lo, hi) in vals:
        if not np.isfinite(x):
            out.append(Finding("FAIL", f"{name} is not a valid number (e.g. ARX pole outside 0..1)"))
        elif x <= 0:
            out.append(Finding("FAIL", f"{name} = {x:.4g} is not positive"))
        elif not (lo <= x <= hi):
            out.append(Finding("FAIL", f"{name} = {x:.4g} outside plausible range [{lo:g}, {hi:g}]"))
        elif _near_edge(x, lo, hi):
            out.append(Finding("WARN", f"{name} = {x:.4g} sits on/near a range limit [{lo:g}, {hi:g}]"))
    if p.n_rc > 1 and np.all(np.isfinite(p.tau)) and not np.all(np.diff(p.tau) > 0):
        out.append(Finding("WARN", "RC time constants are not ordered or are equal"))
    if not out:
        out.append(Finding("OK", "all parameters positive and inside plausible ranges"))
    return out


def rc_identifiability(p: ECMParams, t, I, noise_std):
    """Warn when an RC branch voltage is small compared with the noise."""
    out = []
    for i, (r, tau) in enumerate(zip(p.R, p.tau), start=1):
        if not (np.isfinite(r) and np.isfinite(tau)):
            continue
        v = rc_voltage(t, I, r, tau)
        rms = float(np.sqrt(np.mean(v ** 2)))
        if rms < noise_std:
            out.append(Finding("WARN", f"RC branch {i}: RMS voltage {rms * 1e3:.2f} mV is below the "
                                       f"noise ({noise_std * 1e3:.1f} mV) - weakly identifiable"))
        else:
            out.append(Finding("OK", f"RC branch {i}: RMS voltage {rms * 1e3:.2f} mV (noise "
                                     f"{noise_std * 1e3:.1f} mV)"))
    return out


def fit_quality(m_fit, m_val, noise_std, drift_mV=None, window=100):
    """Compare the errors with the noise level."""
    out = []
    noise_mV = noise_std * 1e3
    if m_fit["rmse_mV"] < 0.9 * noise_mV:
        out.append(Finding("WARN", f"fit RMSE {m_fit['rmse_mV']:.2f} mV is below the noise "
                                   f"({noise_mV:.1f} mV) - possible overfitting"))
    if m_val["rmse_mV"] > 1.3 * noise_mV:
        out.append(Finding("WARN", f"validation RMSE {m_val['rmse_mV']:.2f} mV is clearly above the "
                                   f"noise ({noise_mV:.1f} mV) - something is not modelled"))
    if drift_mV is not None:
        white = noise_mV / np.sqrt(window)
        if drift_mV > 3 * white:
            out.append(Finding("WARN", f"slow residual drift {drift_mV:.2f} mV (white noise would give "
                                       f"~{white:.2f} mV) - check SOC/capacity or slow dynamics"))
    if not out:
        out.append(Finding("OK", "errors are at the noise level and residuals look white"))
    return out


def worst_level(findings):
    order = {"OK": 0, "WARN": 1, "FAIL": 2}
    return max((f.level for f in findings), key=order.get, default="OK")
