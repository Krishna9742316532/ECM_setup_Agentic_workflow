"""Load the measurement data and run basic sanity checks on it."""

from dataclasses import dataclass

import numpy as np
import scipy.io as sio

from .ocv import OCVCurve


@dataclass
class CellData:
    """One test profile. Sign convention: negative current = discharge."""

    t: np.ndarray      # s
    I: np.ndarray      # A
    V: np.ndarray      # V
    ocv: OCVCurve
    name: str = "data"

    @property
    def dt(self):
        return float(np.median(np.diff(self.t)))

    def mask_before(self, t_split):
        return self.t < t_split

    def head(self, t_split):
        """Copy with only the samples t < t_split (used for fitting)."""
        m = self.mask_before(t_split)
        return CellData(self.t[m], self.I[m], self.V[m], self.ocv, self.name + "[fit]")


def load_mat(path, name=None):
    """Load a .mat file with variables Time, Current, Voltage, OCV."""
    d = sio.loadmat(path, squeeze_me=True)
    missing = [k for k in ("Time", "Current", "Voltage", "OCV") if k not in d]
    if missing:
        raise KeyError(f"missing variables in {path}: {missing}")
    return CellData(
        t=np.asarray(d["Time"], dtype=float).ravel(),
        I=np.asarray(d["Current"], dtype=float).ravel(),
        V=np.asarray(d["Voltage"], dtype=float).ravel(),
        ocv=OCVCurve.from_table(d["OCV"]),
        name=name or str(path).split("/")[-1],
    )


def estimate_noise_std(V, I, dI_max=0.5):
    """Model-free estimate of voltage noise std [V].

    Uses sample-to-sample voltage differences where the current hardly
    changes. OCV and RC voltages change slowly, so these differences are
    mostly noise. For white noise, std(diff) = sqrt(2) * std(noise).
    A robust (median-based) std is used so a few steps do not bias it.
    """
    dV = np.diff(V)
    dI = np.diff(I)
    sel = np.abs(dI) < dI_max
    if sel.sum() < 100:
        sel = np.ones_like(dV, dtype=bool)
    x = dV[sel]
    mad = np.median(np.abs(x - np.median(x)))
    return float(1.4826 * mad / np.sqrt(2.0))


def check_data(d: CellData, v_min=2.0, v_max=4.5):
    """Return a list of problems found in the data (empty list = all good).

    Checks: shapes, finite values, time order, uniform sampling, voltage
    range, OCV table, and the current sign convention.
    """
    issues = []
    n = len(d.t)
    if not (len(d.I) == n and len(d.V) == n):
        issues.append("Time, Current and Voltage have different lengths")
        return issues
    for k, x in (("Time", d.t), ("Current", d.I), ("Voltage", d.V)):
        if not np.all(np.isfinite(x)):
            issues.append(f"{k} has NaN or Inf values")
    dt = np.diff(d.t)
    if np.any(dt <= 0):
        issues.append("Time is not strictly increasing")
    elif np.max(np.abs(dt - np.median(dt))) > 1e-6 * max(1.0, np.median(dt)) + 1e-9:
        issues.append("Sampling is not uniform")
    if d.V.min() < v_min or d.V.max() > v_max:
        issues.append(f"Voltage outside [{v_min}, {v_max}] V - wrong unit (mV?) or not a single cell")
    if d.ocv.soc[0] < 0 or d.ocv.soc[-1] > 1:
        issues.append("OCV table SOC is outside [0, 1] - SOC in % instead of fraction?")
    if np.any(np.diff(d.ocv.ocv) <= 0):
        issues.append("OCV is not strictly increasing with SOC")
    # Sign convention check: with V = OCV + R0*I (negative = discharge),
    # a step up in current gives a step up in voltage -> positive slope.
    dI, dV = np.diff(d.I), np.diff(d.V)
    big = np.abs(dI) > 0.2 * np.max(np.abs(dI))
    if big.sum() >= 3:
        slope = np.polyfit(dI[big], dV[big], 1)[0]
        if slope <= 0:
            issues.append("dV/dI at current steps is negative - current sign is probably positive = discharge")
    return issues


def current_step_resistance(d: CellData, frac=0.2):
    """Quick R0 estimate [Ohm] from the voltage jump at large current steps."""
    dI, dV = np.diff(d.I), np.diff(d.V)
    big = np.abs(dI) > frac * np.max(np.abs(dI))
    return float(np.polyfit(dI[big], dV[big], 1)[0])


def check_soc_direction(I, soc):
    """SOC must rise with positive (charge) current and fall on discharge.

    Returns an empty list if OK. Catches a flipped sign in Coulomb counting
    (e.g. z(k+1) = z(k) - Ts/Q*I(k) when negative current = discharge).
    """
    ds = np.diff(soc)
    if np.allclose(ds, 0):
        return ["SOC does not change"]
    c = np.corrcoef(ds, np.asarray(I)[:-1])[0, 1]
    return [] if c > 0 else ["SOC rises during discharge - Coulomb counting sign is flipped"]
