"""SOC-segmented least squares: parameters as a function of SOC.

Python version of the "segmented LSQ" in MATLAB Exercise 4, with two fixes:

1. No ARX lags inside a bin. The MATLAB code takes y(k-1) from the previous
   sample *inside the bin*; when SOC moves back and forth across a bin edge,
   that sample is not the previous one in time, so the regression is wrong.
   Here the regressors are the current I and the RC response x_tau(I),
   computed on the full profile, so time order is always correct.
2. No ARX noise bias (see ls.py). tau is fixed (from a global fit); R0 and
   R1 are fitted per bin with non-negative least squares.

Limit: the RC state x_tau is computed with one tau for all bins, so the
effect of a changing R1 on the RC state history is approximated.
"""

import numpy as np
from scipy.optimize import nnls

from .model import rc_voltage


def segmented_ls(t, y, I, soc, edges, tau=None, mask=None, min_samples=300, min_I_std=5.0):
    """Fit R0 (and R1 if tau is given) in each SOC bin [edges[i], edges[i+1]).

    Bins with too few samples or too little current excitation are skipped.
    Returns a list of dict rows, one per bin.
    """
    mask = np.ones_like(t, dtype=bool) if mask is None else np.asarray(mask, dtype=bool)
    x = rc_voltage(t, I, 1.0, tau) if tau else None
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = mask & (soc >= lo) & (soc < hi)
        row = {"soc_lo": float(lo), "soc_hi": float(hi), "n": int(m.sum()),
               "I_std_A": float(np.std(I[m])) if m.any() else 0.0, "used": False}
        if m.sum() < min_samples or row["I_std_A"] < min_I_std:
            row["reason"] = "too few samples" if m.sum() < min_samples else "too little current"
            rows.append(row)
            continue
        A = I[m][:, None] if x is None else np.column_stack([I[m], x[m]])
        coef, rnorm = nnls(A, y[m])
        row.update(used=True, R0=float(coef[0]), rmse_mV=float(rnorm / np.sqrt(m.sum()) * 1e3))
        if x is not None:
            row["R1"] = float(coef[1])
            row["tau1"] = float(tau)
        rows.append(row)
    return rows


def lookup(rows, soc, key):
    """Per-sample parameter from the bin table (nearest used bin outside)."""
    used = [r for r in rows if r["used"]]
    if not used:
        raise ValueError("no usable SOC bins")
    centers = np.array([(r["soc_lo"] + r["soc_hi"]) / 2 for r in used])
    vals = np.array([r[key] for r in used])
    out = np.empty_like(np.asarray(soc, dtype=float))
    for i, s in enumerate(np.asarray(soc, dtype=float)):
        inside = [j for j, r in enumerate(used) if r["soc_lo"] <= s < r["soc_hi"]]
        out[i] = vals[inside[0]] if inside else vals[np.argmin(np.abs(centers - s))]
    return out
