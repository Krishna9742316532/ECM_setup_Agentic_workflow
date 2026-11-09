"""Recursive least squares (RLS) with a forgetting factor.

RLS updates the parameter estimate at every sample, so it can run online
(for example in a BMS) and it can follow slow parameter changes. The
forgetting factor lam < 1 gives old samples less weight; the effective
memory is about 1/(1-lam) samples.

Two RLS versions, matching the two LS versions in ls.py:
- `arx_rls`: recursive ARX (same noise bias as LS-ARX).
- `fixed_tau_rls`: tau is fixed (for example from separable LS) and RLS
  tracks R0, R1 and an OCV offset c. The regressors are noise-free, so the
  estimates are unbiased. The offset c absorbs small SOC/OCV errors.
"""

import numpy as np

from .ls import arx_regressors, arx_to_physical, block_average
from .model import ECMParams, FitResult, rc_voltage


def rls(Phi, y, lam=1.0, theta0=None, P0=1e2):
    """Plain RLS. Returns the parameter history, shape (N, p)."""
    N, p = Phi.shape
    theta = np.zeros(p) if theta0 is None else np.asarray(theta0, dtype=float).copy()
    P = np.eye(p) * P0 if np.isscalar(P0) else np.asarray(P0, dtype=float).copy()
    hist = np.empty((N, p))
    for k in range(N):
        phi = Phi[k]
        Pphi = P @ phi
        gain = Pphi / (lam + phi @ Pphi)
        theta = theta + gain * (y[k] - phi @ theta)
        P = (P - np.outer(gain, Pphi)) / lam
        P = 0.5 * (P + P.T)  # keep P symmetric (numerical safety)
        hist[k] = theta
    return hist


def arx_rls(t, y, I, lam=0.9998, block=10, t_report=None, Q_Ah=np.nan, soc0=np.nan, P0=1e2):
    """Recursive ARX. lam is per original sample; adjusted for the block.

    P0 is the initial covariance. A finite P0 acts like a weak prior on the
    start values, which matters most when the pole a is close to 1.
    """
    h = float(np.median(np.diff(t))) * block
    yb, Ib = block_average(y, block), block_average(I, block)
    tb = block_average(t, block)
    Phi, target = arx_regressors(yb, Ib)
    hist = rls(Phi, target, lam=lam ** block, theta0=[0.9, 1e-3, -1e-3], P0=P0)
    th = tb[1:]
    phys = np.array([arx_to_physical(a, b0, b1, h) for a, b0, b1 in hist])
    k = len(th) - 1 if t_report is None else int(np.searchsorted(th, t_report) - 1)
    k = max(k, 0)
    R0, R1, tau = (float(v) for v in phys[k])
    return FitResult(
        method=f"RLS-ARX ({h:g} s avg)",
        params=ECMParams(R0=R0, R=(R1,), tau=(tau,), Q_Ah=Q_Ah, soc0=soc0),
        info={"t": th, "R0": phys[:, 0], "R1": phys[:, 1], "tau1": phys[:, 2],
              "lam": lam, "t_report": float(th[k])},
    )


def r0_rls(t, y, I, lam=0.9998, t_report=None, Q_Ah=np.nan, soc0=np.nan):
    """RLS for R0 and an OCV offset c:  y = R0*I + c.

    Useful when the data does not support an RC branch. The offset c shows
    how far the measured voltage drifts away from OCV(SOC) over time.
    """
    Phi = np.column_stack([I, np.ones_like(I)])
    hist = rls(Phi, y, lam=lam, theta0=np.zeros(2))
    k = len(t) - 1 if t_report is None else int(np.searchsorted(t, t_report) - 1)
    k = max(k, 0)
    return FitResult(
        method="RLS R0 + offset",
        params=ECMParams(R0=float(hist[k, 0]), Q_Ah=Q_Ah, soc0=soc0),
        info={"t": t, "R0": hist[:, 0], "offset": hist[:, 1], "lam": lam,
              "t_report": float(t[k]), "offset_at_report_mV": float(hist[k, 1] * 1e3)},
    )


def fixed_tau_rls(t, y, I, tau, lam=0.9998, offset=True, t_report=None,
                  Q_Ah=np.nan, soc0=np.nan):
    """RLS for R0, R1 (and OCV offset c) with a fixed time constant tau."""
    x = rc_voltage(t, I, 1.0, tau)
    cols = [I, x] + ([np.ones_like(I)] if offset else [])
    Phi = np.column_stack(cols)
    hist = rls(Phi, y, lam=lam, theta0=np.zeros(Phi.shape[1]))
    k = len(t) - 1 if t_report is None else int(np.searchsorted(t, t_report) - 1)
    k = max(k, 0)
    R0, R1 = hist[k, 0], hist[k, 1]
    info = {"t": t, "R0": hist[:, 0], "R1": hist[:, 1], "tau1_fixed": float(tau),
            "lam": lam, "t_report": float(t[k])}
    if offset:
        info["offset"] = hist[:, 2]
        info["offset_at_report_mV"] = float(hist[k, 2] * 1e3)
    return FitResult(
        method="RLS-fixed-tau 1RC",
        params=ECMParams(R0=float(R0), R=(float(R1),), tau=(float(tau),), Q_Ah=Q_Ah, soc0=soc0),
        info=info,
    )


def arx2_rls(t, y, I, lam=0.9998, block=10, t_report=None, Q_Ah=np.nan, soc0=np.nan, P0=1e2):
    """Recursive 2nd-order ARX (2RC). Same noise bias as RLS-ARX."""
    from .ls import arx2_regressors, arx2_to_physical
    h = float(np.median(np.diff(t))) * block
    yb, Ib, tb = block_average(y, block), block_average(I, block), block_average(t, block)
    Phi, target = arx2_regressors(yb, Ib)
    hist = rls(Phi, target, lam=lam ** block, theta0=[1.5, -0.55, 1e-3, -1e-3, 0.0], P0=P0)
    th = tb[2:]
    k = len(th) - 1 if t_report is None else max(int(np.searchsorted(th, t_report) - 1), 0)
    R0, R, tau = arx2_to_physical(hist[k], h)
    return FitResult(f"RLS-ARX2 ({h:g} s avg)",
                     ECMParams(R0=R0, R=R, tau=tau, Q_Ah=Q_Ah, soc0=soc0),
                     {"t": th, "theta": hist, "lam": lam, "t_report": float(th[k])})


def r0_rls_sensitivity(t, y, I, lams=(0.97, 0.98, 0.99, 0.999, 0.9998, 1.0),
                       P0s=(1e-8, 1e-7, 1e-6, 1e-5, 1.0, 1e4), lam_for_P0=0.999,
                       t_report=None, t_settle=100.0, offset=False, R0_start=0.01):
    """Effect of forgetting factor and initial covariance on R0-RLS.

    Mirrors the lambda / P0 studies of the MATLAB exercises. For each lambda:
    R0 at t_report and the standard deviation of R0(t) between t_settle and
    t_report (how noisy the online estimate is). Memory = 1/(1-lam) samples.

    For each P0: the time until R0(t) first comes within 5 % of the
    reference (lambda = 1 estimate). P0 only matters relative to 1/I^2:
    with currents of ~100 A, every P0 >= ~1e-3 already means "no trust in
    the start value", so P0 = 1 ... 1e4 all behave the same.
    """
    t_report = t[-1] if t_report is None else t_report
    cols = [I] + ([np.ones_like(I)] if offset else [])
    Phi = np.column_stack(cols)
    th0 = np.r_[R0_start, np.zeros(Phi.shape[1] - 1)]
    win = (t >= t_settle) & (t < t_report)
    k = max(int(np.searchsorted(t, t_report) - 1), 0)
    dt = float(np.median(np.diff(t)))
    lam_rows, P0_rows, hist_lam, hist_P0 = [], [], {}, {}
    for lam in lams:
        h = rls(Phi, y, lam=lam, theta0=th0, P0=1.0)[:, 0]
        hist_lam[lam] = h
        mem = float("inf") if lam >= 1 else dt / (1 - lam)
        lam_rows.append({"lambda": lam, "memory_s": mem, "R0_mOhm": h[k] * 1e3,
                         "R0_std_mOhm": float(np.std(h[win]) * 1e3)})
    ref = rls(Phi, y, lam=1.0, theta0=th0, P0=1.0)[k, 0]
    for P0 in P0s:
        h = rls(Phi, y, lam=lam_for_P0, theta0=th0, P0=P0)[:, 0]
        hist_P0[P0] = h
        close = np.flatnonzero(np.abs(h[: k + 1] - ref) < 0.05 * abs(ref))
        t_conv = float(t[close[0]] - t[0]) if close.size else float("inf")
        P0_rows.append({"P0": P0, "R0_mOhm": h[k] * 1e3, "t_reach_5pct_s": t_conv})
    return {"lambda": lam_rows, "P0": P0_rows, "hist_lambda": hist_lam, "hist_P0": hist_P0,
            "lam_for_P0": lam_for_P0, "R0_ref_mOhm": ref * 1e3}
