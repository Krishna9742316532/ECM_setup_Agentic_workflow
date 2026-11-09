"""Batch least-squares (LS) identification.

Two LS formulations are provided:

1. `arx_ls` - the classic discrete-time ARX form of the 1RC model:
       y[k] = a*y[k-1] + b0*I[k] + b1*I[k-1],   y = V - OCV(SOC)
       R0 = b0,  tau = -dt/ln(a),  R1 = (b1 + a*b0)/(1 - a)
   It is fast and linear, but the regressor y[k-1] contains measurement
   noise. This gives a known bias (a and tau are pulled towards zero). Block
   averaging reduces the noise and the bias, but does not remove it.

2. `separable_ls` - for a FIXED time constant the model is linear in R0 and
   R_i, and the regressors depend only on the current (no noise). So
   ordinary LS gives unbiased R0, R_i. We scan tau on a log grid, keep the
   best one and refine it. Non-negative LS (NNLS) keeps all R >= 0.
"""

import numpy as np
from scipy.optimize import minimize_scalar, nnls

from .model import ECMParams, FitResult, rc_voltage


def block_average(x, n):
    """Average consecutive blocks of n samples (drops the incomplete tail)."""
    if n <= 1:
        return np.asarray(x, dtype=float)
    m = len(x) // n
    return np.asarray(x[: m * n], dtype=float).reshape(m, n).mean(axis=1)


def arx_to_physical(a, b0, b1, h):
    """Convert ARX coefficients to (R0, R1, tau). Invalid pole -> NaN."""
    R0 = float(b0)
    if not (0.0 < a < 1.0):
        return R0, float("nan"), float("nan")
    tau = float(-h / np.log(a))
    R1 = float((b1 + a * b0) / (1.0 - a))
    return R0, R1, tau


def arx_regressors(y, I):
    Phi = np.column_stack([y[:-1], I[1:], I[:-1]])
    return Phi, y[1:]


def arx_ls(t, y, I, block=1, Q_Ah=np.nan, soc0=np.nan):
    """Classic ARX least squares for a 1RC model."""
    h = float(np.median(np.diff(t))) * block
    yb, Ib = block_average(y, block), block_average(I, block)
    Phi, target = arx_regressors(yb, Ib)
    theta, *_ = np.linalg.lstsq(Phi, target, rcond=None)
    a, b0, b1 = theta
    R0, R1, tau = arx_to_physical(a, b0, b1, h)
    label = "LS-ARX" + (f" ({h:g} s avg)" if block > 1 else f" ({1 / h:g} Hz)")
    return FitResult(
        method=label,
        params=ECMParams(R0=R0, R=(R1,), tau=(tau,), Q_Ah=Q_Ah, soc0=soc0),
        info={"theta": theta.tolist(), "h_s": h, "block": block},
    )


def _unit_rc(t, I, tau):
    return rc_voltage(t, I, 1.0, tau)


def _solve(A, y):
    coef, rnorm = nnls(A, y)
    return coef, rnorm ** 2


def separable_ls(t, y, I, n_rc=1, tau_range=(1.0, 200.0), n_grid=60,
                 Q_Ah=np.nan, soc0=np.nan):
    """LS for R0 and R_i on a grid of time constants, then refine (1RC)."""
    if n_rc == 0:
        coef, sse = _solve(np.asarray(I, dtype=float)[:, None], y)
        return FitResult("LS R0-only", ECMParams(R0=float(coef[0]), Q_Ah=Q_Ah, soc0=soc0),
                         {"sse": float(sse)})

    taus = np.geomspace(tau_range[0], tau_range[1], n_grid)
    X = {i: _unit_rc(t, I, tau) for i, tau in enumerate(taus)}
    info = {"tau_grid": taus.tolist()}

    if n_rc == 1:
        sse = np.empty(n_grid)
        for i in range(n_grid):
            _, sse[i] = _solve(np.column_stack([I, X[i]]), y)
        i_best = int(np.argmin(sse))
        lo = taus[max(i_best - 1, 0)]
        hi = taus[min(i_best + 1, n_grid - 1)]

        def f(log_tau):
            return _solve(np.column_stack([I, _unit_rc(t, I, np.exp(log_tau))]), y)[1]

        if hi > lo:
            res = minimize_scalar(f, bounds=(np.log(lo), np.log(hi)), method="bounded",
                                  options={"xatol": 1e-4})
            tau_best = float(np.exp(res.x))
        else:
            tau_best = float(taus[i_best])
        coef, sse_best = _solve(np.column_stack([I, _unit_rc(t, I, tau_best)]), y)
        info.update(sse_grid=sse.tolist(), sse=float(sse_best))
        params = ECMParams(R0=float(coef[0]), R=(float(coef[1]),), tau=(tau_best,),
                           Q_Ah=Q_Ah, soc0=soc0)
        return FitResult("LS-separable 1RC", params, info)

    if n_rc == 2:
        best = (np.inf, None, None)
        for i in range(n_grid):
            for j in range(i + 1, n_grid):
                coef, s = _solve(np.column_stack([I, X[i], X[j]]), y)
                if s < best[0]:
                    best = (s, (i, j), coef)
        s, (i, j), coef = best
        info.update(sse=float(s))
        params = ECMParams(R0=float(coef[0]), R=(float(coef[1]), float(coef[2])),
                           tau=(float(taus[i]), float(taus[j])), Q_Ah=Q_Ah, soc0=soc0)
        return FitResult("LS-separable 2RC", params, info)

    raise ValueError("n_rc must be 0, 1 or 2")


# ---------------------------------------------------------------------------
# Difference method (MATLAB Exercise 3, "Method 2")
# ---------------------------------------------------------------------------

def difference_ls(t, V, I, Q_Ah=np.nan, soc0=np.nan):
    """R0 from voltage and current differences:  V[k]-V[k-1] = R0*(I[k]-I[k-1]).

    OCV changes very little in one sample, so it cancels out: this R0 does
    NOT depend on Q, SOC0 or the OCV table. The regressor (current) has no
    noise, so the estimate is unbiased. Small RC voltage changes are ignored.

    The MATLAB two-column form [I[k], -I[k-1]] is also solved; both columns
    should give the same R0 (a useful consistency check).
    """
    dV, dI = np.diff(V), np.diff(I)
    R0 = float(dI @ dV / (dI @ dI))
    Phi2 = np.column_stack([I[1:], -I[:-1]])
    th2, *_ = np.linalg.lstsq(Phi2, dV, rcond=None)
    return FitResult("LS difference dV = R0*dI", ECMParams(R0=R0, Q_Ah=Q_Ah, soc0=soc0),
                     {"two_column_theta": th2.tolist()})


# ---------------------------------------------------------------------------
# Second-order ARX for a 2RC model (MATLAB Exercise 3, last part)
# ---------------------------------------------------------------------------
#   y[k] = a1*y[k-1] + a2*y[k-2] + b0*I[k] + b1*I[k-1] + b2*I[k-2]
# Exact zero-order-hold mapping with poles p_i = exp(-h/tau_i):
#   a1 = p1 + p2,  a2 = -p1*p2,  b0 = R0
#   b1 = -R0*(p1+p2) + R1*(1-p1) + R2*(1-p2)
#   b2 =  R0*p1*p2 - R1*(1-p1)*p2 - R2*(1-p2)*p1

def arx2_from_physical(R0, R, tau, h):
    p1, p2 = np.exp(-h / np.asarray(tau, dtype=float))
    R1, R2 = R
    return np.array([p1 + p2, -p1 * p2, R0,
                     -R0 * (p1 + p2) + R1 * (1 - p1) + R2 * (1 - p2),
                     R0 * p1 * p2 - R1 * (1 - p1) * p2 - R2 * (1 - p2) * p1])


def arx2_to_physical(theta, h):
    """(R0, (R1, R2), (tau1, tau2)) from ARX2 coefficients; NaN if not physical."""
    a1, a2, b0, b1, b2 = theta
    nan = float("nan")
    disc = a1 ** 2 + 4 * a2
    if disc <= 0:
        return float(b0), (nan, nan), (nan, nan)        # complex/double poles
    p = np.sort([(a1 - np.sqrt(disc)) / 2, (a1 + np.sqrt(disc)) / 2])  # fast first
    if not np.all((p > 0) & (p < 1)):
        return float(b0), (nan, nan), (nan, nan)
    A = np.array([[1 - p[0], 1 - p[1]],
                  [-(1 - p[0]) * p[1], -(1 - p[1]) * p[0]]])
    rhs = np.array([b1 + b0 * (p[0] + p[1]), b2 - b0 * p[0] * p[1]])
    R1, R2 = np.linalg.solve(A, rhs)
    tau = -h / np.log(p)
    return float(b0), (float(R1), float(R2)), (float(tau[0]), float(tau[1]))


def arx2_regressors(y, I):
    Phi = np.column_stack([y[1:-1], y[:-2], I[2:], I[1:-1], I[:-2]])
    return Phi, y[2:]


def arx2_ls(t, y, I, block=1, Q_Ah=np.nan, soc0=np.nan):
    """Classic 2nd-order ARX least squares (same noise bias as arx_ls)."""
    h = float(np.median(np.diff(t))) * block
    yb, Ib = block_average(y, block), block_average(I, block)
    Phi, target = arx2_regressors(yb, Ib)
    theta, *_ = np.linalg.lstsq(Phi, target, rcond=None)
    R0, R, tau = arx2_to_physical(theta, h)
    label = "LS-ARX2" + (f" ({h:g} s avg)" if block > 1 else f" ({1 / h:g} Hz)")
    return FitResult(label, ECMParams(R0=R0, R=R, tau=tau, Q_Ah=Q_Ah, soc0=soc0),
                     {"theta": theta.tolist(), "h_s": h, "block": block})
