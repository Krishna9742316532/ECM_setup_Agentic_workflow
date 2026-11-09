"""Particle swarm optimisation (PSO) for ECM parameters.

Python version of PSO_battery_parameter.m, with these changes:
- the swarm searches R and tau on a LOG scale (they span decades; a uniform
  search in [1 mOhm, 0.5 Ohm] puts almost no particles near 2.5 mOhm);
- time constant tau instead of C (better conditioned, easy bounds);
- fixed seed, several seeds can be compared (reproducibility);
- optional "polish": a local least-squares step from the PSO result
  (hybrid global + local search);
- the same plausibility checks as all other methods.

PSO is a global, gradient-free search. It is useful when the cost has
several minima or no good start values. It is slower than NLS and gives no
standard errors by itself.
"""

import numpy as np
from scipy.optimize import least_squares

from .model import ECMParams, FitResult, simulate


def pso_minimize(f, lb, ub, n_particles=50, iters=100, w=0.8, c1=1.5, c2=1.5, seed=1,
                 v_max_frac=0.2, v0=None):
    """Minimise f(x) inside the box [lb, ub]. Returns (x_best, f_best, history).

    v0: std of the initial velocity (absolute). Default 10 % of the box.
    v_max_frac: velocity limit as a fraction of the box (None = no limit).
    """
    rng = np.random.default_rng(seed)
    lb, ub = np.asarray(lb, float), np.asarray(ub, float)
    span = ub - lb
    X = lb + span * rng.random((n_particles, len(lb)))
    Vel = (0.1 * span if v0 is None else v0) * rng.standard_normal((n_particles, len(lb)))
    cost = np.array([f(x) for x in X])
    pbest, pbest_cost = X.copy(), cost.copy()
    g = int(np.argmin(cost))
    gbest, gbest_cost = X[g].copy(), float(cost[g])
    hist = [gbest_cost]
    for _ in range(iters):
        r1, r2 = rng.random(X.shape), rng.random(X.shape)
        Vel = w * Vel + c1 * r1 * (pbest - X) + c2 * r2 * (gbest - X)
        if v_max_frac is not None:
            Vel = np.clip(Vel, -v_max_frac * span, v_max_frac * span)
        X = np.clip(X + Vel, lb, ub)
        cost = np.array([f(x) for x in X])
        better = cost < pbest_cost
        pbest[better], pbest_cost[better] = X[better], cost[better]
        g = int(np.argmin(pbest_cost))
        if pbest_cost[g] < gbest_cost:
            gbest, gbest_cost = pbest[g].copy(), float(pbest_cost[g])
        hist.append(gbest_cost)
    return gbest, gbest_cost, np.array(hist)


def _space(n_rc, fit_soc, ranges):
    """Search space: ln(R0), [Q, soc0], then ln(R_i), ln(tau_i)."""
    names, lb, ub = ["lnR0"], [np.log(ranges["R0"][0])], [np.log(ranges["R0"][1])]
    if fit_soc:
        names += ["Q_Ah", "soc0"]
        lb += [ranges["Q_Ah"][0], ranges["soc0"][0]]
        ub += [ranges["Q_Ah"][1], ranges["soc0"][1]]
    for i in range(1, n_rc + 1):
        names += [f"lnR{i}", f"lntau{i}"]
        lb += [np.log(ranges["R"][0]), np.log(ranges["tau"][0])]
        ub += [np.log(ranges["R"][1]), np.log(ranges["tau"][1])]
    return names, np.array(lb), np.array(ub)


def _to_params(x, n_rc, fit_soc, Q_Ah, soc0):
    R0 = np.exp(x[0])
    if fit_soc:
        Q_Ah, soc0, rest = x[1], x[2], x[3:]
    else:
        rest = x[1:]
    R = tuple(float(np.exp(v)) for v in rest[0::2])
    tau = tuple(float(np.exp(v)) for v in rest[1::2])
    return ECMParams(R0=float(R0), R=R, tau=tau, Q_Ah=float(Q_Ah), soc0=float(soc0))


def pso_fit(data, mask=None, n_rc=1, ranges=None, fit_soc=True, Q_Ah=None, soc0=None,
            n_particles=50, iters=100, seed=1, polish=True):
    """Fit an ECM by PSO on the samples in `mask` (cost = voltage MSE)."""
    from .config import RANGES
    ranges = ranges or RANGES
    if not fit_soc and (Q_Ah is None or soc0 is None):
        raise ValueError("Q_Ah and soc0 are needed when fit_soc=False")
    mask = np.ones_like(data.t, dtype=bool) if mask is None else np.asarray(mask, bool)
    last = int(np.flatnonzero(mask)[-1]) + 1
    t, I, V, m = data.t[:last], data.I[:last], data.V[:last], mask[:last]
    names, lb, ub = _space(n_rc, fit_soc, ranges)

    def resid(x):
        v, _ = simulate(t, I, _to_params(x, n_rc, fit_soc, Q_Ah, soc0), data.ocv)
        return (v - V)[m]

    def mse(x):
        r = resid(x)
        return float(r @ r / len(r))

    x, f, hist = pso_minimize(mse, lb, ub, n_particles, iters, seed=seed)
    info = {"names": names, "seed": seed, "iters": iters, "n_particles": n_particles,
            "n_evals": n_particles * (iters + 1), "rmse_pso_mV": float(np.sqrt(f) * 1e3),
            "history_rmse_mV": (np.sqrt(hist) * 1e3).tolist()}
    label = f"PSO {'R0-only' if n_rc == 0 else f'{n_rc}RC'}"
    if polish:
        r = least_squares(resid, x, bounds=(lb, ub), method="trf")
        x = r.x
        info["rmse_polished_mV"] = float(np.sqrt(2 * r.cost / m.sum()) * 1e3)
        label += " + LS polish"
    p = _to_params(x, n_rc, fit_soc, Q_Ah, soc0)
    if n_rc > 1:
        o = np.argsort(p.tau)
        p = p.replace(R=tuple(p.R[i] for i in o), tau=tuple(p.tau[i] for i in o))
    edge = 0.01 * (ub - lb)
    info["at_bound"] = [n for n, xi, lo, hi in zip(names, x, lb, ub)
                        if xi - lo < edge[names.index(n)] or hi - xi < edge[names.index(n)]]
    return FitResult(label, p, info)


# ---------------------------------------------------------------------------
# The original MATLAB set-up, for comparison
# ---------------------------------------------------------------------------

def matlab_style_pso(data, Q_Ah=60.0, soc0=0.5, lb=(1e-3, 1e-3, 10.0), ub=(0.5, 0.5, 5000.0),
                     n_particles=50, iters=15, seed=1):
    """PSO exactly as in PSO_battery_parameter.m: [R0, R1, C1] on a linear
    scale, Q and SOC0 fixed, 50 particles, 15 iterations.

    The MATLAB RC update uses I(k) instead of I(k-1); here the project's
    zero-order-hold model is used so results are comparable with the other
    methods (difference: one 0.1 s sample).
    """
    def mse(x):
        R0, R1, C1 = x
        p = ECMParams(R0=R0, R=(R1,), tau=(R1 * C1,), Q_Ah=Q_Ah, soc0=soc0)
        v, _ = simulate(data.t, data.I, p, data.ocv)
        return float(np.mean((v - data.V) ** 2))

    # MATLAB: V = randn*0.01, no velocity limit, bounds by clipping.
    x, f, hist = pso_minimize(mse, lb, ub, n_particles, iters, seed=seed, v_max_frac=None, v0=0.01)
    R0, R1, C1 = x
    p = ECMParams(R0=float(R0), R=(float(R1),), tau=(float(R1 * C1),), Q_Ah=Q_Ah, soc0=soc0)
    at = [n for n, xi, lo, hi in zip(("R0", "R1", "C1"), x, lb, ub)
          if xi - lo < 0.01 * (hi - lo) or hi - xi < 0.01 * (hi - lo)]
    return FitResult(f"PSO MATLAB set-up ({iters} iter)", p,
                     {"rmse_mV": float(np.sqrt(f) * 1e3), "C1_F": float(C1), "at_bound": at, "seed": seed,
                      "lb": list(lb), "ub": list(ub), "history_rmse_mV": (np.sqrt(hist) * 1e3).tolist()})
