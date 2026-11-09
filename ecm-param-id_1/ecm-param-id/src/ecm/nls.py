"""Nonlinear least squares (NLS), output-error form.

The full model is simulated and the voltage error is minimised with
scipy.optimize.least_squares. This is the reference method: no noise bias,
and it can also estimate the capacity Q and the initial SOC. It needs
bounds and start values, and it can be slow for long data.

The parameter vector is [R0, (Q_Ah, soc0), R1, ln(tau1), R2, ln(tau2), ...].
Time constants are fitted on a log scale because they span decades.
"""


import numpy as np
from scipy.optimize import least_squares

from .model import ECMParams, FitResult, simulate


def _unpack(x, n_rc, fit_soc, Q_fixed, soc0_fixed):
    R0 = x[0]
    if fit_soc:
        Q, s0, rest = x[1], x[2], x[3:]
    else:
        Q, s0, rest = Q_fixed, soc0_fixed, x[1:]
    R = tuple(rest[0::2][:n_rc])
    tau = tuple(np.exp(rest[1::2][:n_rc]))
    return ECMParams(R0=float(R0), R=tuple(map(float, R)), tau=tuple(map(float, tau)),
                     Q_Ah=float(Q), soc0=float(s0))


def param_names(n_rc, fit_soc=True):
    names = ["R0"] + (["Q_Ah", "soc0"] if fit_soc else [])
    for i in range(1, n_rc + 1):
        names += [f"R{i}", f"tau{i}"]
    return names


def nls_fit(data, mask=None, n_rc=1, ranges=None, fit_soc=True, Q_Ah=None, soc0=None,
            tau_starts=None, R0_start=2.5e-3, Q_start=54.0, soc0_start=0.5):
    """Fit an n_rc-branch ECM by output-error NLS on the samples in `mask`.

    The model is always simulated from t[0], so states (SOC, RC voltages)
    are consistent; only the samples in `mask` count in the cost.
    """
    from .config import RANGES
    ranges = ranges or RANGES
    t, I, V = data.t, data.I, data.V
    mask = np.ones_like(t, dtype=bool) if mask is None else np.asarray(mask, dtype=bool)
    last = int(np.flatnonzero(mask)[-1]) + 1
    t_s, I_s, V_s, m_s = t[:last], I[:last], V[:last], mask[:last]

    if not fit_soc and (Q_Ah is None or soc0 is None):
        raise ValueError("Q_Ah and soc0 are needed when fit_soc=False")

    lb, ub, scale = [ranges["R0"][0]], [ranges["R0"][1]], [1e-3]
    if fit_soc:
        lb += [ranges["Q_Ah"][0], ranges["soc0"][0]]
        ub += [ranges["Q_Ah"][1], ranges["soc0"][1]]
        scale += [10.0, 0.1]
    for _ in range(n_rc):
        lb += [ranges["R"][0], np.log(ranges["tau"][0])]
        ub += [ranges["R"][1], np.log(ranges["tau"][1])]
        scale += [1e-3, 1.0]
    lb, ub = np.array(lb), np.array(ub)

    if tau_starts is None:
        tau_starts = {0: [()], 1: [(3,), (10,), (30,), (100,)],
                      2: [(2, 30), (5, 100), (10, 150)]}[n_rc]

    def resid(x):
        p = _unpack(x, n_rc, fit_soc, Q_Ah, soc0)
        v, _ = simulate(t_s, I_s, p, data.ocv)
        return (v - V_s)[m_s]

    best, runs = None, []
    for taus in tau_starts:
        x0 = [R0_start] + ([Q_start, soc0_start] if fit_soc else [])
        for tau0 in taus:
            x0 += [1e-3, np.log(tau0)]
        x0 = np.clip(np.array(x0, dtype=float), lb + 1e-12, ub - 1e-12)
        r = least_squares(resid, x0, bounds=(lb, ub), x_scale=scale, method="trf")
        runs.append({"start_tau": list(taus), "cost": float(r.cost)})
        if best is None or r.cost < best.cost:
            best = r

    r = best
    p = _unpack(r.x, n_rc, fit_soc, Q_Ah, soc0)
    # Keep the RC branches ordered by time constant (fast first).
    order = np.argsort(p.tau) if n_rc > 1 else np.arange(n_rc)
    p = p.replace(R=tuple(p.R[i] for i in order), tau=tuple(p.tau[i] for i in order))
    # Old branch number -> new branch number (for the names below).
    rename = {f"{a}{int(old) + 1}": f"{a}{new + 1}"
              for new, old in enumerate(order) for a in ("R", "tau")}

    n, k = int(m_s.sum()), len(r.x)
    sse = float(2 * r.cost)
    names = param_names(n_rc, fit_soc)
    info = {"sse": sse, "n": n, "k": k, "rmse_mV": float(np.sqrt(sse / n) * 1e3),
            "starts": runs, "names": names,
            "at_bound": [rename.get(nm, nm) for nm, a in zip(names, r.active_mask) if a != 0]}

    # Approximate standard errors and correlations from the Jacobian.
    try:
        J = r.jac
        s2 = sse / max(n - k, 1)
        cov = np.linalg.pinv(J.T @ J) * s2
        sd = np.sqrt(np.clip(np.diag(cov), 0, None))
        # Convert ln(tau) std to tau std (delta method).
        sd_nat = sd.copy()
        for i, nm in enumerate(names):
            if nm.startswith("tau"):
                sd_nat[i] = sd[i] * np.exp(r.x[i])
        with np.errstate(invalid="ignore", divide="ignore"):
            corr = cov / np.outer(sd, sd)
        info["std"] = {rename.get(nm, nm): v for nm, v in zip(names, sd_nat.tolist())}
        info["corr"] = np.nan_to_num(corr).round(3).tolist()
        info["corr_names"] = names  # order of the raw optimizer vector
    except np.linalg.LinAlgError:
        info["std"], info["corr"] = {}, []

    label = f"NLS {'R0-only' if n_rc == 0 else f'{n_rc}RC'}" + ("" if fit_soc else " (Q, SOC0 fixed)")
    return FitResult(label, p, info)


def bic(sse, n, k):
    """Bayesian information criterion (lower is better)."""
    return float(n * np.log(sse / n) + k * np.log(n))
