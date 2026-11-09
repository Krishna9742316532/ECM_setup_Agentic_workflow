"""End-to-end identification and validation pipeline.

Run from the project root:
    python -m ecm.pipeline                      # real data + synthetic check
    python -m ecm.pipeline --only real

Outputs go to results/<dataset>/: params.json, summary.md and PNG plots.
"""

import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import scipy

from . import plots
from .checks import (acf, fit_quality, metrics, plausibility, rc_identifiability,
                     slow_drift_mV, worst_level)
from .config import SYNTHETIC_TRUTH, Settings
from .data import (CellData, check_data, current_step_resistance, estimate_noise_std,
                   load_mat)
from .ls import arx_ls, separable_ls
from .model import coulomb_count, simulate
from .nls import bic, nls_fit
from .ocv import OCVCurve, PolyOCV
from .rls import arx_rls, fixed_tau_rls, r0_rls
from .synthetic import make_synthetic


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _r(x, n=6):
    """Round floats for stable, readable JSON."""
    if isinstance(x, dict):
        return {k: _r(v, n) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_r(v, n) for v in x]
    if isinstance(x, (float, np.floating)):
        x = float(x)
        return None if not np.isfinite(x) else float(f"{x:.{n}g}")
    if isinstance(x, np.integer):
        return int(x)
    return x


def evaluate(d, fr, m_fit, m_val, noise, s, check_soc):
    """Simulate a fitted parameter set on the full profile and run all checks."""
    out = {"method": fr.method, "params": fr.params.to_dict(), "valid": fr.valid}
    if not fr.valid:
        out["findings"] = [str(f) for f in plausibility(fr.params, s.ranges, check_soc)]
        out["level"] = "FAIL"
        out["fit"] = out["val"] = {"rmse_mV": float("nan"), "max_abs_mV": float("nan"),
                                   "mean_mV": float("nan")}
        return out, None
    Vm, _ = simulate(d.t, d.I, fr.params, d.ocv)
    res = d.V - Vm
    out["fit"], out["val"] = metrics(d.V, Vm, m_fit), metrics(d.V, Vm, m_val)
    drift = slow_drift_mV(res[m_fit])
    out["residual"] = {"acf": acf(res[m_fit]), "slow_drift_fit_mV": drift}
    findings = (plausibility(fr.params, s.ranges, check_soc)
                + rc_identifiability(fr.params, d.t[m_fit], d.I[m_fit], noise)
                + fit_quality(out["fit"], out["val"], noise, drift))
    out["findings"] = [str(f) for f in findings]
    out["level"] = worst_level(findings)
    return out, Vm


def make_ocv(d, name):
    """'linear' = table interpolation; 'polyN' = polynomial of order N."""
    if name == "linear":
        return d.ocv if not isinstance(d.ocv, PolyOCV) else OCVCurve(d.ocv.soc, d.ocv.ocv)
    if name.startswith("poly"):
        return PolyOCV(d.ocv.soc, d.ocv.ocv, int(name[4:]))
    raise ValueError(f"unknown OCV model {name!r}")


def select_model(d, m_fit, s=Settings()):
    """Fit NLS with 0/1/2 RC branches for every OCV model in s.ocv_models.

    Returns (data with the chosen OCV, NLS fits for that OCV, selection table,
    chosen number of RC branches, chosen OCV model name).
    """
    sel, fits, datas = [], {}, {}
    for name in s.ocv_models:
        dd = CellData(d.t, d.I, d.V, make_ocv(d, name), d.name)
        datas[name] = dd
        fits[name] = {n: nls_fit(dd, m_fit, n_rc=n, ranges=s.ranges) for n in (0, 1, 2)}
        for n, fr in fits[name].items():
            # Not supported: a parameter on/near a bound, or outside its range.
            bad = bool(fr.info["at_bound"]) or worst_level(plausibility(fr.params, s.ranges)) != "OK"
            sel.append({"ocv": name, "model": fr.method, "n_rc": n, "k": fr.info["k"],
                        "rmse_fit_mV": fr.info["rmse_mV"],
                        "bic": bic(fr.info["sse"], fr.info["n"], fr.info["k"]),
                        "at_bound": fr.info["at_bound"], "supported": not bad,
                        "Q_Ah": fr.params.Q_Ah, "soc0": fr.params.soc0})
    ok = [r for r in sel if r["supported"]] or sel
    b = min(ok, key=lambda r: r["bic"])
    return datas[b["ocv"]], fits[b["ocv"]], sel, b["n_rc"], b["ocv"]


def run(d, out_dir, s=Settings(), truth=None, data_sha=None, make_plots=True):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    m_fit, m_val = d.t < s.t_split, d.t >= s.t_split

    # 1. Data checks -------------------------------------------------------
    issues = check_data(d)
    if issues:
        raise ValueError("data checks failed: " + "; ".join(issues))
    noise = estimate_noise_std(d.V, d.I)
    summary = {
        "dataset": d.name, "data_sha256": data_sha,
        "data": {"n": len(d.t), "dt_s": d.dt, "duration_s": float(d.t[-1] - d.t[0]),
                 "I_min_A": float(d.I.min()), "I_max_A": float(d.I.max()),
                 "V_min_V": float(d.V.min()), "V_max_V": float(d.V.max()),
                 "net_charge_Ah": float(np.sum(d.I[:-1] * np.diff(d.t)) / 3600),
                 "noise_std_mV": noise * 1e3,
                 "R0_from_current_steps_mOhm": current_step_resistance(d) * 1e3,
                 "checks": "all passed"},
        "settings": {"t_split_s": s.t_split, "rls_lambda": s.rls_lambda,
                     "ocv_models": list(s.ocv_models),
                     "arx_block": s.arx_block, "ranges": s.ranges},
        "versions": {"python": platform.python_version(), "numpy": np.__version__,
                     "scipy": scipy.__version__},
    }

    # 2. NLS output-error fits; joint choice of OCV model and RC branches ---
    d, nls, sel, best_n, ocv_name = select_model(d, m_fit, s)
    summary["model_selection"] = {
        "table": sel, "recommended": nls[best_n].method, "ocv_model": ocv_name,
        "rule": "lowest BIC among (OCV model, RC branches) whose parameters are all "
                "inside their ranges and not on a bound"}
    best = nls[best_n]

    # 3. SOC set-up for LS / RLS (Q and SOC0 from the recommended NLS model)
    Q, soc0 = best.params.Q_Ah, best.params.soc0
    soc = coulomb_count(d.t, d.I, Q, soc0)
    y = d.V - d.ocv(soc)
    tf, yf, If = d.t[m_fit], y[m_fit], d.I[m_fit]
    summary["soc_setup"] = {"source": f"{best.method}, OCV {ocv_name}", "Q_Ah": Q, "soc0": soc0,
                            "soc_range": [float(soc.min()), float(soc.max())]}

    # 4. LS and RLS --------------------------------------------------------
    ls_runs = [
        separable_ls(tf, yf, If, 0, Q_Ah=Q, soc0=soc0),
        arx_ls(tf, yf, If, block=1, Q_Ah=Q, soc0=soc0),
        arx_ls(tf, yf, If, block=s.arx_block, Q_Ah=Q, soc0=soc0),
        separable_ls(tf, yf, If, 1, s.ranges["tau"], s.tau_grid_n, Q_Ah=Q, soc0=soc0),
        separable_ls(tf, yf, If, 2, s.ranges["tau"], 30, Q_Ah=Q, soc0=soc0),
    ]
    sep1 = ls_runs[3]
    rls_r0 = r0_rls(d.t, y, d.I, lam=s.rls_lambda, t_report=s.t_split, Q_Ah=Q, soc0=soc0)
    rls_arx =arx_rls(d.t, y, d.I, lam=s.rls_lambda, block=s.arx_block,
                      t_report=s.t_split, Q_Ah=Q, soc0=soc0)
    rls_fix = fixed_tau_rls(d.t, y, d.I, sep1.params.tau[0], lam=s.rls_lambda,
                            t_report=s.t_split, Q_Ah=Q, soc0=soc0)

    # 5. Evaluate every method on fit and validation data ------------------
    methods, sims = [], {}
    for fr, check_soc in ([(nls[n], True) for n in (0, 1, 2)]
                          + [(r, False) for r in ls_runs]
                          + [(rls_r0, False), (rls_arx, False), (rls_fix, False)]):
        ev, Vm = evaluate(d, fr, m_fit, m_val, noise, s, check_soc)
        if fr.method.startswith("NLS"):
            ev["std_errors"] = fr.info.get("std", {})
            ev["at_bound"] = fr.info["at_bound"]
        if "offset_at_report_mV" in fr.info:
            ev["ocv_offset_at_report_mV"] = fr.info["offset_at_report_mV"]
        if truth is not None and fr.valid:
            ev["error_vs_truth_pct"] = _truth_error(fr.params, truth, check_soc)
        methods.append(ev)
        sims[fr.method] = Vm
    summary["methods"] = methods

    # 6. Capacity / SOC0 sensitivity (information only, not used) ----------
    sens = []
    for lo, hi in ((0, s.t_split), (0, d.t[-1] + 1), (600, d.t[-1] + 1)):
        m = (d.t >= lo) & (d.t < hi)
        fr = nls_fit(d, m, n_rc=best_n, ranges=s.ranges)
        sens.append({"window_s": [lo, float(min(hi, d.t[-1]))], "Q_Ah": fr.params.Q_Ah,
                     "soc0": fr.params.soc0, "R0_mOhm": fr.params.R0 * 1e3,
                     "rmse_mV": fr.info["rmse_mV"], "Q_std_Ah": fr.info["std"].get("Q_Ah")})
    summary["capacity_sensitivity"] = sens

    if truth is not None:
        summary["truth"] = truth.to_dict()

    # 7. Save --------------------------------------------------------------
    summary = _r(summary)
    (out_dir / "params.json").write_text(json.dumps(summary, indent=2))
    (out_dir / "summary.md").write_text(render_markdown(summary))
    if make_plots:
        plots.data_overview(d, soc, s.t_split, out_dir / "data_overview.png")
        plots.ocv_window(d.ocv, soc, out_dir / "ocv_window.png")
        plots.voltage_fit(d, sims[best.method], s.t_split, best.method, out_dir / "voltage_fit.png")
        plots.residual_stats((d.V - sims[best.method])[m_fit], noise, out_dir / "residuals.png")
        plots.method_comparison([(m["method"], m["fit"]["rmse_mV"], m["val"]["rmse_mV"])
                                 for m in methods if m["valid"]], noise,
                                out_dir / "method_comparison.png")
        plots.rls_tracking(rls_arx, rls_fix, s.t_split, out_dir / "rls_tracking.png", truth)
        plots.rls_r0_offset(rls_r0, s.t_split, out_dir / "rls_r0_offset.png",
                            truth.R0 if truth is not None else None)
        plots.tau_scan(sep1, int(m_fit.sum()), out_dir / "tau_scan.png",
                       truth.tau[0] if truth is not None and truth.n_rc else None)
    return summary


def _truth_error(p, truth, check_soc):
    """Percentage error vs the known truth (synthetic data only)."""
    err = {"R0": 100 * (p.R0 / truth.R0 - 1)}
    # Compare the RC branch closest to the true one (extra branches ignored).
    if truth.n_rc and p.n_rc:
        i = int(np.argmax(p.R))
        err["R1"] = 100 * (p.R[i] / truth.R[0] - 1)
        err["tau1"] = 100 * (p.tau[i] / truth.tau[0] - 1)
    if check_soc:
        err["Q_Ah"] = 100 * (p.Q_Ah / truth.Q_Ah - 1)
        err["soc0_abs_pct_points"] = 100 * (p.soc0 - truth.soc0)
    return err


def _fmt(x, f="{:.3f}"):
    return "-" if x is None else f.format(x)


def render_markdown(s):
    L = [f"# Results: {s['dataset']}", ""]
    d = s["data"]
    L += ["## Data", "",
          f"- {d['n']} samples, dt = {d['dt_s']} s, duration {d['duration_s']} s",
          f"- Current {d['I_min_A']} to {d['I_max_A']} A (negative = discharge), "
          f"net charge {d['net_charge_Ah']:.2f} Ah",
          f"- Voltage {d['V_min_V']:.3f} to {d['V_max_V']:.3f} V",
          f"- Estimated voltage noise: **{d['noise_std_mV']:.2f} mV** (model-free estimate)",
          f"- R0 from current steps: {d['R0_from_current_steps_mOhm']:.3f} mOhm",
          f"- Data checks: {d['checks']}", ""]
    if s.get("data_sha256"):
        L += [f"- SHA-256: `{s['data_sha256']}`", ""]

    L += ["## Model selection (NLS, fitted on t < split)", "",
          "| OCV model | Model | params | Q [Ah] | SOC0 | fit RMSE [mV] | BIC | params on bound | supported |",
          "|---|---|---|---|---|---|---|---|---|"]
    for r in s["model_selection"]["table"]:
        L.append(f"| {r['ocv']} | {r['model']} | {r['k']} | {r['Q_Ah']:.2f} | {r['soc0']:.4f} | "
                 f"{r['rmse_fit_mV']:.3f} | {r['bic']:.1f} | "
                 f"{', '.join(r['at_bound']) or '-'} | {'yes' if r['supported'] else 'no'} |")
    L += ["", f"Recommended: **{s['model_selection']['recommended']}** with OCV model "
              f"**{s['model_selection']['ocv_model']}** ({s['model_selection']['rule']}).", ""]

    so = s["soc_setup"]
    L += ["## SOC set-up used by LS and RLS", "",
          f"Q = {so['Q_Ah']:.2f} Ah, SOC0 = {so['soc0']:.4f} (from {so['source']}); "
          f"SOC range {so['soc_range'][0]:.3f} to {so['soc_range'][1]:.3f}.", ""]

    L += ["## All methods", "",
          "| Method | R0 [mOhm] | R1 [mOhm] | tau1 [s] | R2 [mOhm] | tau2 [s] | Q [Ah] | SOC0 | "
          "fit RMSE [mV] | val RMSE [mV] | check |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    for m in s["methods"]:
        p = m["params"]
        L.append(f"| {m['method']} | {_fmt(p.get('R0_mOhm'))} | {_fmt(p.get('R1_mOhm'))} | "
                 f"{_fmt(p.get('tau1_s'), '{:.2f}')} | {_fmt(p.get('R2_mOhm'))} | "
                 f"{_fmt(p.get('tau2_s'), '{:.2f}')} | {_fmt(p.get('Q_Ah'), '{:.2f}')} | "
                 f"{_fmt(p.get('soc0'), '{:.4f}')} | {_fmt(m['fit']['rmse_mV'], '{:.2f}')} | "
                 f"{_fmt(m['val']['rmse_mV'], '{:.2f}')} | {m['level']} |")
    L += ["", "LS and RLS methods use the Q and SOC0 from the SOC set-up above.", ""]

    if "truth" in s:
        t = s["truth"]
        L += ["## Known truth (synthetic data)", "",
              ", ".join(f"{k} = {v:g}" for k, v in t.items()), "",
              "| Method | R0 err [%] | R1 err [%] | tau1 err [%] | Q err [%] |",
              "|---|---|---|---|---|"]
        for m in s["methods"]:
            e = m.get("error_vs_truth_pct")
            if e:
                L.append(f"| {m['method']} | {_fmt(e.get('R0'), '{:+.2f}')} | "
                         f"{_fmt(e.get('R1'), '{:+.1f}')} | {_fmt(e.get('tau1'), '{:+.1f}')} | "
                         f"{_fmt(e.get('Q_Ah'), '{:+.2f}')} |")
        L.append("")

    L += ["## Findings per method", ""]
    for m in s["methods"]:
        L.append(f"**{m['method']}** ({m['level']})")
        L += [f"- {f}" for f in m["findings"]]
        L.append("")

    L += ["## Capacity sensitivity (information only)", "",
          "Same recommended model fitted on different time windows.", "",
          "| Window [s] | Q [Ah] | Q std [Ah] | SOC0 | R0 [mOhm] | RMSE [mV] |",
          "|---|---|---|---|---|---|"]
    for r in s["capacity_sensitivity"]:
        L.append(f"| {r['window_s'][0]:g}-{r['window_s'][1]:g} | {r['Q_Ah']:.2f} | "
                 f"{_fmt(r['Q_std_Ah'], '{:.2f}')} | {r['soc0']:.4f} | {r['R0_mOhm']:.3f} | "
                 f"{r['rmse_mV']:.2f} |")
    L += ["", f"Versions: {s['versions']}", ""]
    return "\n".join(L)


def main(argv=None):
    ap = argparse.ArgumentParser(description="ECM parameter identification")
    ap.add_argument("--data", default="data/ECM1.mat")
    ap.add_argument("--out", default="results")
    ap.add_argument("--only", choices=["real", "synthetic"], default=None)
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--ocv", nargs="+", default=None,
                    help="OCV models to compare, e.g. linear poly10 (default: from config)")
    a = ap.parse_args(argv)
    s = Settings() if a.ocv is None else Settings(ocv_models=tuple(a.ocv))
    real = load_mat(a.data)
    sha = sha256(a.data)
    if a.only in (None, "synthetic"):
        syn, truth = make_synthetic(real, SYNTHETIC_TRUTH, seed=s.seed)
        r = run(syn, Path(a.out) / "synthetic", s, truth=truth, make_plots=not a.no_plots)
        print(f"[synthetic] recommended: {r['model_selection']['recommended']} "
              f"(OCV: {r['model_selection']['ocv_model']})")
    if a.only in (None, "real"):
        r = run(real, Path(a.out) / "real", s, data_sha=sha, make_plots=not a.no_plots)
        print(f"[real]      recommended: {r['model_selection']['recommended']} "
              f"(OCV: {r['model_selection']['ocv_model']})")
    print(f"Results written to {a.out}/")


if __name__ == "__main__":
    main()
