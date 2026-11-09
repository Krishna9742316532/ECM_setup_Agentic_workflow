"""Python versions of the methods in the MATLAB reference scripts
(Exercise3.m, Exercise4.m), checked with the same validation rules as the
main pipeline.

Run from the project root:
    PYTHONPATH=src python -m ecm.reference_methods

Output: results/reference_methods/ (summary.md, results.json, plots, CSV).
See notes/matlab_review.md for the review of the original MATLAB code.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from . import plots
from .checks import metrics
from .config import (MATLAB_ASSUMPTIONS, SYNTHETIC_2RC, SYNTHETIC_SOC_R0, SYNTHETIC_TRUTH,
                     Settings)
from .data import CellData, check_soc_direction, load_mat
from .export import export_timeseries_csv
from .ls import arx2_ls, arx_ls, difference_ls, separable_ls
from .model import coulomb_count, simulate, simulate_varying
from .nls import nls_fit
from .ocv import PolyOCV
from .pipeline import _r, select_model
from .rls import arx2_rls, arx_rls, r0_rls, r0_rls_sensitivity, rls
from .segmented import lookup, segmented_ls
from .synthetic import make_synthetic, make_synthetic_soc_r0, truth_params


def _rmse(V, Vm, m):
    return metrics(V, Vm, m)["rmse_mV"]


def _fitted(d, m_fit, s):
    p = nls_fit(d, m_fit, n_rc=0, ranges=s.ranges).params
    return p.Q_Ah, p.soc0


def _soc_setup(d, Q, soc0, ocv, sign=+1):
    soc = coulomb_count(d.t, sign * d.I, Q, soc0)
    return soc, d.V - ocv(soc)


def run(out_dir="results/reference_methods", data_path="data/ECM1.mat", s=Settings()):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    d = load_mat(data_path)
    m_fit, m_val = d.t < s.t_split, d.t >= s.t_split
    R = {}

    # Our SOC set-up: the model chosen by the main pipeline's selection.
    d_lin = d
    d, nls_sel, _, best_n, ocv_name = select_model(d_lin, m_fit, s)
    base = nls_sel[best_n]
    Q, soc0 = base.params.Q_Ah, base.params.soc0
    soc, y = _soc_setup(d, Q, soc0, d.ocv)
    soc_window = (float(soc.min()), float(soc.max()))

    # A. OCV model: table interpolation vs 10th-order polynomial ------------
    poly = PolyOCV(d.ocv.soc, d.ocv.ocv, MATLAB_ASSUMPTIONS["ocv_poly_order"])
    d_poly = CellData(d.t, d.I, d.V, poly, d.name + " (poly OCV)")
    fits = {}
    for name, dd in (("linear interpolation", d_lin), (f"polynomial order {poly.order}", d_poly)):
        fr = nls_fit(dd, m_fit, n_rc=0, ranges=s.ranges)
        Vm, _ = simulate(d.t, d.I, fr.params, dd.ocv)
        fits[name] = {"R0_mOhm": fr.params.R0 * 1e3, "Q_Ah": fr.params.Q_Ah,
                      "soc0": fr.params.soc0, "rmse_fit_mV": _rmse(d.V, Vm, m_fit),
                      "rmse_val_mV": _rmse(d.V, Vm, m_val)}
    R["A_ocv_model"] = {
        "poly_table_error_all_mV": {k: v * 1e3 for k, v in poly.table_error().items()},
        "poly_table_error_window_mV": {k: v * 1e3 for k, v in
                                       poly.table_error((soc_window[0] - 0.025,
                                                         soc_window[1] + 0.025)).items()},
        "poly_monotonic": poly.is_monotonic(), "nls_R0_only": fits,
        "selected_by_pipeline": ocv_name}
    plots.ocv_poly_vs_table(d_lin.ocv, poly, soc_window, out / "ocv_poly_vs_table.png")

    # B. MATLAB assumptions vs fitted SOC set-up -----------------------------
    mQ, ms0 = MATLAB_ASSUMPTIONS["Q_Ah"], MATLAB_ASSUMPTIONS["soc0"]
    setups = {
        "MATLAB Ex3 (poly10, Q=60 Ah, SOC0=0.5)": (poly, mQ, ms0, +1),
        "MATLAB Ex3 2RC part (same, flipped sign z-Ts/Q*I)": (poly, mQ, ms0, -1),
        "this project before (linear OCV, Q and SOC0 fitted)": (d_lin.ocv,) + _fitted(d_lin, m_fit, s) + (+1,),
        f"this project now ({ocv_name} OCV, Q={Q:.1f} Ah, SOC0={soc0:.3f}, fitted)": (d.ocv, Q, soc0, +1),
    }
    rows = []
    for name, (ocv, q, s0, sign) in setups.items():
        ss, yy = _soc_setup(d, q, s0, ocv, sign)
        r0 = separable_ls(d.t[m_fit], yy[m_fit], d.I[m_fit], 0).params.R0
        Vm = ocv(ss) + r0 * d.I
        rows.append({"setup": name, "soc_start": float(ss[0]), "soc_end": float(ss[-1]),
                     "soc_direction": check_soc_direction(d.I, ss) or ["OK"],
                     "LS_R0_mOhm": r0 * 1e3, "rmse_fit_mV": _rmse(d.V, Vm, m_fit),
                     "rmse_val_mV": _rmse(d.V, Vm, m_val)})
    R["B_soc_assumptions"] = rows

    # C. R0 from the simple methods of Exercise 3 ---------------------------
    diff = difference_ls(d.t[m_fit], d.V[m_fit], d.I[m_fit])
    rls_l1 = rls(d.I[m_fit][:, None], y[m_fit], lam=1.0, theta0=[0.01], P0=1.0)[-1, 0]
    R["C_R0_methods"] = [
        {"method": "Method 1: LS on V - OCV(SOC)",
         "R0_mOhm": separable_ls(d.t[m_fit], y[m_fit], d.I[m_fit], 0).params.R0 * 1e3,
         "needs_soc": True},
        {"method": "Method 2: difference LS dV = R0*dI (1 column)",
         "R0_mOhm": diff.params.R0 * 1e3, "needs_soc": False},
        {"method": "Method 2: MATLAB 2-column form, theta(1)",
         "R0_mOhm": diff.info["two_column_theta"][0] * 1e3, "needs_soc": False},
        {"method": "Method 2: MATLAB 2-column form, theta(2)",
         "R0_mOhm": diff.info["two_column_theta"][1] * 1e3, "needs_soc": False},
        {"method": "RLS (lambda = 1, P0 = 1) on V - OCV(SOC)", "R0_mOhm": rls_l1 * 1e3,
         "needs_soc": True},
        {"method": f"NLS R0-only, {ocv_name} OCV (reference, also fits Q and SOC0)",
         "R0_mOhm": base.params.R0 * 1e3, "needs_soc": False},
    ]

    # D. RLS forgetting factor and P0 study --------------------------------
    sens = r0_rls_sensitivity(d.t, y, d.I, t_report=s.t_split)
    R["D_rls_sensitivity"] = {"lambda": sens["lambda"], "P0": sens["P0"],
                              "lam_for_P0": sens["lam_for_P0"], "R0_ref_mOhm": sens["R0_ref_mOhm"]}
    plots.rls_lambda_study(d.t, sens, s.t_split, out / "rls_lambda_P0.png")

    # E. ARX family (1st and 2nd order) on real data + synthetic 2RC check --
    arx_rows = []
    for fr in (arx_ls(d.t[m_fit], y[m_fit], d.I[m_fit], 1),
               arx_ls(d.t[m_fit], y[m_fit], d.I[m_fit], 10),
               arx_rls(d.t, y, d.I, s.rls_lambda, 10, s.t_split),
               arx2_ls(d.t[m_fit], y[m_fit], d.I[m_fit], 1),
               arx2_ls(d.t[m_fit], y[m_fit], d.I[m_fit], 10),
               arx2_rls(d.t, y, d.I, s.rls_lambda, 10, s.t_split)):
        arx_rows.append({"data": "ECM1 (real)", "method": fr.method, **fr.params.to_dict(),
                         "physical": fr.valid and all(v > 0 for v in (*fr.params.R, *fr.params.tau))})
    syn_truth = truth_params(SYNTHETIC_2RC)
    for noise in (0.0, SYNTHETIC_2RC["noise_std"]):
        sd, p = make_synthetic(d, dict(SYNTHETIC_2RC, noise_std=noise), seed=s.seed)
        _, ys = _soc_setup(sd, p.Q_Ah, p.soc0, d.ocv)
        runs = [arx2_ls(d.t[m_fit], ys[m_fit], d.I[m_fit], 1),
                arx2_ls(d.t[m_fit], ys[m_fit], d.I[m_fit], 10)]
        if noise > 0:
            runs.append(nls_fit(sd, m_fit, n_rc=2, ranges=s.ranges))
        for fr in runs:
            arx_rows.append({"data": f"synthetic 2RC, noise {noise * 1e3:g} mV", "method": fr.method,
                             **fr.params.to_dict(),
                             "physical": fr.valid and all(v > 0 for v in (*fr.params.R, *fr.params.tau))})
    R["E_arx"] = {"rows": arx_rows, "synthetic_2rc_truth": syn_truth.to_dict()}

    # F. SOC-segmented LS --------------------------------------------------
    edges = np.round(np.arange(0.35, 0.5251, 0.025), 4)
    seg_real = segmented_ls(d.t, y, d.I, soc, edges, tau=None, mask=m_fit)
    R0_lut = lookup(seg_real, soc, "R0")
    V_lut = simulate_varying(d.t, d.I, soc, d.ocv, R0_lut)
    V_glob = d.ocv(soc) + base.params.R0 * d.I
    real_cmp = {"constant R0": {"fit": _rmse(d.V, V_glob, m_fit), "val": _rmse(d.V, V_glob, m_val)},
                "R0(SOC) table": {"fit": _rmse(d.V, V_lut, m_fit), "val": _rmse(d.V, V_lut, m_val)}}

    sd, r0_true = make_synthetic_soc_r0(d, SYNTHETIC_SOC_R0, seed=s.seed)
    T = SYNTHETIC_SOC_R0
    soc_s, ys = _soc_setup(sd, T["Q_Ah"], T["soc0"], d.ocv)
    tau_hat = separable_ls(d.t[m_fit], ys[m_fit], d.I[m_fit], 1).params.tau[0]
    seg_syn = segmented_ls(d.t, ys, d.I, soc_s, edges, tau=tau_hat, mask=m_fit)
    for r in seg_syn:
        if r["used"]:
            mid = (r["soc_lo"] + r["soc_hi"]) / 2
            r["R0_true_mid"] = float(r0_true(mid))
            r["R0_err_pct"] = 100 * (r["R0"] / r["R0_true_mid"] - 1)
    R0s, R1s = lookup(seg_syn, soc_s, "R0"), lookup(seg_syn, soc_s, "R1")
    V_lut_s = simulate_varying(d.t, d.I, soc_s, d.ocv, R0s, [R1s], [np.full_like(R1s, tau_hat)])
    sep = separable_ls(d.t[m_fit], ys[m_fit], d.I[m_fit], 1)
    V_glob_s, _ = simulate(d.t, d.I, sep.params.replace(Q_Ah=T["Q_Ah"], soc0=T["soc0"]), d.ocv)
    syn_cmp = {"constant R0, R1": {"fit": _rmse(sd.V, V_glob_s, m_fit), "val": _rmse(sd.V, V_glob_s, m_val)},
               "R0(SOC), R1(SOC) table": {"fit": _rmse(sd.V, V_lut_s, m_fit),
                                          "val": _rmse(sd.V, V_lut_s, m_val)}}
    R["F_segmented"] = {"edges": edges.tolist(), "real_rows": seg_real, "real_rmse_mV": real_cmp,
                        "synthetic_rows": seg_syn, "synthetic_tau_hat_s": tau_hat,
                        "synthetic_rmse_mV": syn_cmp}
    plots.soc_segments([("ECM1 (real): R0 per SOC bin", seg_real, None, base.params.R0),
                        ("Synthetic: R0 rises at low SOC", seg_syn, r0_true, sep.params.R0)],
                       out / "soc_segments.png")

    # G. Particle swarm optimisation (PSO_battery_parameter.m) -------------
    from .data import estimate_noise_std
    from .pso import matlab_style_pso, pso_fit
    syn1, p1 = make_synthetic(d_lin, SYNTHETIC_TRUTH, seed=s.seed)
    pso_rows, curves = [], []

    def _row(data_name, fr, truth=None):
        p = fr.params
        row = {"data": data_name, "method": fr.method, "seed": fr.info.get("seed"),
               **p.to_dict(), "at_bound": fr.info.get("at_bound", []),
               "rmse_mV": fr.info.get("rmse_polished_mV", fr.info.get("rmse_pso_mV",
                                                                        fr.info.get("rmse_mV")))}
        if truth is not None and p.n_rc:
            row["err_pct"] = {"R0": 100 * (p.R0 / truth.R0 - 1), "R1": 100 * (p.R[0] / truth.R[0] - 1),
                              "tau1": 100 * (p.tau[0] / truth.tau[0] - 1)}
        return row

    for seed in (1, 2, 3):
        fr = pso_fit(syn1, m_fit, n_rc=1, ranges=s.ranges, iters=60, seed=seed, polish=False)
        pso_rows.append(_row("synthetic 1RC", fr, p1))
        curves.append(("this project, synthetic (seeds 1-3)" if seed == 1 else "_nolegend_",
                       fr.info["history_rmse_mV"], 0))
    pso_rows.append(_row("synthetic 1RC", pso_fit(syn1, m_fit, n_rc=1, ranges=s.ranges, iters=60,
                                                    seed=1, polish=True), p1))
    for n_rc in (0, 1):
        for seed in (1, 2, 3):
            fr = pso_fit(d, m_fit, n_rc=n_rc, ranges=s.ranges, iters=60, seed=seed, polish=True)
            pso_rows.append(_row(f"ECM1 real ({ocv_name} OCV)", fr))
            if n_rc == 1 and seed == 1:
                curves.append(("this project, real 1RC, seed 1", fr.info["history_rmse_mV"], 1))
    head = lambda dd: CellData(dd.t[m_fit], dd.I[m_fit], dd.V[m_fit], dd.ocv, dd.name)  # noqa: E731
    mr = matlab_style_pso(head(d_poly), Q_Ah=mQ, soc0=ms0)
    pso_rows.append(_row("ECM1 real (MATLAB set-up: poly10, Q=60, SOC0=0.5)", mr))
    curves.append(("MATLAB set-up, real", mr.info["history_rmse_mV"], 2))
    ms = matlab_style_pso(head(syn1), Q_Ah=p1.Q_Ah, soc0=p1.soc0)
    pso_rows.append(_row("synthetic 1RC (MATLAB set-up, true Q and SOC0)", ms, p1))
    curves.append(("MATLAB set-up, synthetic", ms.info["history_rmse_mV"], 3))
    R["G_pso"] = {"rows": pso_rows, "synthetic_truth": p1.to_dict(),
                  "matlab_bounds": {"lb [R0, R1, C1]": mr.info["lb"], "ub [R0, R1, C1]": mr.info["ub"]}}
    plots.pso_convergence(curves, estimate_noise_std(d.V, d.I), out / "pso_convergence.png")

    # H. CSV export --------------------------------------------------------
    online = r0_rls(d.t, y, d.I, lam=s.rls_lambda)
    export_timeseries_csv(out / "timeseries_ECM1.csv", {
        "time_s": d.t, "current_A": d.I, "voltage_V": d.V, "soc": soc, "ocv_V": d.ocv(soc),
        "v_model_V": V_glob, "residual_mV": (d.V - V_glob) * 1e3,
        "R0_rls_online_mOhm": online.info["R0"] * 1e3, "R0_soc_table_mOhm": R0_lut * 1e3,
    }, notes=[
        "ECM1.mat (public). Negative current = discharge.",
        f"soc: Coulomb counting with Q = {Q:.2f} Ah, SOC0 = {soc0:.4f} (fitted, {base.method}, t < {s.t_split:g} s)",
        f"ocv_V: OCV model '{ocv_name}' of the OCV table. v_model_V: OCV + R0*I with the global R0.",
        f"R0_rls_online_mOhm: RLS R0 + offset, lambda = {s.rls_lambda}. Estimates, NOT ground truth.",
        "R0_soc_table_mOhm: SOC-segmented LS (2.5 % bins, fit window only). Estimates, NOT ground truth.",
    ])
    R["H_export"] = {"file": "timeseries_ECM1.csv", "rows": len(d.t)}
    R["soc_setup"] = {"Q_Ah": Q, "soc0": soc0, "soc_range": soc_window}

    R = _r(R)
    (out / "results.json").write_text(json.dumps(R, indent=2))
    (out / "summary.md").write_text(render(R))
    return R


def _f(x, f="{:.3f}"):
    return "-" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f.format(x)


def render(R):
    L = ["# MATLAB reference methods - Python results (ECM1)", "",
         "Methods from `Exercise3.m` / `Exercise4.m`, re-implemented and checked. "
         "Review of the MATLAB code: `notes/matlab_review.md`.", ""]
    A = R["A_ocv_model"]
    L += ["## A. OCV model: interpolation vs polynomial (MATLAB polyfit order 10)", "",
          f"- Polynomial error at the table points: max {A['poly_table_error_all_mV']['max_abs_V']:.1f} mV, "
          f"RMS {A['poly_table_error_all_mV']['rms_V']:.1f} mV",
          f"- In this test's SOC window: max {A['poly_table_error_window_mV']['max_abs_V']:.1f} mV "
          "(noise is about 10 mV)", f"- Monotonic: {A['poly_monotonic']}", "",
          "| OCV model | R0 [mOhm] | Q [Ah] | SOC0 | fit RMSE [mV] | val RMSE [mV] |", "|---|---|---|---|---|---|"]
    for k, v in A["nls_R0_only"].items():
        L.append(f"| {k} | {v['R0_mOhm']:.3f} | {v['Q_Ah']:.2f} | {v['soc0']:.4f} | "
                 f"{v['rmse_fit_mV']:.2f} | {v['rmse_val_mV']:.2f} |")
    L += ["", f"OCV model chosen by the main pipeline: **{A['selected_by_pipeline']}**.", "",
          "## B. SOC assumptions: MATLAB scripts vs this project", "",
          "| Set-up | SOC start -> end | SOC direction | LS R0 [mOhm] | fit RMSE [mV] | val RMSE [mV] |",
          "|---|---|---|---|---|---|"]
    for r in R["B_soc_assumptions"]:
        L.append(f"| {r['setup']} | {r['soc_start']:.3f} -> {r['soc_end']:.3f} | "
                 f"{'; '.join(r['soc_direction'])} | {r['LS_R0_mOhm']:.3f} | "
                 f"{r['rmse_fit_mV']:.2f} | {r['rmse_val_mV']:.2f} |")
    L += ["", "## C. R0 from the simple methods (fit window)", "",
          "| Method | R0 [mOhm] | needs Q / SOC0 / OCV |", "|---|---|---|"]
    for r in R["C_R0_methods"]:
        L.append(f"| {r['method']} | {r['R0_mOhm']:.4f} | {'yes' if r['needs_soc'] else 'no'} |")
    D = R["D_rls_sensitivity"]
    L += ["", "## D. RLS for R0: forgetting factor and initial covariance", "",
          "| lambda | memory [s] | R0 at split [mOhm] | std of R0(t) [mOhm] |", "|---|---|---|---|"]
    for r in D["lambda"]:
        mem = "inf" if r["memory_s"] is None else f"{r['memory_s']:.1f}"
        L.append(f"| {r['lambda']:g} | {mem} | {r['R0_mOhm']:.4f} | {r['R0_std_mOhm']:.4f} |")
    L += ["", f"P0 study (lambda = {D['lam_for_P0']:g}, start value 10 mOhm, reference R0 = "
              f"{D['R0_ref_mOhm']:.4f} mOhm):", "",
          "| P0 | R0 at split [mOhm] | time to first reach +-5 % of reference [s] |", "|---|---|---|"]
    for r in D["P0"]:
        L.append(f"| {r['P0']:g} | {r['R0_mOhm']:.4f} | {_f(r['t_reach_5pct_s'], '{:.1f}')} |")
    E = R["E_arx"]
    L += ["", "## E. ARX methods (1RC and 2RC)", "",
          "Synthetic 2RC truth: " + ", ".join(f"{k} = {v:g}" for k, v in E["synthetic_2rc_truth"].items()), "",
          "| Data | Method | R0 [mOhm] | R1 [mOhm] | tau1 [s] | R2 [mOhm] | tau2 [s] | physical |",
          "|---|---|---|---|---|---|---|---|"]
    for r in E["rows"]:
        L.append(f"| {r['data']} | {r['method']} | {_f(r.get('R0_mOhm'))} | {_f(r.get('R1_mOhm'))} | "
                 f"{_f(r.get('tau1_s'), '{:.2f}')} | {_f(r.get('R2_mOhm'))} | {_f(r.get('tau2_s'), '{:.2f}')} | "
                 f"{'yes' if r['physical'] else 'NO'} |")
    F = R["F_segmented"]
    L += ["", "## F. SOC-segmented LS", "", "Real data (R0-only, 2.5 % SOC bins, fit window):", "",
          "| SOC bin | samples | R0 [mOhm] | RMSE [mV] |", "|---|---|---|---|"]
    for r in F["real_rows"]:
        last = _f(r.get("rmse_mV"), "{:.2f}") if r["used"] else f"skipped: {r.get('reason')}"
        L.append(f"| {r['soc_lo']:.3f}-{r['soc_hi']:.3f} | {r['n']} | "
                 f"{_f(r.get('R0') and r['R0'] * 1e3, '{:.4f}')} | {last} |")
    L += ["", "| Model | fit RMSE [mV] | val RMSE [mV] |", "|---|---|---|"]
    for k, v in F["real_rmse_mV"].items():
        L.append(f"| {k} | {v['fit']:.2f} | {v['val']:.2f} |")
    L += ["", f"Synthetic data with R0(SOC) (tau fixed at {F['synthetic_tau_hat_s']:.1f} s from separable LS):", "",
          "| SOC bin | R0 [mOhm] | true R0 at bin centre | error [%] | R1 [mOhm] |", "|---|---|---|---|---|"]
    for r in F["synthetic_rows"]:
        if r["used"]:
            L.append(f"| {r['soc_lo']:.3f}-{r['soc_hi']:.3f} | {r['R0'] * 1e3:.4f} | "
                     f"{r['R0_true_mid'] * 1e3:.4f} | {r['R0_err_pct']:+.2f} | {r['R1'] * 1e3:.3f} |")
    L += ["", "| Model | fit RMSE [mV] | val RMSE [mV] |", "|---|---|---|"]
    for k, v in F["synthetic_rmse_mV"].items():
        L.append(f"| {k} | {v['fit']:.2f} | {v['val']:.2f} |")
    H = R["G_pso"]
    L += ["", "## G. Particle swarm optimisation (PSO)", "",
          "Synthetic truth: " + ", ".join(f"{k} = {v:g}" for k, v in H["synthetic_truth"].items()), "",
          f"MATLAB bounds: {H['matlab_bounds']}", "",
          "| Data | Method | seed | R0 [mOhm] | R1 [mOhm] | tau1 [s] | C1 [F] | Q [Ah] | RMSE [mV] | on bound | "
          "error vs truth (R0 / R1 / tau1) |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in H["rows"]:
        e = r.get("err_pct")
        es = "-" if not e else f"{e['R0']:+.1f} % / {e['R1']:+.1f} % / {e['tau1']:+.1f} %"
        L.append(f"| {r['data']} | {r['method']} | {r['seed'] if r['seed'] is not None else '-'} | "
                 f"{_f(r.get('R0_mOhm'))} | {_f(r.get('R1_mOhm'))} | {_f(r.get('tau1_s'), '{:.2f}')} | "
                 f"{_f(r.get('C1_F'), '{:.0f}')} | {_f(r.get('Q_Ah'), '{:.2f}')} | {_f(r.get('rmse_mV'), '{:.2f}')} | "
                 f"{', '.join(r['at_bound']) or '-'} | {es} |")
    L += ["", "## H. CSV export", "",
          f"`{R['H_export']['file']}` ({R['H_export']['rows']} rows). The header explains every column. "
          "Identified parameters are estimates, not ground truth.", ""]
    return "\n".join(L)


def main(argv=None):
    ap = argparse.ArgumentParser(description="MATLAB reference methods in Python")
    ap.add_argument("--data", default="data/ECM1.mat")
    ap.add_argument("--out", default="results/reference_methods")
    a = ap.parse_args(argv)
    run(a.out, a.data)
    print(f"Results written to {a.out}/")


if __name__ == "__main__":
    main()
