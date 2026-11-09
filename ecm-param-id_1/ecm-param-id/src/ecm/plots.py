"""Plots saved as PNG files. One measure per axis; split time marked."""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# Categorical colours in fixed order; measured data in neutral grey.
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
GREY, INK, MUTED = "#8a8a85", "#2b2b29", "#6b6b66"

plt.rcParams.update({
    "figure.dpi": 110, "savefig.dpi": 140, "font.size": 9,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "axes.titlesize": 10,
    "axes.titleweight": "bold", "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": "#e4e3de", "grid.linewidth": 0.6,
    "xtick.color": MUTED, "ytick.color": MUTED, "legend.frameon": False,
    "lines.linewidth": 1.2,
})


def _split(ax, t_split):
    ax.axvline(t_split, color=MUTED, ls="--", lw=0.8)


def _save(fig, path):
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def data_overview(d, soc, t_split, path):
    fig, ax = plt.subplots(3, 1, figsize=(9, 6.5), sharex=True)
    ax[0].plot(d.t, d.I, color=BLUE, lw=0.5)
    ax[0].set_ylabel("Current [A]\n(neg. = discharge)")
    ax[1].plot(d.t, d.V, color=BLUE, lw=0.5)
    ax[1].set_ylabel("Voltage [V]")
    ax[2].plot(d.t, soc * 100, color=BLUE)
    ax[2].set_ylabel("SOC [%]")
    ax[2].set_xlabel("Time [s]")
    for a in ax:
        _split(a, t_split)
    ax[0].set_title(f"{d.name}: input data (dashed line = fit | validation split)")
    _save(fig, path)


def ocv_window(ocv, soc, path):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.axvspan(soc.min() * 100, soc.max() * 100, color=BLUE, alpha=0.12, lw=0,
               label="SOC range used in this test")
    ax.plot(ocv.soc * 100, ocv.ocv, "-o", color=BLUE, ms=3, label="OCV table")
    ax.set_xlabel("SOC [%]")
    ax.set_ylabel("OCV [V]")
    ax.set_title("OCV curve and operating window")
    ax.legend(loc="lower right")
    _save(fig, path)


def voltage_fit(d, Vm, t_split, label, path, zoom=(800, 1000)):
    res = (d.V - Vm) * 1e3
    fig, ax = plt.subplots(3, 1, figsize=(9, 7.5))
    ax[0].plot(d.t, d.V, color=GREY, lw=0.5, label="measured")
    ax[0].plot(d.t, Vm, color=BLUE, lw=0.8, label=f"model: {label}")
    ax[0].set_ylabel("Voltage [V]")
    ax[0].legend(loc="lower left", ncol=2)
    ax[0].set_title("Measured vs model voltage")
    _split(ax[0], t_split)
    m = (d.t >= zoom[0]) & (d.t <= zoom[1])
    ax[1].plot(d.t[m], d.V[m], color=GREY, lw=0.6, label="measured")
    ax[1].plot(d.t[m], Vm[m], color=BLUE, lw=1.0, label="model")
    ax[1].set_ylabel("Voltage [V]")
    ax[1].set_title(f"Zoom {zoom[0]:g}-{zoom[1]:g} s")
    ax[2].plot(d.t, res, color=GREY, lw=0.4, label="residual")
    k = np.ones(100) / 100
    ax[2].plot(d.t, np.convolve(res, k, mode="same"), color=ORANGE, lw=1.4,
               label="10 s moving average")
    ax[2].axhline(0, color=INK, lw=0.6)
    ax[2].set_ylabel("Residual V_meas - V_model [mV]")
    ax[2].set_xlabel("Time [s]")
    ax[2].legend(loc="upper left", ncol=2)
    _split(ax[2], t_split)
    _save(fig, path)


def residual_stats(res_V, noise_std, path, max_lag=300):
    res = res_V * 1e3
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    ax[0].hist(res, bins=80, density=True, color=BLUE, alpha=0.8, label="residuals")
    x = np.linspace(res.min(), res.max(), 200)
    s = noise_std * 1e3
    ax[0].plot(x, np.exp(-0.5 * (x / s) ** 2) / (s * np.sqrt(2 * np.pi)), color=ORANGE,
               lw=1.6, label=f"Gaussian, std = noise ({s:.1f} mV)")
    ax[0].set_xlabel("Residual [mV]")
    ax[0].set_ylabel("Density")
    ax[0].set_title("Residual distribution")
    ax[0].legend(loc="upper left", fontsize=8)
    r = res - res.mean()
    lags = np.arange(1, max_lag + 1)
    c = np.array([r[:-k] @ r[k:] for k in lags]) / (r @ r)
    ax[1].plot(lags * 0.1, c, color=BLUE)
    band = 2 / np.sqrt(len(r))
    ax[1].axhspan(-band, band, color=GREY, alpha=0.25, lw=0, label="white-noise band (95%)")
    ax[1].set_xlabel("Lag [s]")
    ax[1].set_ylabel("Autocorrelation")
    ax[1].set_title("Residual autocorrelation")
    ax[1].legend(loc="upper right", fontsize=8)
    _save(fig, path)


def method_comparison(rows, noise_std, path):
    """rows: list of (method, rmse_fit, rmse_val); invalid methods skipped."""
    rows = [r for r in rows if np.isfinite(r[1]) and np.isfinite(r[2])]
    names = [r[0] for r in rows]
    y = np.arange(len(rows))
    h = 0.38
    fig, ax = plt.subplots(figsize=(8, 0.45 * len(rows) + 1.6))
    ax.barh(y - h / 2, [r[1] for r in rows], height=h - 0.04, color=BLUE, label="fit (t < split)")
    ax.barh(y + h / 2, [r[2] for r in rows], height=h - 0.04, color=ORANGE, label="validation (t >= split)")
    ax.axvline(noise_std * 1e3, color=INK, ls=":", lw=1.2, label=f"noise level {noise_std * 1e3:.1f} mV")
    for yi, r in zip(y, rows):
        ax.text(r[1] + 0.2, yi - h / 2, f"{r[1]:.2f}", va="center", fontsize=7, color=INK)
        ax.text(r[2] + 0.2, yi + h / 2, f"{r[2]:.2f}", va="center", fontsize=7, color=INK)
    ax.set_yticks(y, names)
    ax.invert_yaxis()
    ax.set_xlabel("Voltage RMSE [mV]")
    ax.set_title("Voltage error per method")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.09 - 1.0 / (len(rows) + 4)),
              ncol=3, fontsize=8)
    ax.grid(axis="y", visible=False)
    _save(fig, path)


def rls_tracking(arx, fixed, t_split, path, truth=None):
    fig, ax = plt.subplots(3, 1, figsize=(9, 7), sharex=True)
    ax[0].plot(arx.info["t"], np.asarray(arx.info["R0"]) * 1e3, color=ORANGE, label=arx.method)
    ax[0].plot(fixed.info["t"], np.asarray(fixed.info["R0"]) * 1e3, color=BLUE, label=fixed.method)
    ax[0].set_ylabel("R0 [mOhm]")
    ax[1].plot(arx.info["t"], np.asarray(arx.info["R1"]) * 1e3, color=ORANGE, label=arx.method)
    ax[1].plot(fixed.info["t"], np.asarray(fixed.info["R1"]) * 1e3, color=BLUE, label=fixed.method)
    ax[1].set_ylabel("R1 [mOhm]")
    ax[2].plot(arx.info["t"], arx.info["tau1"], color=ORANGE, label=arx.method)
    ax[2].axhline(fixed.info["tau1_fixed"], color=BLUE, label=f"{fixed.method} (fixed)")
    ax[2].set_ylabel("tau1 [s]")
    ax[2].set_xlabel("Time [s]")
    ax[2].set_yscale("log")
    if truth is not None:
        for a, v in zip(ax, (truth.R0 * 1e3, truth.R[0] * 1e3, truth.tau[0])):
            a.axhline(v, color=INK, ls=":", lw=1.2, label="true value")
    # Robust y-limits (RLS start-up transients are large).
    for a, key, s in ((ax[0], "R0", 1e3), (ax[1], "R1", 1e3)):
        vals = np.concatenate([np.asarray(arx.info[key]) * s, np.asarray(fixed.info[key]) * s])
        vals = vals[np.isfinite(vals)]
        lo, hi = np.percentile(vals, [1, 99])
        pad = 0.15 * (hi - lo + 1e-9)
        a.set_ylim(lo - pad, hi + pad)
    ax[2].set_ylim(0.1, 1000)
    for a in ax:
        _split(a, t_split)
    ax[2].legend(loc="lower right", fontsize=8, ncol=2)
    ax[0].set_title("RLS parameter tracking (value at the dashed line is reported)")
    _save(fig, path)


def rls_r0_offset(fr, t_split, path, truth_R0=None):
    t = np.asarray(fr.info["t"])
    fig, ax = plt.subplots(2, 1, figsize=(9, 5.5), sharex=True)
    r0 = np.asarray(fr.info["R0"]) * 1e3
    ax[0].plot(t, r0, color=BLUE, label=fr.method)
    if truth_R0 is not None:
        ax[0].axhline(truth_R0 * 1e3, color=INK, ls=":", lw=1.2, label="true value")
    lo, hi = np.percentile(r0[len(r0) // 20:], [1, 99])
    pad = 0.3 * (hi - lo + 1e-6)
    ax[0].set_ylim(lo - pad, hi + pad)
    ax[0].set_ylabel("R0 [mOhm]")
    ax[0].legend(loc="upper right", fontsize=8)
    ax[0].set_title("Online R0 and OCV offset (RLS, runs over the whole profile)")
    ax[1].plot(t, np.asarray(fr.info["offset"]) * 1e3, color=ORANGE, label="offset c = V - OCV - R0*I")
    ax[1].axhline(0, color=INK, lw=0.6)
    ax[1].set_ylim(-30, 30)
    ax[1].set_ylabel("Offset [mV]")
    ax[1].set_xlabel("Time [s]")
    ax[1].legend(loc="upper left", fontsize=8)
    for a in ax:
        _split(a, t_split)
    _save(fig, path)


def tau_scan(sep, n_fit, path, truth_tau=None):
    """Fit error vs time constant. A flat curve means tau is not identifiable."""
    taus = np.asarray(sep.info["tau_grid"])
    rmse = np.sqrt(np.asarray(sep.info["sse_grid"]) / n_fit) * 1e3
    rise = rmse - rmse.min()
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.plot(taus, rise, "-o", color=BLUE, ms=3, label="best fit for each tau")
    ax.axvline(sep.params.tau[0], color=ORANGE, lw=1.4, label=f"selected tau = {sep.params.tau[0]:.1f} s")
    if truth_tau is not None:
        ax.axvline(truth_tau, color=INK, ls=":", lw=1.2, label=f"true tau = {truth_tau:g} s")
    ax.set_ylim(-0.05, max(1.0, rise.max() * 1.1))
    if rise.max() < 0.05:
        ax.text(0.5, 0.55, "Flat curve: the data cannot tell the time constant\n"
                           f"(whole range changes RMSE by only {rise.max() * 1e3:.2f} uV)",
                transform=ax.transAxes, ha="center", color=INK, fontsize=9)
    ax.set_xscale("log")
    ax.set_xlabel("RC time constant tau1 [s]")
    ax.set_ylabel("Fit RMSE above best value [mV]")
    ax.set_title(f"Separable LS: error vs time constant (best RMSE {rmse.min():.2f} mV)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3, fontsize=8)
    _save(fig, path)


def ocv_poly_vs_table(table_ocv, poly_ocv, soc_window, path):
    s = np.linspace(0, 1, 1001)
    fig, ax = plt.subplots(2, 1, figsize=(7, 6), sharex=True,
                           gridspec_kw={"height_ratios": [2, 1]})
    ax[0].plot(table_ocv.soc * 100, table_ocv.ocv, "o", color=GREY, ms=4, label="OCV table")
    ax[0].plot(s * 100, table_ocv(s), color=BLUE, label="linear interpolation (used)")
    ax[0].plot(s * 100, poly_ocv(s), color=ORANGE, ls="--",
               label=f"polynomial order {poly_ocv.order} (MATLAB)")
    ax[0].set_ylabel("OCV [V]")
    ax[0].legend(loc="lower right", fontsize=8)
    ax[0].set_title("OCV model: table interpolation vs polynomial")
    err = (poly_ocv(s) - table_ocv(s)) * 1e3
    ax[1].axvspan(soc_window[0] * 100, soc_window[1] * 100, color=BLUE, alpha=0.12, lw=0,
                  label="SOC range of this test")
    ax[1].plot(s * 100, err, color=ORANGE, label="polynomial - interpolation")
    ax[1].axhline(0, color=INK, lw=0.6)
    ax[1].set_ylim(-40, 40)
    ax[1].set_ylabel("Difference [mV]")
    ax[1].set_xlabel("SOC [%]")
    ax[1].legend(loc="upper right", fontsize=8)
    _save(fig, path)


def rls_lambda_study(t, sens, t_split, path):
    colors = [BLUE, ORANGE, AQUA, YELLOW, "#e87ba4", "#008300"]
    fig, ax = plt.subplots(2, 1, figsize=(9, 6.5))
    for c, (lam, h) in zip(colors, sens["hist_lambda"].items()):
        ax[0].plot(t, h * 1e3, color=c, lw=0.8, label=f"lambda = {lam:g}")
    ax[0].set_ylim(1.5, 3.5)
    ax[0].set_ylabel("R0 estimate [mOhm]")
    ax[0].set_xlabel("Time [s]")
    ax[0].set_title("RLS for R0: effect of the forgetting factor lambda (P0 = 1)")
    ax[0].legend(loc="lower right", ncol=3, fontsize=8)
    _split(ax[0], t_split)
    for c, (P0, h) in zip(colors, sens["hist_P0"].items()):
        ax[1].plot(t, h * 1e3, color=c, lw=1.0, label=f"P0 = {P0:g}")
    ax[1].set_xscale("log")
    ax[1].set_xlim(0.1, t_split)
    ax[1].set_ylim(0, 12)
    ax[1].set_ylabel("R0 estimate [mOhm]")
    ax[1].set_xlabel("Time [s] (log scale, up to the split)")
    ax[1].set_title(f"Effect of the initial covariance P0 (lambda = {sens['lam_for_P0']:g}, start R0 = 10 mOhm)")
    ax[1].legend(loc="upper right", ncol=3, fontsize=8)
    _save(fig, path)


def soc_segments(panels, path):
    """panels: list of (title, rows, truth_fn or None, global_R0 or None)."""
    fig, axs = plt.subplots(1, len(panels), figsize=(5.2 * len(panels), 4), squeeze=False)
    for ax, (title, rows, truth, g) in zip(axs[0], panels):
        used = [r for r in rows if r["used"]]
        for r in used:
            ax.hlines(r["R0"] * 1e3, r["soc_lo"] * 100, r["soc_hi"] * 100, color=BLUE, lw=3)
        ax.plot([], [], color=BLUE, lw=3, label="R0 per SOC bin (segmented LS)")
        lo = min(r["soc_lo"] for r in used)
        hi = max(r["soc_hi"] for r in used)
        if truth is not None:
            s = np.linspace(lo, hi, 50)
            ax.plot(s * 100, truth(s) * 1e3, color=INK, ls=":", lw=1.4, label="true R0(SOC)")
        if g is not None:
            ax.axhline(g * 1e3, color=ORANGE, lw=1.2, label=f"global R0 = {g * 1e3:.3f} mOhm")
        ax.set_xlabel("SOC [%]")
        ax.set_ylabel("R0 [mOhm]")
        ax.set_title(title)
        ax.legend(loc="upper right", fontsize=8)
        vals = [r["R0"] * 1e3 for r in used]
        pad = max(0.1, 0.4 * (max(vals) - min(vals)))
        ax.set_ylim(min(vals) - pad, max(vals) + 2 * pad)
    _save(fig, path)


def pso_convergence(curves, noise_std, path):
    """curves: list of (label, rmse history in mV, colour index)."""
    colors = [BLUE, ORANGE, AQUA, YELLOW, "#e87ba4", "#008300"]
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    for label, h, ci in curves:
        ax.plot(np.arange(len(h)), h, color=colors[ci % len(colors)], lw=1.4, label=label)
    ax.axhline(noise_std * 1e3, color=INK, ls=":", lw=1.2, label=f"noise level {noise_std * 1e3:.1f} mV")
    ax.set_yscale("log")
    ax.set_xlabel("PSO iteration")
    ax.set_ylabel("Best voltage RMSE so far [mV]")
    ax.set_title("PSO convergence")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, fontsize=8)
    _save(fig, path)
