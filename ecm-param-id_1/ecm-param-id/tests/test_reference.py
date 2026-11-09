"""Tests for the methods ported from the MATLAB reference scripts."""

import numpy as np
import pytest

from conftest import y_and_masks
from ecm.config import SYNTHETIC_2RC, SYNTHETIC_SOC_R0
from ecm.data import CellData, check_soc_direction
from ecm.ls import (arx2_from_physical, arx2_ls, arx2_to_physical, difference_ls,
                    separable_ls)
from ecm.model import ECMParams, coulomb_count, simulate, simulate_varying
from ecm.nls import nls_fit
from ecm.ocv import PolyOCV
from ecm.rls import r0_rls_sensitivity
from ecm.segmented import lookup, segmented_ls
from ecm.synthetic import make_synthetic, make_synthetic_soc_r0


def rel(a, b):
    return abs(a / b - 1)


# --- OCV polynomial ---------------------------------------------------------

def test_poly_ocv_is_smooth_but_misses_table(real):
    p = PolyOCV(real.ocv.soc, real.ocv.ocv, 10)
    assert p.is_monotonic()
    err = p.table_error()
    assert 0.010 < err["max_abs_V"] < 0.020        # up to ~15 mV off the table
    s = np.linspace(0.1, 0.9, 9)
    assert np.allclose(p.inverse(p(s)), s, atol=1e-3)


def test_real_data_fits_poly_ocv_with_stable_capacity(real):
    """With the polynomial OCV, Q no longer depends on the fit window."""
    dp = CellData(real.t, real.I, real.V, PolyOCV(real.ocv.soc, real.ocv.ocv, 10))
    qs = [nls_fit(dp, m, n_rc=0).params.Q_Ah
          for m in (real.t < 1400, real.t >= 600, np.ones_like(real.t, bool))]
    assert max(qs) - min(qs) < 0.5
    assert all(abs(q - 60.0) < 0.5 for q in qs)


# --- Coulomb-counting sign (MATLAB Ex3, 2RC part) ---------------------------

def test_flipped_soc_sign_is_detected(real):
    ok = coulomb_count(real.t, real.I, 60, 0.5)
    flipped = coulomb_count(real.t, -real.I, 60, 0.5)
    assert check_soc_direction(real.I, ok) == []
    assert check_soc_direction(real.I, flipped) != []


# --- Difference method (MATLAB Ex3 Method 2) --------------------------------

def test_difference_ls_exact_on_r0_only_data(real):
    p = ECMParams(R0=2.5e-3, Q_Ah=60, soc0=0.5)
    v, _ = simulate(real.t, real.I, p, real.ocv)
    r = difference_ls(real.t, v, real.I)
    assert rel(r.params.R0, 2.5e-3) < 1e-3
    a, b = r.info["two_column_theta"]
    assert rel(a, 2.5e-3) < 2e-3 and rel(b, 2.5e-3) < 2e-3


def test_difference_ls_independent_of_soc_setup(real):
    """No Q, SOC0 or OCV needed: result is near the NLS R0 (2.50 mOhm)."""
    r = difference_ls(real.t, real.V, real.I)
    assert abs(r.params.R0 * 1e3 - 2.50) < 0.03


# --- 2nd-order ARX (MATLAB Ex3, 2RC part) -----------------------------------

def test_arx2_mapping_roundtrip():
    th = arx2_from_physical(2.5e-3, (1e-3, 1.5e-3), (5.0, 60.0), 0.1)
    R0, R, tau = arx2_to_physical(th, 0.1)
    assert R0 == pytest.approx(2.5e-3)
    assert R == pytest.approx((1e-3, 1.5e-3), rel=1e-6)
    assert tau == pytest.approx((5.0, 60.0), rel=1e-6)


def test_arx2_complex_poles_give_nan():
    R0, R, tau = arx2_to_physical([1.0, -0.5, 1e-3, 0, 0], 0.1)
    assert np.isnan(R[0]) and np.isnan(tau[0])


def test_arx2_ls_exact_without_noise(real):
    d, p = make_synthetic(real, dict(SYNTHETIC_2RC, noise_std=0.0), seed=0)
    y, m, t, I = y_and_masks(d, p)
    r = arx2_ls(t[m], y[m], I[m], block=1).params
    assert rel(r.R0, p.R0) < 1e-6
    assert r.R == pytest.approx(p.R, rel=1e-4) and r.tau == pytest.approx(p.tau, rel=1e-4)


def test_arx2_ls_breaks_down_with_noise(real):
    """Known limitation: 10 mV noise makes the ARX2 poles non-physical."""
    d, p = make_synthetic(real, SYNTHETIC_2RC, seed=42)
    y, m, t, I = y_and_masks(d, p)
    r = arx2_ls(t[m], y[m], I[m], block=1)
    assert not r.valid


def test_nls_2rc_recovers_2rc_truth(real):
    d, p = make_synthetic(real, SYNTHETIC_2RC, seed=42)
    r = nls_fit(d, d.t < 1400, n_rc=2).params
    assert rel(r.R0, p.R0) < 0.01
    assert r.tau == pytest.approx(p.tau, rel=0.05) and r.R == pytest.approx(p.R, rel=0.05)


# --- Time-varying simulation and SOC-segmented LS (MATLAB Ex4) -------------

def test_simulate_varying_equals_simulate_for_constant_params(real):
    p = ECMParams(2.5e-3, (1.5e-3,), (30.0,), 54, 0.5)
    v1, soc = simulate(real.t, real.I, p, real.ocv)
    n = len(soc)
    v2 = simulate_varying(real.t, real.I, soc, real.ocv, np.full(n, p.R0),
                          [np.full(n, p.R[0])], [np.full(n, p.tau[0])])
    assert np.max(np.abs(v1 - v2)) < 1e-12


def test_segmented_ls_recovers_soc_dependent_r0(real):
    d, r0_true = make_synthetic_soc_r0(real, SYNTHETIC_SOC_R0, seed=42)
    T = SYNTHETIC_SOC_R0
    soc = coulomb_count(d.t, d.I, T["Q_Ah"], T["soc0"])
    y = d.V - d.ocv(soc)
    edges = np.arange(0.35, 0.5251, 0.025)
    rows = segmented_ls(d.t, y, d.I, soc, edges, tau=T["tau"][0], mask=d.t < 1400)
    used = [r for r in rows if r["used"]]
    assert len(used) >= 4
    for r in used:
        assert rel(r["R0"], r0_true((r["soc_lo"] + r["soc_hi"]) / 2)) < 0.03
    # Look-up table returns the bin value inside a bin.
    s_mid = (used[0]["soc_lo"] + used[0]["soc_hi"]) / 2
    assert lookup(rows, np.array([s_mid]), "R0")[0] == pytest.approx(used[0]["R0"])


def test_segmented_ls_skips_empty_bins(real):
    soc = coulomb_count(real.t, real.I, 60, 0.5)
    rows = segmented_ls(real.t, real.V - real.ocv(soc), real.I, soc, [0.9, 0.95, 1.0])
    assert not any(r["used"] for r in rows)


# --- RLS lambda / P0 study --------------------------------------------------

def test_smaller_lambda_gives_noisier_r0(real, synthetic):
    d, p = synthetic
    y, m, t, I = y_and_masks(d, p)
    # R0-only regression on 1RC data still shows the variance trend.
    s = r0_rls_sensitivity(t, y, I, lams=(0.97, 0.99, 0.9998), P0s=(1.0,), t_report=1400)
    stds = [r["R0_std_mOhm"] for r in s["lambda"]]
    assert stds[0] > stds[1] > stds[2]


def test_tiny_p0_converges_slowly(real):
    soc = coulomb_count(real.t, real.I, 60, 0.5)
    y = real.V - PolyOCV(real.ocv.soc, real.ocv.ocv, 10)(soc)
    s = r0_rls_sensitivity(real.t, y, real.I, lams=(1.0,), P0s=(1e-7, 1.0), t_report=1400)
    t_small, t_big = (r["t_reach_5pct_s"] for r in s["P0"])
    assert t_small > 10 * t_big


def test_separable_ls_r0_only_on_poly_setup(real):
    soc = coulomb_count(real.t, real.I, 60, 0.5)
    y = real.V - PolyOCV(real.ocv.soc, real.ocv.ocv, 10)(soc)
    m = real.t < 1400
    r = separable_ls(real.t[m], y[m], real.I[m], 0).params
    assert abs(r.R0 * 1e3 - 2.50) < 0.01


# --- PSO (PSO_battery_parameter.m) -----------------------------------------

def test_pso_minimize_finds_simple_minimum():
    from ecm.pso import pso_minimize
    x, f, hist = pso_minimize(lambda x: float(np.sum((x - [1.0, -2.0]) ** 2)),
                              [-5, -5], [5, 5], iters=80, seed=0)
    assert np.allclose(x, [1.0, -2.0], atol=1e-2)
    assert np.all(np.diff(hist) <= 0)            # best cost never gets worse


def test_pso_recovers_synthetic_1rc(synthetic):
    from ecm.pso import pso_fit
    d, p = synthetic
    r = pso_fit(d, d.t < 1400, n_rc=1, iters=60, seed=1, polish=True).params
    assert rel(r.R0, p.R0) < 0.01 and rel(r.R[0], p.R[0]) < 0.05 and rel(r.tau[0], p.tau[0]) < 0.05


def test_pso_is_reproducible_with_same_seed(synthetic):
    from ecm.pso import pso_fit
    d, _ = synthetic
    a = pso_fit(d, d.t < 1400, n_rc=0, iters=10, seed=7, polish=False).params
    b = pso_fit(d, d.t < 1400, n_rc=0, iters=10, seed=7, polish=False).params
    assert a == b


def test_matlab_pso_setup_gets_stuck_on_bounds(synthetic):
    """Linear search over [1 mOhm, 0.5 Ohm] and C1 <= 5000 F: the swarm
    collapses on the lower bounds and misses the true parameters."""
    from ecm.pso import matlab_style_pso
    d, p = synthetic
    m = d.t < 1400
    r = matlab_style_pso(CellData(d.t[m], d.I[m], d.V[m], d.ocv), Q_Ah=p.Q_Ah, soc0=p.soc0)
    assert r.info["at_bound"]                     # result sits on a bound
    assert r.info["rmse_mV"] > 2 * 10.0           # far above the 10 mV noise


# --- End-to-end: the reference study runs and writes all outputs ----------

def test_reference_study_runs(tmp_path):
    from conftest import DATA
    from ecm.reference_methods import run
    R = run(tmp_path, DATA)
    for f in ("summary.md", "results.json", "timeseries_ECM1.csv", "ocv_poly_vs_table.png",
              "rls_lambda_P0.png", "soc_segments.png", "pso_convergence.png"):
        assert (tmp_path / f).exists(), f
    assert R["A_ocv_model"]["selected_by_pipeline"] == "poly10"
    text = (tmp_path / "summary.md").read_text()
    assert "## G. Particle swarm" in text and "## H. CSV export" in text
