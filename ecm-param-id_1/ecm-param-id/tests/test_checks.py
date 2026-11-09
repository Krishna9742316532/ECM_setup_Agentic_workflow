"""The plausibility and quality checks must flag bad results."""

import numpy as np

from ecm.checks import (fit_quality, metrics, plausibility, rc_identifiability,
                        slow_drift_mV, worst_level)
from ecm.config import RANGES
from ecm.model import ECMParams

GOOD = ECMParams(R0=2.5e-3, R=(1.5e-3,), tau=(30.0,), Q_Ah=54, soc0=0.5)


def test_good_parameters_pass():
    assert worst_level(plausibility(GOOD, RANGES)) == "OK"


def test_negative_resistance_fails():
    assert worst_level(plausibility(GOOD.replace(R=(-1e-4,)), RANGES)) == "FAIL"


def test_out_of_range_tau_fails():
    assert worst_level(plausibility(GOOD.replace(tau=(0.03,)), RANGES)) == "FAIL"


def test_nan_fails():
    assert worst_level(plausibility(GOOD.replace(tau=(float("nan"),)), RANGES)) == "FAIL"


def test_value_on_bound_warns():
    f = plausibility(GOOD.replace(tau=(200.0,)), RANGES)
    assert worst_level(f) == "WARN" and any("tau1" in x.message for x in f)


def test_small_rc_branch_is_flagged():
    t = np.arange(0, 100, 0.1)
    I = 50 * np.sin(t / 5)
    tiny = GOOD.replace(R=(1e-6,))
    assert worst_level(rc_identifiability(tiny, t, I, noise_std=0.01)) == "WARN"


def test_metrics_and_quality():
    V = np.zeros(1000)
    rng = np.random.default_rng(0)
    Vm = 0.01 * rng.standard_normal(1000)
    m = np.ones(1000, dtype=bool)
    mt = metrics(V, Vm, m)
    assert abs(mt["rmse_mV"] - 10) < 1
    assert worst_level(fit_quality(mt, mt, 0.01)) == "OK"
    big = dict(mt, rmse_mV=20.0)
    assert worst_level(fit_quality(mt, big, 0.01)) == "WARN"


def test_slow_drift_detects_trend():
    rng = np.random.default_rng(1)
    white = 0.01 * rng.standard_normal(5000)
    trend = white + np.linspace(0, 0.02, 5000)
    assert slow_drift_mV(white) < 2
    assert slow_drift_mV(trend) > 5
