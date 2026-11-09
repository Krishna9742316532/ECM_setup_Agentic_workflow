"""Unit tests of the ECM equations (physics and sign convention)."""

import numpy as np
import pytest

from ecm.model import ECMParams, coulomb_count, rc_voltage, simulate
from ecm.ocv import OCVCurve

OCV = OCVCurve([0.0, 1.0], [3.0, 4.0])   # simple linear OCV: 3 V + 1 V * SOC


def test_zero_current_gives_ocv():
    t = np.arange(0, 100, 0.1)
    p = ECMParams(R0=2e-3, R=(1e-3,), tau=(10.0,), Q_Ah=50, soc0=0.5)
    v, soc = simulate(t, np.zeros_like(t), p, OCV)
    assert np.allclose(v, 3.5) and np.allclose(soc, 0.5)


def test_discharge_lowers_soc_and_voltage():
    t = np.arange(0, 360, 0.1)
    I = -50.0 * np.ones_like(t)       # 50 A discharge
    p = ECMParams(R0=2e-3, Q_Ah=50, soc0=0.5)
    v, soc = simulate(t, I, p, OCV)
    assert soc[-1] < soc[0] and v[-1] < 3.5
    # 50 A for ~360 s = 5 Ah = 10 % of 50 Ah
    assert soc[-1] == pytest.approx(0.4, abs=1e-3)


def test_coulomb_count_units():
    t = np.arange(0, 3600.1, 1.0)
    soc = coulomb_count(t, -10.0 * np.ones_like(t), Q_Ah=100.0, soc0=1.0)
    assert soc[-1] == pytest.approx(0.9, abs=1e-6)   # 10 Ah out of 100 Ah


def test_rc_step_response_matches_analytic():
    t = np.arange(0, 200, 0.1)
    I = np.where(t > 0, -20.0, 0.0)   # step at t = 0+
    R, tau = 1.5e-3, 30.0
    v = rc_voltage(t, I, R, tau)
    # Zero-order hold: the current set at sample 1 acts from t[1] onwards,
    # so the exact solution is shifted by one sample (0.1 s).
    analytic = -20.0 * R * (1 - np.exp(-(t[1:] - t[1]) / tau))
    assert np.allclose(v[1:], analytic, rtol=0, atol=1e-12)
    # Sign: discharge current gives a negative RC voltage (V below OCV).
    assert v[-1] < 0


def test_nonuniform_time_matches_uniform():
    t = np.arange(0, 50, 0.1)
    I = np.sin(t)
    v1 = rc_voltage(t, I, 1e-3, 5.0)
    t2 = t.copy()
    t2[-1] += 1e-3                      # makes sampling non-uniform -> loop path
    v2 = rc_voltage(t2, I, 1e-3, 5.0)
    assert np.allclose(v1[:-1], v2[:-1], atol=1e-12)


def test_capacitance_is_tau_over_r():
    p = ECMParams(R0=1e-3, R=(2e-3,), tau=(40.0,))
    assert p.C[0] == pytest.approx(20000.0)


def test_mismatched_rc_lengths_raise():
    with pytest.raises(ValueError):
        ECMParams(R0=1e-3, R=(1e-3,), tau=())
