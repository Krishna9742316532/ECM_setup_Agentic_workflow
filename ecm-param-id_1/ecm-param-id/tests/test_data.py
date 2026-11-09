"""Checks on the raw data: units, sign convention, sampling, OCV table."""

import hashlib

import numpy as np

from conftest import DATA
from ecm.data import (CellData, check_data, current_step_resistance,
                      estimate_noise_std)

EXPECTED_SHA = "c0c2b75d87b5e133bcb6821e860f65a61255b96db4421f008278815563e37097"


def test_data_file_unchanged():
    assert hashlib.sha256(DATA.read_bytes()).hexdigest() == EXPECTED_SHA


def test_shapes_and_sampling(real):
    assert len(real.t) == len(real.I) == len(real.V) == 20005
    assert np.isclose(real.dt, 0.1)
    assert np.all(np.diff(real.t) > 0)


def test_all_checks_pass(real):
    assert check_data(real) == []


def test_voltage_is_single_cell_in_volts(real):
    assert 2.5 < real.V.min() and real.V.max() < 4.3


def test_sign_convention_negative_is_discharge(real):
    # Net charge is negative and the voltage falls -> negative = discharge.
    net_Ah = np.sum(real.I[:-1] * np.diff(real.t)) / 3600
    assert net_Ah < -5
    assert real.V[-200:].mean() < real.V[:200].mean()
    # Voltage jumps in the same direction as the current.
    assert current_step_resistance(real) > 0


def test_flipped_sign_is_detected(real):
    flipped = CellData(real.t, -real.I, real.V, real.ocv, "flipped")
    assert any("sign" in msg for msg in check_data(flipped))


def test_voltage_in_millivolts_is_detected(real):
    mv = CellData(real.t, real.I, real.V * 1000, real.ocv, "mV")
    assert any("Voltage outside" in msg for msg in check_data(mv))


def test_ocv_table(real):
    assert np.isclose(real.ocv.soc[0], 0.0) and np.isclose(real.ocv.soc[-1], 1.0)
    assert np.all(np.diff(real.ocv.ocv) > 0)
    # Inverse lookup is consistent.
    s = np.linspace(0.05, 0.95, 19)
    assert np.allclose(real.ocv.inverse(real.ocv(s)), s, atol=1e-9)


def test_noise_estimate_on_known_noise(synthetic):
    d, _ = synthetic
    assert abs(estimate_noise_std(d.V, d.I) - 0.010) < 0.001


def test_real_noise_level(real):
    assert 0.008 < estimate_noise_std(real.V, real.I) < 0.012
