"""Parameter recovery on synthetic data with KNOWN parameters.

If a method fails here, its results on real data cannot be trusted.
Truth: R0 = 2.5 mOhm, R1 = 1.5 mOhm, tau1 = 30 s, Q = 54 Ah, SOC0 = 0.5.
"""

import pytest

from conftest import y_and_masks
from ecm.ls import arx_ls, separable_ls
from ecm.nls import nls_fit
from ecm.rls import arx_rls, fixed_tau_rls, r0_rls


def rel(a, b):
    return abs(a / b - 1)


# --- Noise-free data: every method must be (almost) exact ----------------

def test_arx_ls_exact_without_noise(synthetic_clean):
    d, p = synthetic_clean
    y, m, t, I = y_and_masks(d, p)
    r = arx_ls(t[m], y[m], I[m], block=1).params
    assert rel(r.R0, p.R0) < 1e-6 and rel(r.R[0], p.R[0]) < 1e-4 and rel(r.tau[0], p.tau[0]) < 1e-4


def test_separable_ls_exact_without_noise(synthetic_clean):
    d, p = synthetic_clean
    y, m, t, I = y_and_masks(d, p)
    r = separable_ls(t[m], y[m], I[m], 1).params
    assert rel(r.R0, p.R0) < 1e-5 and rel(r.R[0], p.R[0]) < 1e-3 and rel(r.tau[0], p.tau[0]) < 1e-3


def test_arx_rls_converges_without_noise(synthetic_clean):
    d, p = synthetic_clean
    y, m, t, I = y_and_masks(d, p)
    r = arx_rls(t, y, I, lam=1.0, block=1, t_report=1400, P0=1e8).params
    assert rel(r.R0, p.R0) < 1e-3 and rel(r.R[0], p.R[0]) < 1e-2 and rel(r.tau[0], p.tau[0]) < 1e-2


def test_nls_exact_without_noise(synthetic_clean):
    d, p = synthetic_clean
    r = nls_fit(d, d.t < 1400, n_rc=1).params
    assert rel(r.R0, p.R0) < 1e-4 and rel(r.R[0], p.R[0]) < 1e-3
    assert rel(r.tau[0], p.tau[0]) < 1e-3 and rel(r.Q_Ah, p.Q_Ah) < 1e-3
    assert abs(r.soc0 - p.soc0) < 1e-4


# --- 10 mV noise: unbiased methods must stay close ------------------------

def test_nls_with_noise(synthetic):
    d, p = synthetic
    r = nls_fit(d, d.t < 1400, n_rc=1).params
    assert rel(r.R0, p.R0) < 0.01 and rel(r.R[0], p.R[0]) < 0.05
    assert rel(r.tau[0], p.tau[0]) < 0.05 and rel(r.Q_Ah, p.Q_Ah) < 0.01


def test_separable_ls_with_noise(synthetic):
    d, p = synthetic
    y, m, t, I = y_and_masks(d, p)
    r = separable_ls(t[m], y[m], I[m], 1).params
    assert rel(r.R0, p.R0) < 0.01 and rel(r.R[0], p.R[0]) < 0.05 and rel(r.tau[0], p.tau[0]) < 0.05


def test_fixed_tau_rls_with_noise(synthetic):
    d, p = synthetic
    y, m, t, I = y_and_masks(d, p)
    r = fixed_tau_rls(t, y, I, tau=p.tau[0], t_report=1400).params
    assert rel(r.R0, p.R0) < 0.01 and rel(r.R[0], p.R[0]) < 0.05


def test_r0_rls_on_r0_only_data(real):
    """Pure-resistance data: RLS R0 + offset must find R0 and zero offset."""
    from ecm.model import ECMParams, simulate
    p = ECMParams(R0=2.5e-3, Q_Ah=54, soc0=0.5)
    v, soc = simulate(real.t, real.I, p, real.ocv)
    y = v - real.ocv(soc) + 0.004          # add a 4 mV OCV offset
    r = r0_rls(real.t, y, real.I, t_report=1400)
    assert rel(r.params.R0, 2.5e-3) < 1e-3
    assert abs(r.info["offset_at_report_mV"] - 4.0) < 0.05


# --- Known weakness, documented as a test --------------------------------

def test_arx_ls_is_biased_by_voltage_noise(synthetic):
    """Equation-error ARX: noise in y[k-1] pulls tau down (known bias).

    This is why the separable LS / fixed-tau RLS / NLS methods exist.
    """
    d, p = synthetic
    y, m, t, I = y_and_masks(d, p)
    r10hz = arx_ls(t[m], y[m], I[m], block=1).params
    r1s = arx_ls(t[m], y[m], I[m], block=10).params
    assert r10hz.tau[0] < 0.2 * p.tau[0]                 # strongly biased at 10 Hz
    assert 0.2 * p.tau[0] < r1s.tau[0] < p.tau[0]       # block averaging helps, still low


@pytest.mark.parametrize("n_rc", [0, 1, 2])
def test_nls_r0_is_robust_on_real_data(real, n_rc):
    """R0 must not depend on the model structure on the real data."""
    r = nls_fit(real, real.t < 1400, n_rc=n_rc).params
    assert abs(r.R0 * 1e3 - 2.50) < 0.02
