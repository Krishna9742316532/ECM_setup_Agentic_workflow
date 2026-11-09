"""Synthetic test data with known parameters.

Uses the real current profile and OCV table, a known ECM and white voltage
noise with a fixed seed. If a method cannot recover the known parameters
here, its results on real data cannot be trusted.
"""

import numpy as np

from .data import CellData
from .model import ECMParams, simulate


def truth_params(truth: dict) -> ECMParams:
    return ECMParams(R0=truth["R0"], R=tuple(truth["R"]), tau=tuple(truth["tau"]),
                     Q_Ah=truth["Q_Ah"], soc0=truth["soc0"])


def make_synthetic_soc_r0(base: CellData, truth: dict, seed=42, name="synthetic-soc-R0"):
    """Synthetic data where R0 depends on SOC: R0(soc) = R0_ref*(1 + k*(soc_ref - soc)).

    Used to test SOC-segmented LS. Returns (data, truth dict with R0 function).
    """
    from .model import coulomb_count, simulate_varying
    soc = coulomb_count(base.t, base.I, truth["Q_Ah"], truth["soc0"])

    def r0_of_soc(s):
        return truth["R0_ref"] * (1 + truth["k"] * (truth["soc_ref"] - np.asarray(s)))

    v = simulate_varying(base.t, base.I, soc, base.ocv, r0_of_soc(soc),
                         [np.full_like(soc, r) for r in truth["R"]],
                         [np.full_like(soc, tau) for tau in truth["tau"]])
    rng = np.random.default_rng(seed)
    v = v + truth.get("noise_std", 0.0) * rng.standard_normal(len(v))
    return CellData(base.t.copy(), base.I.copy(), v, base.ocv, name), r0_of_soc


def make_synthetic(base: CellData, truth: dict, seed=42, name="synthetic"):
    p = truth_params(truth)
    v, _ = simulate(base.t, base.I, p, base.ocv)
    rng = np.random.default_rng(seed)
    v = v + truth.get("noise_std", 0.0) * rng.standard_normal(len(v))
    return CellData(base.t.copy(), base.I.copy(), v, base.ocv, name), p
