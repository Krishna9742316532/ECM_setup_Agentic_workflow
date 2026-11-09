import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ecm.config import SYNTHETIC_TRUTH  # noqa: E402
from ecm.data import load_mat  # noqa: E402
from ecm.synthetic import make_synthetic  # noqa: E402

DATA = ROOT / "data" / "ECM1.mat"


@pytest.fixture(scope="session")
def real():
    return load_mat(DATA)


@pytest.fixture(scope="session")
def synthetic(real):
    """Synthetic 1RC data with 10 mV noise (fixed seed) and its truth."""
    return make_synthetic(real, SYNTHETIC_TRUTH, seed=42)


@pytest.fixture(scope="session")
def synthetic_clean(real):
    """Same as `synthetic` but without noise."""
    truth = dict(SYNTHETIC_TRUTH, noise_std=0.0)
    return make_synthetic(real, truth, seed=0)


def y_and_masks(d, p, t_split=1400.0):
    """V - OCV(SOC) using the given parameters' Q and SOC0."""
    from ecm.model import coulomb_count
    soc = coulomb_count(d.t, d.I, p.Q_Ah, p.soc0)
    y = d.V - d.ocv(soc)
    m = d.t < t_split
    return y, m, np.asarray(d.t), np.asarray(d.I)
