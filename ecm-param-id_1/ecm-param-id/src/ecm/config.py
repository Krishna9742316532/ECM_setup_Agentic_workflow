"""Project settings: plausible parameter ranges, split time and method settings.

All values are in SI units (Ohm, s, A, V) unless the name says otherwise.
These ranges are the same as in CLAUDE.md. They are used both as optimizer
bounds and as plausibility checks.
"""

from dataclasses import dataclass, field

# Plausible ranges for a large automotive Li-ion cell (about 50 Ah).
RANGES = {
    "R0": (0.5e-3, 10e-3),    # Ohm
    "R": (0.1e-3, 10e-3),     # Ohm, each RC branch
    "tau": (1.0, 200.0),      # s, each RC branch
    "Q_Ah": (30.0, 100.0),    # Ah
    "soc0": (0.40, 0.60),     # -
}

# Fit on t < T_SPLIT, validate on t >= T_SPLIT.
T_SPLIT = 1400.0  # s


@dataclass(frozen=True)
class Settings:
    t_split: float = T_SPLIT
    # Recursive least squares forgetting factor. Memory ~ 1/(1-lam) samples.
    # 0.9998 at 10 Hz -> about 500 s memory.
    rls_lambda: float = 0.9998
    # Block-average length for the ARX methods (10 samples = 1 s).
    arx_block: int = 10
    # Time-constant grid for the separable least squares scan.
    tau_grid_n: int = 60
    # Random seed for anything random (synthetic noise only).
    seed: int = 42
    # OCV models compared in model selection: table interpolation and a
    # 10th-order polynomial of the table (as in the MATLAB course scripts).
    ocv_models: tuple = ("linear", "poly10")
    ranges: dict = field(default_factory=lambda: dict(RANGES))


# Known truth for the synthetic check dataset (same current profile and
# OCV table as the real data). C1 = tau1 / R1 = 20 000 F.
SYNTHETIC_TRUTH = {
    "R0": 2.5e-3,
    "R": (1.5e-3,),
    "tau": (30.0,),
    "Q_Ah": 54.0,
    "soc0": 0.50,
    "noise_std": 0.010,  # V
}

# Synthetic 2RC dataset (tests the 2nd-order ARX methods).
SYNTHETIC_2RC = {
    "R0": 2.5e-3,
    "R": (1.0e-3, 1.5e-3),
    "tau": (5.0, 60.0),
    "Q_Ah": 54.0,
    "soc0": 0.50,
    "noise_std": 0.010,
}

# Synthetic dataset with SOC-dependent R0 (tests SOC-segmented LS).
# R0(soc) = R0_ref * (1 + k*(soc_ref - soc)): 2.5 mOhm at 50 %, ~3.2 at 36 %.
SYNTHETIC_SOC_R0 = {
    "R0_ref": 2.5e-3,
    "k": 2.0,
    "soc_ref": 0.50,
    "R": (1.5e-3,),
    "tau": (30.0,),
    "Q_Ah": 54.0,
    "soc0": 0.50,
    "noise_std": 0.010,
}

# Assumptions used in the MATLAB reference scripts (for comparison only).
MATLAB_ASSUMPTIONS = {"Q_Ah": 60.0, "soc0": 0.5, "ocv_poly_order": 10}
