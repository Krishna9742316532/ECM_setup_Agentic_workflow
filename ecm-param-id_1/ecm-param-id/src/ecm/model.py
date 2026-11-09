"""Equivalent circuit model (Thevenin type) with 0, 1 or 2 RC branches.

    V = OCV(SOC) + R0*I + sum_i V_i
    dV_i/dt = -V_i/(R_i*C_i) + I/C_i,   tau_i = R_i*C_i
    SOC = soc0 + (1/(3600*Q_Ah)) * integral(I dt)

Sign convention: negative current = discharge (V drops below OCV, SOC falls).
The RC branches use the exact zero-order-hold discretisation:
    V_i[k] = a_i*V_i[k-1] + R_i*(1-a_i)*I[k-1],   a_i = exp(-dt/tau_i)
"""

from dataclasses import dataclass, asdict, field

import numpy as np
from scipy.signal import lfilter


@dataclass(frozen=True)
class ECMParams:
    R0: float                 # Ohm
    R: tuple = ()             # Ohm, one per RC branch
    tau: tuple = ()           # s, one per RC branch
    Q_Ah: float = 50.0        # Ah
    soc0: float = 0.5         # -

    def __post_init__(self):
        if len(self.R) != len(self.tau):
            raise ValueError("R and tau must have the same length")

    @property
    def n_rc(self):
        return len(self.R)

    @property
    def C(self):
        """Capacitances [F], C_i = tau_i / R_i."""
        return tuple(t / r if r > 0 else float("inf") for r, t in zip(self.R, self.tau))

    def replace(self, **kw):
        d = asdict(self)
        d.update(kw)
        return ECMParams(**d)

    def to_dict(self):
        """Readable dict with both SI and engineering units."""
        out = {"R0_mOhm": self.R0 * 1e3}
        for i, (r, t, c) in enumerate(zip(self.R, self.tau, self.C), start=1):
            out[f"R{i}_mOhm"] = r * 1e3
            out[f"tau{i}_s"] = t
            out[f"C{i}_F"] = c
        out["Q_Ah"] = self.Q_Ah
        out["soc0"] = self.soc0
        return out


@dataclass
class FitResult:
    """Output of any identification method."""

    method: str
    params: ECMParams
    info: dict = field(default_factory=dict)

    @property
    def valid(self):
        vals = [self.params.R0, *self.params.R, *self.params.tau]
        return all(np.isfinite(v) for v in vals)


def charge_throughput(t, I):
    """Integral of I dt [As], left-rectangle rule (matches zero-order hold)."""
    dt = np.diff(t)
    return np.concatenate([[0.0], np.cumsum(I[:-1] * dt)])


def coulomb_count(t, I, Q_Ah, soc0):
    return soc0 + charge_throughput(t, I) / (3600.0 * Q_Ah)


def rc_voltage(t, I, R, tau):
    """Voltage over one RC branch, starting from 0 V at t[0]."""
    dt = np.diff(t)
    if np.allclose(dt, dt[0], rtol=1e-6, atol=1e-9):
        a = np.exp(-dt[0] / tau)
        return lfilter([0.0, R * (1.0 - a)], [1.0, -a], I)
    # Non-uniform sampling: plain loop.
    v = np.zeros_like(I, dtype=float)
    for k in range(1, len(I)):
        a = np.exp(-dt[k - 1] / tau)
        v[k] = a * v[k - 1] + R * (1.0 - a) * I[k - 1]
    return v


def simulate(t, I, params: ECMParams, ocv):
    """Return (V_model, soc) for the current profile I(t)."""
    soc = coulomb_count(t, I, params.Q_Ah, params.soc0)
    v = ocv(soc) + params.R0 * I
    for r, tau in zip(params.R, params.tau):
        v = v + rc_voltage(t, I, r, tau)
    return v, soc


def simulate_varying(t, I, soc, ocv, R0_k, R_k=(), tau_k=()):
    """Simulate with time-varying parameters (one value per sample).

    Used for SOC-dependent look-up tables. R_k and tau_k are lists of arrays,
    one array per RC branch. With constant arrays this equals `simulate`.
    """
    I = np.asarray(I, dtype=float)
    v = ocv(soc) + np.asarray(R0_k, dtype=float) * I
    dt = np.diff(t)
    for Rb, taub in zip(R_k, tau_k):
        Rb = np.broadcast_to(np.asarray(Rb, dtype=float), I.shape)
        taub = np.broadcast_to(np.asarray(taub, dtype=float), I.shape)
        a = np.exp(-dt / taub[:-1])
        u = Rb[:-1] * (1.0 - a) * I[:-1]
        x = np.zeros_like(I)
        for k in range(1, len(I)):
            x[k] = a[k - 1] * x[k - 1] + u[k - 1]
        v = v + x
    return v
