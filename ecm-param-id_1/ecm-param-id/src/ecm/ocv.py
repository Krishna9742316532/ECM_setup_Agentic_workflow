"""Open-circuit voltage curve OCV(SOC)."""

import numpy as np


class OCVCurve:
    """Piecewise-linear OCV(SOC) lookup.

    The table is sorted so that SOC is ascending (np.interp needs this).
    Outside the table range the value is clamped to the end points; use
    `in_range` to check for this.
    """

    def __init__(self, soc, ocv):
        soc = np.asarray(soc, dtype=float)
        ocv = np.asarray(ocv, dtype=float)
        if soc.shape != ocv.shape or soc.ndim != 1:
            raise ValueError("soc and ocv must be 1-D arrays of the same length")
        order = np.argsort(soc)
        self.soc = soc[order]
        self.ocv = ocv[order]
        if np.any(np.diff(self.soc) <= 0):
            raise ValueError("SOC points in the OCV table must be unique")

    @classmethod
    def from_table(cls, table):
        """Build from an (n, 2) table: column 0 = SOC, column 1 = OCV [V]."""
        table = np.asarray(table, dtype=float)
        return cls(table[:, 0], table[:, 1])

    def __call__(self, soc):
        return np.interp(soc, self.soc, self.ocv)

    def inverse(self, v):
        """SOC for a given OCV. Needs a strictly increasing OCV curve."""
        if np.any(np.diff(self.ocv) <= 0):
            raise ValueError("OCV curve is not strictly increasing; inverse is not unique")
        return np.interp(v, self.ocv, self.soc)

    def slope(self, soc):
        """dOCV/dSOC [V per unit SOC] of the linear segment at `soc`."""
        s = np.diff(self.ocv) / np.diff(self.soc)
        idx = np.clip(np.searchsorted(self.soc, soc) - 1, 0, len(s) - 1)
        return s[idx]

    def in_range(self, soc):
        soc = np.asarray(soc)
        return bool(np.all((soc >= self.soc[0]) & (soc <= self.soc[-1])))


class PolyOCV(OCVCurve):
    """OCV(SOC) as a least-squares polynomial (like MATLAB polyfit(...,10)).

    numpy's Polynomial.fit scales SOC to [-1, 1] internally, so a high order
    is numerically safe. A polynomial is smooth, but it does NOT pass through
    the table points; check `table_error()` before using it.
    """

    def __init__(self, soc, ocv, order=10):
        super().__init__(soc, ocv)
        self.order = order
        self.poly = np.polynomial.Polynomial.fit(self.soc, self.ocv, order)
        self._grid = np.linspace(self.soc[0], self.soc[-1], 4001)
        self._grid_v = self.poly(self._grid)

    @classmethod
    def from_table(cls, table, order=10):
        table = np.asarray(table, dtype=float)
        return cls(table[:, 0], table[:, 1], order)

    def __call__(self, soc):
        return self.poly(np.asarray(soc, dtype=float))

    def inverse(self, v):
        if np.any(np.diff(self._grid_v) <= 0):
            raise ValueError("polynomial OCV is not monotonic; inverse is not unique")
        return np.interp(v, self._grid_v, self._grid)

    def slope(self, soc):
        return self.poly.deriv()(np.asarray(soc, dtype=float))

    def table_error(self, soc_window=None):
        """Max and RMS error at the table points [V] (optionally in a SOC window)."""
        e = self.poly(self.soc) - self.ocv
        if soc_window is not None:
            m = (self.soc >= soc_window[0]) & (self.soc <= soc_window[1])
            e = e[m]
        return {"max_abs_V": float(np.max(np.abs(e))), "rms_V": float(np.sqrt(np.mean(e ** 2)))}

    def is_monotonic(self):
        return bool(np.all(np.diff(self._grid_v) > 0))
