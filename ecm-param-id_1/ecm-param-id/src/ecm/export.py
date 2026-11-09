"""Export time series to CSV (for example as input data for an ANN).

Note: identified parameters are ESTIMATES with known biases and limits, not
ground truth. The CSV header says how each column was produced.
"""

import numpy as np


def export_timeseries_csv(path, columns: dict, notes=()):
    """Write equal-length 1-D arrays as CSV. Lines in `notes` go in the header."""
    names = list(columns)
    data = np.column_stack([np.asarray(columns[k], dtype=float) for k in names])
    header = "\n".join([*(f"# {n}" for n in notes), ",".join(names)])
    np.savetxt(path, data, delimiter=",", header=header, comments="", fmt="%.8g")
    return path
