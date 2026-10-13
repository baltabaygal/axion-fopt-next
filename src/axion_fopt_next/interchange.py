"""Versioned, solver-independent arrival-map interchange (little endian)."""
from pathlib import Path
import struct
import numpy as np

MAGIC = b"AFOPT001"


def write_arrival_map(path, arrival, box_size):
    """Write cell-centred cubic (x,y,z) map, C order (z fastest).

    Values and box_size use the same conformal/comoving inverse-mass unit.
    Positive infinity means never converted. Negative finite times are allowed.
    """
    arrival = np.asarray(arrival, dtype=float)
    if arrival.ndim != 3 or len(set(arrival.shape)) != 1 or arrival.shape[0] < 1:
        raise ValueError("Expected a nonempty cubic arrival map")
    if np.any(np.isnan(arrival)) or np.any(np.isneginf(arrival)):
        raise ValueError("Arrival times must be finite or positive infinity")
    if not np.isfinite(box_size) or box_size <= 0:
        raise ValueError("box_size must be positive and finite")
    n = arrival.shape[0]
    if n > 2**32 - 1:
        raise ValueError("Grid dimension exceeds format limit")
    with Path(path).open("wb") as stream:
        stream.write(MAGIC + struct.pack("<Id", n, box_size))
        stream.write(arrival.astype("<f8").tobytes(order="C"))
