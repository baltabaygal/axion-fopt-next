"""Run from repo root: PYTHONPATH=src python examples/export_background.py"""
import hashlib
import json
from pathlib import Path
import numpy as np
from axion_fopt_next import sample_catalogue

out = Path("outputs"); out.mkdir(exist_ok=True)
cat = sample_catalogue([0., 1.], [4., 4.], 3., .5, seed=42)
arrival = cat.lattice_arrival_times(32)
np.savez_compressed(out / "background.npz", arrival=arrival, birth_times=cat.times,
                    centres=cat.positions, box_size=cat.box_size, wall_speed=cat.wall_speed)
metadata = {"schema_version": 1, "seed": 42, "rng": "numpy.default_rng/PCG64",
            "numpy_version": np.__version__, "coordinates": "comoving",
            "time": "conformal", "units": "c=1; common inverse-mass unit for time and length",
            "axes": ["x", "y", "z"], "grid": "cell-centred", "n": 32,
            "candidate_count": cat.candidate_count, "accepted_count": len(cat.times),
            "false_vacuum_fraction_at_tau_1": float(np.mean(arrival > 1)),
            "sha256_npz": hashlib.sha256((out / "background.npz").read_bytes()).hexdigest()}
(out / "background.json").write_text(json.dumps(metadata, indent=2) + "\n")
print(json.dumps(metadata, indent=2))
