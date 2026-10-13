"""Cheap ensemble validation in static space; no lattice field evolution."""
import json
from pathlib import Path
import numpy as np
from axion_fopt_next import sample_catalogue

probe_seed = 902
probes = np.random.default_rng(probe_seed).uniform(0, 3, (1000, 3))
fractions, counts = [], []
for seed in range(100):
    cat = sample_catalogue([0., 1.], [4., 4.], 3., .5, seed)
    fractions.append(float(np.mean(~cat.phase(1., probes))))
    counts.append(cat.candidate_count)
theory = float(np.exp(-np.pi*4*.5**3/3))
stderr = float(np.std(fractions, ddof=1)/10)
report = {"status": "provisional numeric limit check", "model": "constant q; static background",
          "box_size": 3, "wall_speed": .5, "tau_end": 1, "q": 4,
          "catalogue_seeds": [0, 99], "probe_seed": probe_seed, "probes_per_realization": 1000,
          "realizations": 100, "numpy_version": np.__version__,
          "analytic_P_FV": theory, "measured_P_FV": float(np.mean(fractions)),
          "ensemble_standard_error": stderr, "difference_in_standard_errors": float((np.mean(fractions)-theory)/stderr),
          "expected_candidates": 108, "mean_candidates": float(np.mean(counts))}
Path("outputs").mkdir(exist_ok=True)
Path("outputs/statistics.json").write_text(json.dumps(report, indent=2)+"\n")
print(json.dumps(report, indent=2))
