"""Deterministic adapter fixture: PYTHONPATH=src python examples/export_one_bubble.py."""
import hashlib
import json
from pathlib import Path
from axion_fopt_next import Catalogue
from axion_fopt_next.interchange import write_arrival_map

out = Path("outputs")
out.mkdir(exist_ok=True)
cat = Catalogue(4., .7, [1.], [[.2, 1., 2.]])
path = out / "one_bubble.afopt"
write_arrival_map(path, cat.lattice_arrival_times(16), cat.box_size)
metadata = {
    "schema": "AFOPT001", "n": 16, "box_size": 4.,
    "grid": "cell-centred", "axes": ["x", "y", "z"], "ordering": "C; z fastest",
    "coordinates": "comoving", "time": "conformal", "units": "c=1",
    "birth_times": cat.times.tolist(), "centres": cat.positions.tolist(),
    "wall_speed": cat.wall_speed,
    "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
}
(out / "one_bubble.json").write_text(json.dumps(metadata, indent=2) + "\n")
print(path)
