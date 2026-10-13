# axion-fopt-next

A second-generation framework for scalar fields in prescribed stochastic
cosmological first-order phase transitions. The published axion-fopt code is
the reference implementation and is not modified by this project.

## Current state

Implemented: seeded Poisson candidate sampling for a tabulated conformal-time
intensity, false-vacuum exclusion, periodic constant-speed bubble geometry,
first-arrival grid export, and an axion potential/force API independent of any
solver. Only NumPy is needed. This is a correctness-first prototype: rejection
and projection use brute-force geometry. Spatial indexing is a next milestone.

Not yet implemented: CosmoLattice field evolution, physical collision response,
variable bubble speeds, dynamical order parameters, published abundance
reproduction, or minicluster spectra. Overlap means conversion by the earliest
wall; it does not simulate wall/plasma collision dynamics.

## Run

```bash
python -m pip install -e .
python -m unittest discover -s tests -v
python examples/export_background.py
```

Without installing:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
PYTHONPATH=src python examples/export_background.py
```

Tests cover periodic geometry, birth-time causality, exclusion, reproducibility,
Poisson mean/variance, linear-rate time distributions, force derivatives, and
an ensemble comparison to the analytic false-vacuum fraction. A passing test
suite establishes these limited checks, not the published physics benchmark.

See [physics conventions](docs/physics.md), [adapter plan](docs/cosmolattice.md),
and [roadmap](docs/roadmap.md). Licensing is intentionally undecided; no uploaded
source code or third-party solver code is redistributed in this repository.

## C++ adapter prototype

A versioned arrival-map exporter and standalone C++20 reader/force evaluator are
now available. The tests compile the probe and compare its one-bubble results
with Python at multiple stage times. This does not yet evolve fields through
CosmoLattice. See [prototype and source findings](docs/adapter-prototype.md).
