"""Geometry-only driver with frozen paper classes/kernels; no field evolution."""
import ast
import contextlib
import hashlib
import io
from pathlib import Path
from typing import Any, Dict, Optional, Tuple
import numpy as np
from axion_fopt_next import Catalogue

COMMIT="1b67fecf3d40b707a9dd6f4c2eabefc8f9bb3f1b"
BLOB="ace7473a90f2a88f4e4c5121e6345b86b16a24e3"
DEFAULT_SOURCE=Path("/Users/baltabay/Desktop/miniclusters/axion-fopt/sim_core/axion_sim.py")


def load_paper(source=DEFAULT_SOURCE):
    """Load only four pinned declarations, disabling on-disk Numba JIT caching.

    The original module is never imported: it configures threads and FFTW at
    import time. The source is verified before any declarations are compiled.
    """
    import numba as nb
    raw=Path(source).read_bytes()
    blob=hashlib.sha1(f"blob {len(raw)}\0".encode()+raw).hexdigest()
    if blob!=BLOB:
        raise ValueError("Paper source differs from pinned GitHub blob")
    names={"evolve_bubbles_vectorized","find_surviving_tracers_simple",
           "OptimizedBubbleManager","NucleationTracerManager"}
    tree=ast.parse(raw)
    definitions=[node for node in tree.body if isinstance(node,(ast.FunctionDef,ast.ClassDef)) and node.name in names]
    if {node.name for node in definitions}!=names: raise ValueError("Missing paper declarations")
    for node in definitions:
        if isinstance(node,ast.FunctionDef):
            for decorator in node.decorator_list:
                if isinstance(decorator,ast.Call):
                    for keyword in decorator.keywords:
                        if keyword.arg=="cache": keyword.value=ast.Constant(False)
    module=ast.fix_missing_locations(ast.Module(body=definitions,type_ignores=[]))
    namespace=dict(np=np,nb=nb,FLOAT_TYPE=np.float32,Any=Any,Dict=Dict,Optional=Optional,Tuple=Tuple)
    exec(compile(module,str(source),"exec"),namespace)
    return namespace,dict(commit=COMMIT,git_blob=blob,sha256=hashlib.sha256(raw).hexdigest(),
                          source_path=str(source),numba=nb.__version__,jit_disk_cache=False)


def simulate(paper, case, seed, ntracers, outer_dt, times):
    """Keep the paper's left-end rate, birth batching and end-step engulfment.

    outer_dt is imposed here instead of invoking the field's CFL controller.
    This isolates its influence; this is not an unmodified coupled paper run.
    """
    saved_rng=np.random.get_state()
    try:
        np.random.seed(seed)
        with contextlib.redirect_stdout(io.StringIO()):
            tracers=paper["NucleationTracerManager"](ntracers,case.box)
        bubbles=paper["OptimizedBubbleManager"](case.box)
        internal=[1.]
        tau=0.
        sampled=[0.]
        attempted=0
        while tau<case.end-1e-12:
            dt=min(outer_dt,case.end-tau)
            sub=0.
            while sub<dt-1e-14:
                step=dt-sub
                if tracers.num_active_tracers>0:
                    p=float(case.rate(tau+sub))*case.box**3/ntracers*step
                    if p>.1:
                        step*=.1/p; p=.1
                    if p>0.:
                        coords=tracers.get_active_tracer_coords()
                        mask=np.random.rand(len(coords)).astype(np.float32)<p
                        attempted+=int(mask.sum())
                        for centre in coords[mask]: bubbles.add_bubble(centre,case.speed,tau+sub)
                        tracers.active_indices=tracers.active_indices[~mask]
                sub+=step
            tau+=dt
            bubbles.evolve_bubbles(tau)
            centres,radii=bubbles.get_centers_and_radii_comoving()
            tracers.update_state(centres,radii)
            sampled.append(tau); internal.append(tracers.get_false_vacuum_fraction())
        catalogue=Catalogue(case.box,case.speed,bubbles.birth_times[:bubbles.count],bubbles.centers[:bubbles.count],attempted)
        # Tracer snapshots are step-end observables. Interpolate only for display;
        # shared requested times are checked to coincide with outer-step endpoints.
        curve=np.interp(times,sampled,internal)
        for t in times:
            if np.min(np.abs(np.asarray(sampled)-t))>1e-9:
                raise ValueError("Requested F(t) must be a tracer step endpoint")
        return catalogue,curve
    finally:
        np.random.set_state(saved_rng)
