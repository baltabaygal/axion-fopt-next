"""Prepare frozen physical Fig.2 catalogue and checked arrival maps."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import numba
from run_field import ROOT
from analytic import Case
from tracer_reference import load_paper,simulate,DEFAULT_SOURCE
from spatial_engulfment import accelerated_paper
from axion_fopt_next import Catalogue
from arrival_projection import arrival_map


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--paper-source',type=Path,default=DEFAULT_SOURCE)
    args=parser.parse_args();root=args.output.resolve();root.mkdir(parents=True,exist_ok=False)
    numba.set_num_threads(min(8,numba.config.NUMBA_NUM_THREADS))
    paper,provenance=load_paper(args.paper_source)
    start=time.perf_counter()
    case=Case('fig2',8*np.pi,.15,55.,h=.05,beta=.2)
    tau=np.arange(111)*.5
    cat,internal=simulate(accelerated_paper(paper),case,0,100000,.01,tau)
    seconds=time.perf_counter()-start
    np.savez_compressed(root/'catalogue.npz',times=cat.times,positions=cat.positions,tau=tau,internal_F=internal)
    (root/'catalogue.json').write_text(json.dumps(dict(seconds=seconds,bubbles=len(cat.times),first_birth=float(cat.times.min()),
                   last_birth=float(cat.times.max()),seed=0,tracers=100000,dt=.01,source=provenance),indent=2)+'\n')
    reference=cat.lattice_arrival_times(8);indexed=arrival_map(cat,8)
    error=float(np.max(np.abs(reference-indexed)))
    np.testing.assert_allclose(indexed,reference,atol=2e-12,rtol=0)
    start=time.perf_counter();a=arrival_map(cat,128);projection_seconds=time.perf_counter()-start
    np.save(root/'arrival128.npy',a)
    cells=np.random.default_rng(991).integers(0,128,(256,3))
    direct=cat.arrival_times((cells+.5)*cat.box_size/128)
    spot_error=float(np.max(np.abs(a[tuple(cells.T)]-direct)))
    np.testing.assert_allclose(a[tuple(cells.T)],direct,atol=2e-12,rtol=0)
    (root/'projection.json').write_text(json.dumps(dict(seconds=projection_seconds,N8_max_error=error,N128_spot_error=spot_error),indent=2)+'\n')
    np.save(root/'arrival64.npy',arrival_map(cat,64))
    print(f'Prepared {len(cat.times)} bubbles; oracle errors {error:g}, {spot_error:g}',flush=True)


if __name__=='__main__':main()
