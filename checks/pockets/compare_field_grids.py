"""Fresh paired-grid history comparison using predeclared stiff-wall gates."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from model import history


def compare(a,b,box=None):
    a,b=np.atleast_1d(a),np.atleast_1d(b)
    for h in (a,b):
        if not len(h) or not np.isfinite(h['u']).all() or abs(h['u'][0])>1e-9 or np.any(np.diff(np.round(h['u'],9))<=0):
            raise ValueError('Require finite, increasing histories starting at zero')
    shared,i,j=np.intersect1d(np.round(a['u'],9),np.round(b['u'],9),return_indices=True)
    ea=lambda h:sum(h[k] for k in ('axion_kinetic','axion_gradient','coupling'))
    aa,ab=ea(a),ea(b)
    if aa[0]<=0 or ab[0]<=0:raise ValueError('Require positive initial gas energies')
    dr=float(np.max(np.abs(a['fv_radius'][i]-b['fv_radius'][j])))
    df=float(np.max(np.abs(a['axion_far'][i]/aa[0]-b['axion_far'][j]/ab[0])))
    depletion_difference=float(abs((1-aa[i[-1]]/aa[0])-(1-ab[j[-1]]/ab[0])))
    report=dict(
                final_shared_time=float(shared[-1]),max_radius_difference=dr,
                max_far_fraction_difference=df,depletion_difference=depletion_difference,
                history_gates=dict(radius=dr<=.25,far_fraction=df<=.005,depletion=depletion_difference<=.01,
                                   full_duration=bool(shared[-1]>=80)),
                monitor_retention=dict(coarse=float(a['monitor_energy'][i[-1]]/a['monitor_energy'][0]),
                                       fine=float(b['monitor_energy'][j[-1]]/b['monitor_energy'][0])),
                interpretation='Two-grid history consistency only. Individual wall/accounting gates remain unchanged; no lifetime claim.')
    if 'tracking_energy' in a.dtype.names and 'tracking_energy' in b.dtype.names:
        if box is None or not np.isfinite(box) or box<=0:raise ValueError('Moving centroid comparison requires common positive box')
        def cumulative(h,key,indices):
            # Integrate at each run's own saved cadence, then sample shared times.
            f=h[key];result=np.r_[0.,np.cumsum(.5*(f[1:]+f[:-1])*np.diff(h['u']))]
            return result[indices]/h['tracking_energy'][0]
        storage=float(np.max(np.abs(a['tracking_energy'][i]/a['tracking_energy'][0]-b['tracking_energy'][j]/b['tracking_energy'][0])))
        flux_a=sum(cumulative(a,key,i) for key in ('tracking_phi_flux','tracking_axion_flux'))
        flux_b=sum(cumulative(b,key,j) for key in ('tracking_phi_flux','tracking_axion_flux'))
        current=float(np.max(np.abs(flux_a-flux_b)))
        exchange=float(np.max(np.abs(cumulative(a,'tracking_exchange',i)-cumulative(b,'tracking_exchange',j))))
        delta=np.stack([a['tracking_c'+axis][i]-b['tracking_c'+axis][j] for axis in 'xyz'],-1)
        delta=(delta+box/2)%box-box/2
        centre=float(np.linalg.norm(delta,axis=1).max())
        report['moving_comparison']=dict(max_storage_fraction_difference=storage,
            max_cumulative_current_fraction_difference=current,max_cumulative_exchange_fraction_difference=exchange,
            max_centroid_distance=centre,initial_window_energy_relative_difference=float(abs(a['tracking_energy'][0]-b['tracking_energy'][0])/a['tracking_energy'][0]),
            gates=dict(storage=storage<=.001,current=current<=.001,exchange=exchange<=.001,centroid=centre<=.25),
            interpretation='Additional prerefined-grid gates in MOVING_GRID_COMPARISON.md; do not replace original candidate gates.')
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--coarse',type=Path,required=True)
    p.add_argument('--fine',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    configs=[json.loads((root/'input.json').read_text()) for root in (args.coarse,args.fine)]
    if configs[0]['lSide']!=configs[1]['lSide']:raise ValueError('Grid comparison requires a common box')
    report=compare(history(args.coarse),history(args.fine),configs[0]['lSide'])
    report.update(sources={str((root/'native/history.csv').resolve()):hashlib.sha256((root/'native/history.csv').read_bytes()).hexdigest()
             for root in (args.coarse,args.fine)},analysis_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             moving_contract_sha256=hashlib.sha256(Path(__file__).with_name('MOVING_GRID_COMPARISON.md').read_bytes()).hexdigest())
    (args.output/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
