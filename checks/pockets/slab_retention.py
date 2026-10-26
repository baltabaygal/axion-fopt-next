"""Wall-following endpoint windows; no cumulative leakage without histories."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from planar_diagnostics import wall_geometry
from planar_stencils import window_budget
from paper_reference import Parameters


def weights(x,box,centre,radius,centre_velocity,radius_velocity,padding=4.):
    if not 0<padding<radius or radius+padding>=box/2:raise ValueError('Invalid slab window')
    d=(x-centre+box/2)%box-box/2
    q=np.clip((abs(d)-radius+padding)/(2*padding),0,1)
    w=1-q**3*(10-15*q+6*q*q)
    wd=30*q*q*(1-q)**2*(np.sign(d)*centre_velocity+radius_velocity)/(2*padding)
    return w,wd


def measure(state,box,p,padding=4.):
    right=wall_geometry(state,box);reflected=wall_geometry(state[::-1],box)
    left=box-reflected['position'];lv=-reflected['phase_velocity'];rv=right['phase_velocity']
    centre=(left+right['position'])/2;radius=(right['position']-left)/2
    x=(np.arange(len(state))+.5)*box/len(state)
    w,wd=weights(x,box,centre,radius,(lv+rv)/2,(rv-lv)/2,padding)
    inside=window_budget(state,box,w,wd,p,4)
    total=window_budget(state,box,np.ones(len(state)),np.zeros(len(state)),p,4)['axion_energy']
    return dict(left=left,right=right['position'],radius=radius,centre=centre,
        centre_velocity=(lv+rv)/2,radius_velocity=(rv-lv)/2,padding=padding,
        inside_axion_energy=inside['axion_energy'],outside_axion_energy=total-inside['axion_energy'],
        outside_fraction_of_current_energy=(total-inside['axion_energy'])/total if total else 0.,
        instantaneous_net_outward_axion_flux=inside['axion_flux'],
        instantaneous_lab_flux=inside['axion_lab_flux'],instantaneous_advection=inside['axion_advection'])


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--case',type=Path,nargs='+',required=True)
    ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    rows=[]
    for case in args.case:
        report=json.loads((case/'report.json').read_text())
        if not report['numerical_gates_passed']:raise ValueError('Failed source case')
        p=Parameters(dv=report['contract']['delta_v']);row=dict(case=str(case),
            report_sha256=hashlib.sha256((case/'report.json').read_bytes()).hexdigest())
        for endpoint in ('initial','final'):
            state=np.loadtxt(case/'filled/native'/f'{endpoint}.txt')
            row[endpoint]=[measure(state,report['contract']['box'],p,padding) for padding in (2.,4.,8.)]
        rows.append(row)
    result=dict(scope='Endpoint wall-following energy and instantaneous flux only; no cumulative escape or lifetime inference',
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),rows=rows,solver_runs=0)
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    for row in rows:print(row['case'],[(r['padding'],r['outside_fraction_of_current_energy']) for r in row['final']])

if __name__=='__main__':main()
