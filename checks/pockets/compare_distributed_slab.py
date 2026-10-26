"""Fixed-scale coarse/baseline slab comparisons; two levels do not prove order."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from refine_planar_scattering import coarse_sample


def compare(base,case):
    a=json.loads((base/'report.json').read_text());b=json.loads((case/'report.json').read_text())
    if not a['numerical_gates_passed'] or not b['numerical_gates_passed']:raise ValueError('Failed pilot numerics')
    for key in ('box','duration','delta_v','seed','window'):
        if a['contract'][key]!=b['contract'][key]:raise ValueError('Physical contract differs: '+key)
    if a['build']['executable_sha256']!=b['build']['executable_sha256']:raise ValueError('Different solver')
    rows=[]
    for label in ('empty','filled'):
        x=np.loadtxt(base/label/'native/final.txt');y=np.loadtxt(case/label/'native/final.txt')
        if len(x)>len(y):x=coarse_sample(x,len(y))
        elif len(y)>len(x):y=coarse_sample(y,len(x))
        h=np.genfromtxt(base/label/'native/history.csv',delimiter=',',names=True)
        j=np.genfromtxt(case/label/'native/history.csv',delimiter=',',names=True)
        if not np.array_equal(h['t'],j['t']):raise ValueError('History times differ')
        rows.append(dict(label=label,max_scaled_endpoint_field_difference=float(np.max(abs(x-y)/[1.,.04,1.,.02])),
            max_right_wall_position_difference=float(np.max(abs(h['wall_x']-j['wall_x']))),
            final_length_difference=float(2*(j['wall_x'][-1]-h['wall_x'][-1]))))
    return dict(case=str(case),n=b['contract']['n'],dt=b['contract']['dt'],rows=rows,
        report_sha256=hashlib.sha256((case/'report.json').read_bytes()).hexdigest())


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--baseline',type=Path,required=True)
    ap.add_argument('--case',type=Path,nargs='+',required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    r=dict(scope='Two-level sensitivity only; convergence order and fine-grid size gates pending',
        fixed_scales=[1.,.04,1.,.02],baseline=str(args.baseline),
        baseline_report_sha256=hashlib.sha256((args.baseline/'report.json').read_bytes()).hexdigest(),
        cases=[compare(args.baseline,c) for c in args.case],convergence_established=False,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (args.output/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))

if __name__=='__main__':main()
