"""Describe sampled recrossings and physical below-threshold volume."""
import argparse
import json
from pathlib import Path
import numpy as np
from model import read_state
from expanding_model import physical_terms
from build_backend import sha

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True)
    p.add_argument('--case',default='main',choices=['main','space','time-half']);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    native=args.source/args.case/'native';h=np.genfromtxt(native/'history.csv',delimiter=',',names=True)
    m,y=read_state(native/'final.bin');n=m['n'];box=m['box'];a=float(h['a'][-1]);s=y[...,0];mask=s<.5
    values=physical_terms(y,box/n,a);ea=values['K_A']+values['G_A']+values['V_A']
    best=int(np.argmin(h['F_half']));physical=h['F_half']*(h['a']*box)**3;averages=[]
    for start in np.arange(24,float(m['t'])-7.99,8):
        end=start+8;selected=(h['u']>=start-1e-10)&(h['u']<=end+1e-10)
        t=h['t_cosmic'][selected]
        if len(t)>1:averages.append(dict(u_start=float(start),u_end=float(end),cosmic_weighted_physical_volume=float(np.trapezoid(physical[selected],t)/(t[-1]-t[0]))))
    report=dict(passed=bool(abs(float(mask.mean())-h['F_half'][-1])<1e-15 and abs(float(np.mean(abs(s-1)<.1))-h['near_true_fraction'][-1])<1e-15),
        case=args.case,n=n,u_final=m['t'],below_threshold_count=int(mask.sum()),physical_cell_spacing=a*box/n,
        physical_below_threshold_volume=float(mask.sum()*(a*box/n)**3),minimum_sampled_below_count=int(round(float(h['F_half'][best])*n**3)),
        u_at_minimum_below_count=float(h['u'][best]),maximum_sampled_global_minimum=float(h['S_min'].max()),
        axion_energy_fraction_in_below_cells=float(ea[mask].sum()/ea.sum()),physical_volume_averages=averages,
        history_sha256=sha(native/'history.csv'),state_sha256=sha(native/'final.bin'),script_sha256=sha(__file__),
        interpretation='Instantaneous threshold regions and sampled histories; physical volume separates field evolution from dilution by an expanding box. No claim of stable pockets, oscillons, pressure causality or a continuum volume.')
    args.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()
