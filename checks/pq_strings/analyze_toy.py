"""Quantitative finite-box toy observables; never a paper abundance comparison."""
import argparse,json,csv
from pathlib import Path
import numpy as np
from model import read_state,defects,radial_coupling,energies
from spectra import currents,mode_spectrum,normalized_yield,coherent_zero_mode

def analyze(root):
    cfg=json.loads((root/'run/input.json').read_text());box=cfg['lSide'];star=cfg['tau_star'];origin=cfg['tau_origin']
    final_meta,_=read_state(root/'run/native/final.bin')
    if abs(final_meta['u']-cfg['tMax'])>1e-9:raise ValueError('Incomplete run')
    rows=[]
    paths=sorted((root/'run/native').glob('frame_*.bin'))
    # Spatial refinement needs the same fixed late window, not repeated early FFTs.
    if cfg['N']>=128:paths=paths[-int(round(60/(cfg['dt']*cfg['snapshot_every'])))-1:]
    for file in paths:
        meta,s=read_state(file);d,_=defects(s,box,meta['tau']);ma=np.sqrt(meta['chi']);ms=np.sqrt(2*radial_coupling(cfg,meta['u']))
        j0,ji,parts=currents(s,box,meta['a']);e=energies(s,box,meta['a'],radial_coupling(cfg,meta['u']),meta['chi'])
        row=dict(tau=meta['tau'],a=meta['a'],plaquettes=d['plaquettes'],walls=d['wall_links'],
          core_cells=1/(meta['a']*ms*(box/cfg['N'])),wall_cells=None if ma==0 else 1/(meta['a']*ma*(box/cfg['N'])),
          theta_max=float(np.max(abs(np.arctan2(s[...,1],s[...,0])))),radial_rms=float(np.sqrt(np.mean((np.linalg.norm(s[...,:2],axis=-1)-1)**2))),
          rho_total=e['total'],rho_tilt=e['tilt'],rho_angular_kinetic=parts['kinetic_angular'],rho_angular_gradient=parts['gradient_angular'])
        zero=coherent_zero_mode(s,meta['a'],ma) if ma>0 else dict(eligible=False,number_density=None)
        row['harmonic_eligible']=bool(zero['eligible'] and not d['plaquettes'] and not d['wall_links'])
        for scheme in ('continuum','lattice'):
            summary=mode_spectrum(j0,ji,box,meta['a'],ma,dispersion=scheme)['summary']
            nonzero=summary['number_density'];n0=zero['number_density'] if row['harmonic_eligible'] else None
            row['nu_nonzero_'+scheme]=normalized_yield(nonzero,meta['a'],star/origin,origin/star**2)
            row['nu_total_'+scheme]=None if n0 is None else normalized_yield(nonzero+n0,meta['a'],star/origin,origin/star**2)
            if scheme=='lattice':row['zero_number_fraction']=None if n0 is None else n0/(nonzero+n0);row['rho_axion_harmonic']=None if n0 is None else summary['harmonic_axion_energy_proxy']+zero['energy_density']
        rows.append(row)
    tail=[]
    for r in reversed(rows):
        if not r['harmonic_eligible']:break
        tail.append(r)
    tail.reverse()
    # Fixed late window declared in the run contract: final 60 conformal units.
    late=[r for r in tail if r['tau']>=rows[-1]['tau']-60]
    v=np.array([r['nu_total_lattice'] for r in late],dtype=float)
    phase=0. if len(late)<2 else float(np.trapezoid([r['a']*cfg['mass_cap'] for r in late],[r['tau'] for r in late]))
    drift=None if len(v)<2 else float(np.ptp(v)/np.mean(v))
    topology=json.loads((root/'topology.json').read_text())
    defect_rows=[r for r in topology if r['plaquettes'] or r['walls']]
    minimum_core=min(r['core_cells'] for r in defect_rows) if defect_rows else None
    summary=dict(paper_reproduction=False,model='smooth capped susceptibility, finite periodic volume',
       collapse_bracket=None,minimum_core_cells_with_defects=minimum_core,
       late_window=[r['tau'] for r in late],late_phase_advance=phase,late_yield_peak_to_peak=drift,
       late_mean_nu_toy=None if not len(v) else float(np.mean(v)),
       late_half_window_mean_change=None if len(v)<4 else float(abs(np.mean(v[:len(v)//2])-np.mean(v[len(v)//2:]))/np.mean(v)),
       final=rows[-1],gates=dict(defects_resolved=minimum_core is not None and minimum_core>=4,
       harmonic_tail=len(late)>1,two_mass_periods=bool(phase>=4*np.pi),plateau=drift is not None and drift<=.05,
       temporal_refinement=False,spatial_refinement=False))
    summary['gates']={k:bool(v) for k,v in summary['gates'].items()}
    for i in range(1,len(topology)):
        if topology[i-1]['plaquettes'] and not topology[i]['plaquettes']:
            summary['collapse_bracket']=[topology[i-1]['tau'],topology[i]['tau']]
    (root/'toy_analysis.json').write_text(json.dumps(dict(summary=summary,rows=rows),indent=2))
    with (root/'toy_history.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    print(json.dumps(summary,indent=2));return summary
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source',type=Path);a=p.parse_args();analyze(a.source)
