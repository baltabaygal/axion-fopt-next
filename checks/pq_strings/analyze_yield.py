"""Audit saved PQ states with current spectra; never label a transient as relic yield."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import numpy as np
from model import ROOT,read_state,defects,energies,radial_coupling
from spectra import currents,mode_spectrum,shell_spectrum,normalized_yield

CONTRACT=dict(
    source='arXiv:2609.20801v1 Eqs.7,21,24; arXiv:2512.13653v1 Eqs.80-81',
    current='J_mu=Im(conj(Phi/|Phi|)*d_mu Phi); midpoint-link spatial discretization',
    zero_policy='raise on exact zero projections; homogeneous k=0 excluded and reported',
    primary='continuum frequency/denominator with lattice K projection',
    companion='fully lattice frequency/denominator for stencil normalization',
    numerical_gates=dict(relative_parseval=2e-12,relative_partition=2e-12),
    comparison_requirements=dict(minimum_core_cells=4.,minimum_wall_cells=4.,
        maximum_ma_ms=.16,maximum_kmin_tau_star=1.5,minimum_tau_ratio=3.,
        maximum_adiabatic_parameter=.1,minimum_post_defect_phase_advance=4*np.pi),
    interpretation='Existing short, finite-box runs are estimator diagnostics, not a reproduction claim')


def write_csv(path,columns):
    with path.open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(columns.keys());writer.writerows(zip(*columns.values()))


def audit(source,output):
    cfg=json.loads((source/'input.json').read_text())
    output.mkdir(parents=True,exist_ok=False)
    a_star=cfg['tau_star']/cfg['tau_origin'];h_star=cfg['tau_origin']/cfg['tau_star']**2
    rows=[];hashes={};checks=[]
    paths=sorted((source/'native').glob('frame_*.bin'))
    if not paths:
        paths=[source/'native/initial.bin',source/'native/final.bin']
    for path in paths:
        meta,state=read_state(path);hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        d,_=defects(state,meta['box'],meta['tau'])
        j0,ji,partition=currents(state,meta['box'],meta['a'])
        e=energies(state,meta['box'],meta['a'],radial_coupling(cfg,meta['u']),meta['chi'])
        row=dict(snapshot=path.name,lambda_actual=radial_coupling(cfg,meta['u']),**meta,plaquettes=d['plaquettes'],wall_links=d['wall_links'],
                 theta_max=float(np.max(np.abs(np.angle(state[...,0]+1j*state[...,1])))),
                 radial_rms_from_vacuum=float(np.sqrt(np.mean((np.linalg.norm(state[...,:2],axis=-1)-1)**2))),
                 kmin_tau_star=2*np.pi/meta['box']*cfg['tau_star'],tau_ratio=meta['tau']/cfg['tau_star'],
                 rho_pq_total=e['total'],rho_radial_potential=e['radial'],rho_tilt=e['tilt'],**partition)
        denominator=max(e['kinetic']+e['gradient'],1e-30)
        partition_error=max(abs(partition['kinetic_total']-e['kinetic']),
                            abs(partition['gradient_total']-e['gradient']),
                            abs(partition['kinetic_unassigned']),abs(partition['gradient_unassigned']))/denominator
        row['partition_relative_error']=partition_error
        numerical=partition_error<CONTRACT['numerical_gates']['relative_partition']
        ms=np.sqrt(2*radial_coupling(cfg,meta['u']));ma=np.sqrt(meta['chi']);dx=meta['box']/meta['n']
        row.update(core_cells=1/(meta['a']*ms*dx),wall_cells=None if ma==0 else 1/(meta['a']*ma*dx),ma_ms=ma/ms)
        for scheme in ('continuum','lattice'):
            modes=mode_spectrum(j0,ji,meta['box'],meta['a'],ma,dispersion=scheme)
            shell=shell_spectrum(modes,meta['box'])
            shell['k_tau_star']=shell['k']*cfg['tau_star']
            shell['nu_per_dlnk']=(meta['a']/a_star)**3*shell['number_per_dlnk']/h_star
            write_csv(output/f"{path.stem}_{scheme}.csv",shell)
            summary=modes['summary']
            row[scheme]=summary
            row[scheme]['nu_diagnostic']=normalized_yield(summary['number_density'],meta['a'],a_star,h_star)
            parseval_scale=max(float(np.mean(j0*j0)),float(np.mean(np.sum(ji*ji,axis=0))),1e-30)
            residual=max(summary['temporal_parseval_error'],summary['spatial_parseval_error'])/parseval_scale
            row[scheme]['parseval_relative_error']=residual
            shell_error=abs(float(np.sum(shell['number_per_dlnk']*shell['dlnk']))-summary['number_density'])/max(summary['number_density'],1e-30)
            row[scheme]['shell_sum_relative_error']=shell_error
            numerical=numerical and residual<CONTRACT['numerical_gates']['relative_parseval'] and shell_error<2e-12
        fundamental=2*np.pi/meta['box']
        omega0=meta['a']*ma;omega_min=np.hypot(fundamental,omega0)
        row['mass_phase_frequency']=omega0
        row['adiabatic_parameter_kmin']=(1+cfg['n_tau']/2)*omega0**2/(meta['tau']*omega_min**3)
        rows.append(row);checks.append(bool(numerical))
    last=rows[-1];criteria=CONTRACT['comparison_requirements']
    defect_free_tail=[]
    for row in reversed(rows):
        if row['plaquettes'] or row['wall_links']:break
        defect_free_tail.append(row)
    defect_free_tail.reverse()
    phase_advance=0. if len(defect_free_tail)<2 else float(np.trapezoid(
        [r['mass_phase_frequency'] for r in defect_free_tail],[r['tau'] for r in defect_free_tail]))
    qualification=dict(
        paper_mass_law=cfg.get('mass_schedule','powerlaw')=='powerlaw',
        physical_coupling_stage=last['tau']>=cfg.get('core_growth_until',0.),
        growing_mass_eq_normalized=cfg['chi_star']>0 and abs(cfg['chi_star']-h_star**2)<1e-14,
        after_paper_collapse_epoch=last['tau_ratio']>=criteria['minimum_tau_ratio'],
        defects_absent=not last['plaquettes'] and not last['wall_links'],
        two_post_defect_mass_periods=phase_advance>=criteria['minimum_post_defect_phase_advance'],
        core_resolved=last['core_cells']>=criteria['minimum_core_cells'],
        wall_resolved=last['wall_cells'] is not None and last['wall_cells']>=criteria['minimum_wall_cells'],
        mass_hierarchy=last['ma_ms']<=criteria['maximum_ma_ms'],
        ir_peak_sampled=last['kmin_tau_star']<=criteria['maximum_kmin_tau_star'],
        adiabatic_kmin=last['adiabatic_parameter_kmin']<=criteria['maximum_adiabatic_parameter'],
        small_angles=last['theta_max']<=.3,
        radial_near_vacuum=last['radial_rms_from_vacuum']<=.1,
        plateau_tested=False,ensemble_tested=False)
    report=dict(source=str(source.resolve()),parameters=cfg,contract=CONTRACT,
                snapshot_sha256=hashes,rows=rows,numerical_checks_passed=all(checks),
                comparison_checks={k:bool(v) for k,v in qualification.items()},
                qualified_for_paper_comparison=all(qualification.values()),
                post_defect_phase_advance=phase_advance,
                normalization=dict(a_star=a_star,h_star=h_star,f_a=1.,
                    control_note='For chi=0 this is the growing run fiducial normalization, not an actual mass-H equality'))
    (output/'analysis.json').write_text(json.dumps(report,indent=2)+'\n')
    write_csv(output/'history.csv',dict(
        tau=[r['tau'] for r in rows],tau_ratio=[r['tau_ratio'] for r in rows],
        a=[r['a'] for r in rows],chi=[r['chi'] for r in rows],
        plaquettes=[r['plaquettes'] for r in rows],wall_links=[r['wall_links'] for r in rows],
        nu_continuum=[r['continuum']['nu_diagnostic'] for r in rows],
        nu_lattice=[r['lattice']['nu_diagnostic'] for r in rows],
        n_temporal=[r['continuum']['temporal_estimator'] for r in rows],
        n_spatial=[r['continuum']['spatial_estimator'] for r in rows],
        rho_pq=[r['rho_pq_total'] for r in rows],
        rho_axion_harmonic_proxy=[r['continuum']['harmonic_axion_energy_proxy'] for r in rows]))
    print(source.name,'nu diagnostic',last['continuum']['nu_diagnostic'],
          'numerics',all(checks),'paper-qualified',report['qualified_for_paper_comparison'],flush=True)
    return report


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--source',type=Path,nargs='+',required=True)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'contract.json').write_text(json.dumps(CONTRACT,indent=2)+'\n')
    reports={}
    for i,source in enumerate(args.source):reports[f'run_{i:02d}_{source.name}']=audit(source,args.output/f'run_{i:02d}_{source.name}')
    (args.output/'summary.json').write_text(json.dumps(dict(
        sources=reports,python=platform.python_version(),numpy=np.__version__,
        implementation_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (ROOT/'checks/pq_strings/spectra.py',Path(__file__),ROOT/'checks/pq_strings/model.py')}),indent=2)+'\n')
    raise SystemExit(0 if all(r['numerical_checks_passed'] for r in reports.values()) else 1)
if __name__=='__main__':main()
