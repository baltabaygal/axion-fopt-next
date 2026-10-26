"""Choose detector geometry with continuum packet data, without solver runs."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from paper_reference import Parameters
from planar_initial import PlanarSeed,modes
from formation_scattering_reference import scatter


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    # Declare geometry and success requirements before evaluating either detector.
    contract=dict(n=2048,box=512.,time_end=220.,reference_output_interval=.5,
        original=dict(centre=96.,inner=80.,outer=112.),
        proposed=dict(centre=104.,inner=88.,outer=120.),
        gates=dict(initial_fraction_min=.999,final_fraction_min=.98,
                   minimum_during_encounter_max=.01,flux_balance_relative=2e-3))
    (args.output/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    seed=PlanarSeed();p=Parameters(dv=0.);n=contract['n'];box=seed.box;dx=box/n
    x=(np.arange(n)+.5)*dx;z=x-box/2
    indices,k,omega,weights=modes(seed,p)
    wave=np.empty((n,len(k)),dtype=complex);gradient=np.empty_like(wave)
    component_rows=[]
    for j,(kk,ww) in enumerate(zip(k,omega)):
        row=scatter(float(ww),extent=16.,samples=513)
        r=complex(row['reflection_real'],row['reflection_imag']);kr=np.sqrt(complex(ww**2-1.))
        left=z< -16.;right=z>16.;middle=~(left|right)
        wave[left,j]=np.exp(1j*kk*z[left])+r*np.exp(-1j*kk*z[left])
        gradient[left,j]=1j*kk*(np.exp(1j*kk*z[left])-r*np.exp(-1j*kk*z[left]))
        wave[right,j]=row['wave'][-1]*np.exp(1j*kr*(z[right]-16.));gradient[right,j]=1j*kr*wave[right,j]
        wave[middle,j]=np.interp(z[middle],row['x'],row['wave']);gradient[middle,j]=np.interp(z[middle],row['x'],row['derivative'])
        component_rows.append(dict(k=float(kk),omega=float(ww),r_real=r.real,r_imag=r.imag,
                                   flux_residual=row['flux_residual']))
    phi=.5*(1+np.tanh(z));mass2=p.mass_false**2*np.exp(p.kappa*phi)
    times=np.arange(0,contract['time_end']+.25,.5)
    energy0=float(box/2*seed.amplitude**2*np.sum(omega**2*weights**2))
    histories={name:np.empty((len(times),3)) for name in ('original','proposed')}
    window_data={}
    for name in histories:
        spec=contract[name];d=(x-spec['centre']+box/2)%box-box/2
        s=np.clip((abs(d)-spec['inner'])/(spec['outer']-spec['inner']),0,1)
        w=1-s**3*(10+s*(-15+6*s))
        wx=-30*s*s*(1-s)**2/(spec['outer']-spec['inner'])*np.sign(d)
        window_data[name]=(w,wx)
    phases=seed.amplitude*weights*np.exp(1j*k*(box/2-seed.centre))
    for start in range(0,len(times),8):
        t=times[start:start+8]
        coefficients=phases[:,None]*np.exp(-1j*omega[:,None]*t)
        a=(wave@coefficients).real;pa=(wave@(-1j*omega[:,None]*coefficients)).real;grad=(gradient@coefficients).real
        rho=.5*(pa*pa+grad*grad+mass2[:,None]*a*a);current=-pa*grad
        for name,(w,wx) in window_data.items():
            histories[name][start:start+len(t),:]=np.column_stack((t,np.sum(w[:,None]*rho,axis=0)*dx,
                                                               -np.sum(wx[:,None]*current,axis=0)*dx))
    rows={}
    for name,h in histories.items():
        np.savetxt(args.output/(name+'.csv'),h,delimiter=',',header='t,window_energy,outward_flux',comments='')
        flux=h[:,2];integral=np.r_[0.,np.cumsum(.5*(flux[1:]+flux[:-1])*np.diff(h[:,0]))]
        rows[name]=dict(initial_fraction=float(h[0,1]/energy0),final_fraction=float(h[-1,1]/energy0),
            minimum_during_encounter=float(h[:,1].min()/energy0),minimum_time=float(h[np.argmin(h[:,1]),0]),
            outward_crossing=float(np.trapezoid(np.maximum(flux,0.),h[:,0])/energy0),
            returned_crossing=float(np.trapezoid(np.maximum(-flux,0.),h[:,0])/energy0),
            flux_balance_relative=float(np.max(abs(h[:,1]-h[0,1]+integral))/energy0))
    r=rows['proposed'];g=contract['gates']
    checks=dict(initial=r['initial_fraction']>g['initial_fraction_min'],returned_coverage=r['final_fraction']>g['final_fraction_min'],
        incident_exit=r['minimum_during_encounter']<g['minimum_during_encounter_max'],
        energy_balance=r['flux_balance_relative']<g['flux_balance_relative'])
    result=dict(kind=__doc__,scope='Detector design from fixed-wall continuum ODE modes; not a new compiled result',
        contract=contract,reference_initial_energy=energy0,rows=rows,checks=checks,passed=all(checks.values()),
        components=component_rows,original_pilot_gate='Still failed; no old artifact or threshold changed',
        source_sha256={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in
            ('design_planar_detector.py','formation_scattering_reference.py','planar_initial.py','paper_reference.py')})
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(rows=rows,checks=checks,passed=result['passed']),indent=2))
    if not result['passed']:raise SystemExit(1)

if __name__=='__main__':main()
