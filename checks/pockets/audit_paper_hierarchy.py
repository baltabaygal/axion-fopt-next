"""Reproducible model/resource audit before new formation experiments."""
import argparse
import hashlib
import json
from pathlib import Path
import math
from paper_reference import head_on_reflection


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    source=Path('outputs/pockets/global-stiffness-01/n384_state.json')
    state=json.loads(source.read_text());spec=state['spec'];lam=state['lambda_'];w=spec['width']
    sigma=math.sqrt(lam/2)/6;rc=3*sigma/spec['delta_v'];thickness=2*math.atanh(.8)*w
    rows=[]
    for speed in (0.,.6,.9,.999):
        kinematics=head_on_reflection(spec['omega_max'],0.,speed)
        rows.append(dict(speed=speed,**kinematics,
            initial_thickness_cells=thickness/(spec['box']/spec['n'])/kinematics['gamma'],
            necessary_head_on_reflection=kinematics['wall_frame_frequency']<spec['mass_tv']))
    target=[]
    for rc_over_w in (10.,80.,400.):
        # Assumed thin-wall hierarchy and instantaneous boosted tanh profile.
        radius_over_w=rc_over_w/.01;gamma=1/math.sqrt(1-.999**2)
        minimum_n=math.ceil(3.2*radius_over_w*4*gamma/(2*math.atanh(.8)))
        target.append(dict(assumed_rc_over_w=rc_over_w,quartic_bias_over_degenerate_barrier=8/rc_over_w,
            assumed_rc_over_radius=.01,initial_speed=.999,box_over_radius=3.2,
            minimum_n_for_four_wall_cells=minimum_n,
            four_field_storage_gib=32*float(minimum_n)**3/2**30))
    # Independently substitute pressure balance into the full energy budget.
    dv=spec['delta_v'];R=spec['r0'];volume=4*math.pi/3*R**3;area=4*math.pi*R**2
    pressure=dv+2*sigma/R;direct=sigma*area+(dv+3*pressure)*volume
    derived=dv*volume*(4+3*rc/R);printed=dv*volume*(4+rc/R)
    report=dict(sources={str(source):hashlib.sha256(source.read_bytes()).hexdigest(),
        'paper_reference.py':hashlib.sha256(Path(__file__).with_name('paper_reference.py').read_bytes()).hexdigest(),
        'audit_paper_hierarchy.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        paper='https://arxiv.org/html/2407.01676v2',
        old_toy=dict(rc=rc,rc_over_radius=rc/R,rc_over_width=rc/w,phi_mass_coupling='phi^4, not exponential',initial_wall_speed=0.),
        moving_wall_checks=rows,uniform_grid_estimates=target,
        finite_tension_algebra=dict(full_budget=direct,pressure_substitution=derived,
             printed_pdf_eq_5_2=printed,substitution_relative_error=abs(direct-derived)/direct,
             printed_equality_relative_difference=abs(direct-printed)/direct,
             conclusion='Printed finite-tension equality is inconsistent with its unsimplified budget and Eq. (5.1); leading bulk limit unaffected. Use full budget, not a silently corrected paper formula.'),
        interpretation='Planning/reference audit only. Assumed tanh boosts and thin-wall ratios are not measured field dynamics; no new formation or stability result.')
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
