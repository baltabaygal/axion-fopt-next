"""Necessary scale-hierarchy bounds for uniform fixed-resolution PQ production runs.

This audit launches no simulations. It exposes assumptions rather than choosing
parameters to make a desired yield appear. Core growth used only in preparation
cannot change these physical-production bounds.
"""
import argparse
import csv
import json
import math
from pathlib import Path


def bounds(ratio,exponent=8.,core_cells=4.,wall_cells=4.,mass_ratio=.16,
           peak=3.,modes_to_peak=2.,horizon_factor=2.):
    if ratio<=0 or min(core_cells,wall_cells,mass_ratio,peak,modes_to_peak,horizon_factor)<=0 or exponent<0:
        raise ValueError('Require positive scale criteria and nonnegative exponent')
    # a*m_a*tau_star=r^(1+n/2), independent of the a=1 normalization.
    power=ratio**(1+exponent/2)
    B_ir=2*math.pi*modes_to_peak/peak
    B_horizon=horizon_factor*ratio
    B=max(B_ir,B_horizon)
    required_core=core_cells*B*power/mass_ratio
    required_wall=wall_cells*B*power
    required=math.ceil(max(required_core,required_wall))
    return dict(tau_ratio=ratio,n_tau=exponent,minimum_L_over_tau_star=B,
                ir_box_minimum=B_ir,horizon_box_minimum=B_horizon,
                minimum_N=required,core_bound_N=required_core,wall_bound_N=required_wall,
                four_component_state_bytes=32*required**3,
                assumptions=dict(core_cells=core_cells,wall_cells=wall_cells,
                    max_ma_ms=mass_ratio,peak_k_tau_star=peak,
                    fundamental_intervals_to_peak=modes_to_peak,
                    L_over_tau_end_minimum=horizon_factor))


def best_ratio(n,**criteria):
    # Monotone necessary inequalities: find largest allowed r at optimal B.
    low,high=1e-6,10.
    for _ in range(80):
        mid=(low+high)/2
        if bounds(mid,**criteria)['minimum_N']<=n:low=mid
        else:high=mid
    return low


def evaluate(n,box,origin,star,lam,exponent,ratio):
    tau=ratio*star;a=tau/origin;dx=box/n
    ma=origin/star**2*ratio**(exponent/2);ms=math.sqrt(2*lam)
    return dict(N=n,tau_ratio=ratio,tau=tau,a=a,ma=ma,ms=ms,ma_ms=ma/ms,
                core_cells=1/(a*ms*dx),wall_cells=1/(a*ma*dx),
                kmin_tau_star=2*math.pi*star/box,L_over_tau_end=box/tau)


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    contract=dict(target_tau_ratio=3.,n_tau=8.,core_cells=4.,wall_cells=4.,mass_ratio=.16,
                  peak=3.,modes_to_peak=2.,horizon_factor=2.,
                  scope='necessary lower bounds; not sufficient conditions for scaling or accurate yield',
                  preparation_note='lambda may vary during preparation but remains constant in physical production',
                  memory_note='four double component arrays only; excludes ghosts, work arrays and FFTs')
    (args.output/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    target=bounds(3.)
    reference_relaxed=bounds(3.,core_cells=1.,mass_ratio=.2)
    grid_table=[dict(N=n,best_possible_tau_ratio=best_ratio(n),
                     target_feasible=n>=target['minimum_N'],state_bytes=32*n**3)
                for n in (32,64,128,256,512,1024,2048,4096,8192,12288,32768,65536)]
    index_table=[bounds(3.,exponent=n) for n in (0.,2.,4.,6.,7.,8.)]
    current=[evaluate(64,24.,20.,27.,.08,8.,r) for r in (32/27,1.5,2.,2.5,3.)]
    result=dict(contract=contract,target=target,one_cell_core_relaxed=reference_relaxed,
                grids=grid_table,mass_growth_indices=index_table,current_setup=current,
                conclusion='No 32^3,64^3 or128^3 fixed-comoving full-PQ run meets these post-collapse criteria; changing units, lambda or tau_star cannot remove the combined bound.')
    (args.output/'feasibility.json').write_text(json.dumps(result,indent=2)+'\n')
    for filename,rows in [('grids.csv',grid_table),('current_setup.csv',current)]:
        with (args.output/filename).open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
