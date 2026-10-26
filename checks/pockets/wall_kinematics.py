"""Contour velocity versus measured width; no resolution-gate changes.

The level-set speed -pi_phi/|grad phi| is not necessarily a material speed.
Values >=1 are retained as diagnostic failures of a simple boosted-kink
interpretation, not clipped into a Lorentz factor.
"""
import numpy as np
from scipy.ndimage import map_coordinates
from wall_diagnostics import local_widths


def measure(phi,pi,dx,reach=8.):
    records=[]
    def observe(point,normal,gradient,width):
        momentum=float(map_coordinates(pi,point[:,None],order=1,mode='grid-wrap')[0])
        records.append([width,-momentum/gradient,*point])
    wall=local_widths(phi,dx,reach=reach,observer=observe)
    data=np.array(records,dtype=float).reshape((-1,2+phi.ndim))
    return wall,data


def summary(data,dx):
    if not len(data):return dict(measured=0)
    width,speed=data[:,0],np.abs(data[:,1]);thin=width<4;valid=speed<1
    def group(mask):
        selected=speed[mask]
        return dict(count=int(mask.sum()),median_abs_level_speed=float(np.median(selected)) if len(selected) else None,
                    fraction_abs_speed_ge_one=float(np.mean(selected>=1)) if len(selected) else None)
    return dict(measured=len(data),all=group(np.ones(len(data),bool)),thin=group(thin),resolved=group(~thin),
                minimum_width_cells=float(width.min()),minimum_width_abs_level_speed=float(speed[np.argmin(width)]),
                median_boost_corrected_width=float(np.median(width[valid]*dx/np.sqrt(1-speed[valid]**2))) if valid.any() else None,
                interpretation='Boost-corrected width is conditional on a local planar material-wall interpretation; no resolution or stability certificate.')
