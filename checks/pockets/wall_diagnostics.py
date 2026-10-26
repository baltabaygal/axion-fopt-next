"""Measure local wall widths at every half-vacuum contour edge in a slice.

The original central-row gate is retained separately. Extra threshold
crossings are reported, and missing/ambiguous local walls are never dropped
from a pass decision. Slice normals omit out-of-plane gradients: these are
diagnostics, not proof of three-dimensional wall resolution.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import map_coordinates

from model import read_slice
from run_prototype import crossings, wall_cells


def sampled_gradient(phi,points,dx):
    """Periodic central gradients interpolated at points, without grid-sized gradients.

    Linear interpolation commutes with an integer lattice shift and difference.
    Only the arithmetic order differs from interpolating a full gradient grid.
    """
    components=[]
    for axis in range(phi.ndim):
        shift=np.zeros((phi.ndim,1));shift[axis,0]=1.
        plus=map_coordinates(phi,points.T+shift,order=1,mode='grid-wrap')
        minus=map_coordinates(phi,points.T-shift,order=1,mode='grid-wrap')
        components.append((plus-minus)/(2*dx))
    return np.stack(components,-1)


def local_widths(phi, dx, reach=8., samples_per_cell=8, details_limit=0, observer=None):
    """Normals point toward increasing phi; widths use nearest .1/.9 crossings.

    Coordinates are cell indices, not positions. All level-.5 edges are
    included, including disconnected contours. No spherical assumption.
    """
    points = []
    for axis in range(phi.ndim):
        following = np.roll(phi, -1, axis)
        mask = (phi < .5) != (following < .5)
        indices = np.array(np.nonzero(mask)).T
        if not len(indices):
            continue
        key = tuple(indices.T)
        interpolation = (.5-phi[key])/(following[key]-phi[key])
        q = indices.astype(float)
        q[:, axis] += interpolation
        points.append(q)
    if not points:
        return dict(contour_edges=0, measured=0, ambiguous=0, min_cells=None,
                    median_cells=None, fraction_below_four_cells=None,
                    diagnostic_pass=False,ambiguity_counts={})
    points = np.concatenate(points)
    normals = sampled_gradient(phi,points,dx)
    norms = np.linalg.norm(normals, axis=1)
    valid_gradient = norms > 1e-12
    normals = np.divide(normals, norms[:, None], out=np.zeros_like(normals), where=valid_gradient[:, None])
    half = int(np.ceil(reach/dx*samples_per_cell))
    offsets = np.linspace(-reach, reach, 2*half+1)
    widths = []
    ambiguous = 0
    details=[]
    reasons={}
    def record(point,normal,reason,values=None,left=None,right=None):
        reasons[reason]=reasons.get(reason,0)+1
        if len(details)<details_limit:
            details.append(dict(point_index=point.tolist(),normal=normal.tolist(),reason=reason,
                offsets=offsets.tolist() if values is not None else None,
                values=values.tolist() if values is not None else None,left=left,right=right))
    for point, normal, valid, norm in zip(points, normals, valid_gradient, norms):
        if not valid:
            ambiguous += 1
            record(point,normal,'vanishing_gradient')
            continue
        values = map_coordinates(phi, point[:, None]+normal[:, None]*offsets/dx,
                                 order=1, mode='grid-wrap')
        positions = []
        for level in (.1, .9):
            edges = np.flatnonzero((values[:-1] < level) & (values[1:] >= level))
            roots = offsets[edges]+(level-values[edges])/(values[edges+1]-values[edges])*np.diff(offsets)[edges]
            eligible = roots[roots <= 0] if level == .1 else roots[roots >= 0]
            positions.append((eligible[-1] if level == .1 else eligible[0]) if len(eligible) else None)
        left, right = positions
        if left is None or right is None or right <= left:
            ambiguous += 1
            record(point,normal,'missing_threshold_bracket',values,left,right)
            continue
        # More than one half-height crossing between the chosen endpoints
        # means a folded/multi-wall region, not a single measured wall.
        middle = (offsets >= left) & (offsets <= right)
        v = values[middle]
        half_crossings = np.count_nonzero((v[:-1] < .5) != (v[1:] < .5))
        if half_crossings != 1:
            ambiguous += 1
            record(point,normal,'multiple_half_crossings' if half_crossings>1 else 'no_sampled_half_crossing',values,left,right)
            continue
        widths.append(float((right-left)/dx))
        if observer is not None:
            observer(point,normal,float(norm),widths[-1])
    result=dict(contour_edges=len(points), measured=len(widths), ambiguous=ambiguous,
                ambiguity_counts=reasons,
                min_cells=min(widths) if widths else None,
                median_cells=float(np.median(widths)) if widths else None,
                fraction_below_four_cells=float(np.mean(np.array(widths) < 4)) if widths else None,
                diagnostic_pass=bool(widths) and ambiguous == 0 and min(widths) >= 4)
    if details_limit:
        result['ambiguous_details']=details
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True, help='Native output directory')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output must be fresh')
    results = []
    for path in sorted(args.source.glob('slice_*.bin')):
        meta, s = read_slice(path)
        dx = meta['box']/meta['n']; phi = s[..., 0]
        width = wall_cells(s, dx)
        result = dict(time=meta['u'], source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                      original_central_row_width=float(width) if np.isfinite(width) else None,
                      central_row_crossing_counts={str(v): len(crossings(phi[:,len(phi)//2], v, dx))
                                                   for v in (.1, .5, .9)},
                      local_normal=local_widths(phi, dx))
        results.append(result)
    if not results:
        raise ValueError('No slices')
    args.output.mkdir(parents=True, exist_ok=False)
    result = dict(source=str(args.source.resolve()), slices=results,
                  method_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  interpretation='Post-hoc diagnosis; does not replace original pilot gate or establish 3D convergence.')
    (args.output/'wall_audit.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    for row in results:
        print(row['time'], row['central_row_crossing_counts'], row['local_normal'])


if __name__ == '__main__':
    main()
