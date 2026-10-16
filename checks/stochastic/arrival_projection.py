"""Exact periodic first-arrival projection with conservative cell bounds."""
import numpy as np
import numba as nb


@nb.njit(parallel=True,cache=True)
def project(points,times,centres,box,speed,cells,starts,order,minbirth):
    result=np.empty(len(points));width=box/cells
    for p in nb.prange(len(points)):
        best=np.inf
        x,y,z=points[p]
        near=(int(x/width)*cells+int(y/width))*cells+int(z/width)
        # Seed a finite upper bound with bubbles in the point's own cell.
        for at in range(starts[near],starts[near+1]):
            b=order[at]
            dx=abs(x-centres[b,0]);dy=abs(y-centres[b,1]);dz=abs(z-centres[b,2])
            dx=min(dx,box-dx);dy=min(dy,box-dy);dz=min(dz,box-dz)
            best=min(best,times[b]+np.sqrt(dx*dx+dy*dy+dz*dz)/speed)
        for cell in range(cells**3):
            if starts[cell]==starts[cell+1] or cell==near:continue
            iz=cell%cells;iy=(cell//cells)%cells;ix=cell//(cells*cells)
            dx=abs(x-(ix+.5)*width);dy=abs(y-(iy+.5)*width);dz=abs(z-(iz+.5)*width)
            dx=max(0.,min(dx,box-dx)-width/2);dy=max(0.,min(dy,box-dy)-width/2);dz=max(0.,min(dz,box-dz)-width/2)
            bound=minbirth[cell]+np.sqrt(dx*dx+dy*dy+dz*dz)/speed
            # Conservatively pad the lower bound to protect rounding at faces.
            if bound>best+1e-12:continue
            for at in range(starts[cell],starts[cell+1]):
                b=order[at]
                if times[b]>=best:break
                dx=abs(x-centres[b,0]);dy=abs(y-centres[b,1]);dz=abs(z-centres[b,2])
                dx=min(dx,box-dx);dy=min(dy,box-dy);dz=min(dz,box-dz)
                best=min(best,times[b]+np.sqrt(dx*dx+dy*dy+dz*dz)/speed)
        result[p]=best
    return result


def arrival_map(catalogue,n,cells=8):
    centres=np.ascontiguousarray(catalogue.positions,dtype=np.float64)
    times=np.ascontiguousarray(catalogue.times,dtype=np.float64)
    labels=np.minimum((centres/catalogue.box_size*cells).astype(np.int64),cells-1)
    ids=(labels[:,0]*cells+labels[:,1])*cells+labels[:,2]
    # Stable ordering preserves increasing birth times inside each cell.
    order=np.argsort(ids,kind='stable').astype(np.int64)
    starts=np.r_[0,np.cumsum(np.bincount(ids,minlength=cells**3))].astype(np.int64)
    first=np.full(cells**3,np.inf)
    for cell in range(cells**3):
        if starts[cell]<starts[cell+1]:first[cell]=times[order[starts[cell]]]
    x=(np.arange(n)+.5)*catalogue.box_size/n
    points=np.stack(np.meshgrid(x,x,x,indexing='ij'),axis=-1).reshape(-1,3)
    return project(points,times,centres,catalogue.box_size,catalogue.wall_speed,cells,starts,order,first).reshape(n,n,n)
