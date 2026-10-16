"""Periodic bubble-to-cell CSR index; exact paper-kernel distance predicate.

Numba is an experiment dependency, not a dependency of the core package.
The unchanged paper all-pairs kernel remains the fallback and oracle.
"""
import numpy as np
import numba as nb


@nb.njit(cache=True)
def build_index(centres,radii,box,ncells,max_references):
    counts=np.zeros(ncells**3,dtype=np.int64)
    lows=np.empty((len(centres),3),dtype=np.int64)
    spans=np.empty((len(centres),3),dtype=np.int64)
    # Bounds use float64 and a conservative float32-roundoff margin. The
    # final distance predicate, rather than cell overlap, decides conversion.
    pad=32*np.finfo(np.float32).eps*box
    total=0
    for b in range(len(centres)):
        for axis in range(3):
            c=float(centres[b,axis])+box/2
            low=int(np.floor((c-float(radii[b])-pad)*ncells/box))
            high=int(np.floor((c+float(radii[b])+pad)*ncells/box))
            lows[b,axis]=low
            spans[b,axis]=min(ncells,high-low+1)
        total+=spans[b,0]*spans[b,1]*spans[b,2]
        if total>max_references:
            return np.empty(0,dtype=np.int64),np.empty(0,dtype=np.int32)
        for i in range(spans[b,0]):
            x=(lows[b,0]+i)%ncells
            for j in range(spans[b,1]):
                y=(lows[b,1]+j)%ncells
                for k in range(spans[b,2]):
                    z=(lows[b,2]+k)%ncells
                    counts[(x*ncells+y)*ncells+z]+=1
    offsets=np.empty(len(counts)+1,dtype=np.int64);offsets[0]=0
    for i in range(len(counts)):offsets[i+1]=offsets[i]+counts[i]
    entries=np.empty(total,dtype=np.int32)
    cursor=offsets[:-1].copy()
    for b in range(len(centres)):
        for i in range(spans[b,0]):
            x=(lows[b,0]+i)%ncells
            for j in range(spans[b,1]):
                y=(lows[b,1]+j)%ncells
                for k in range(spans[b,2]):
                    z=(lows[b,2]+k)%ncells
                    cell=(x*ncells+y)*ncells+z
                    entries[cursor[cell]]=b;cursor[cell]+=1
    return offsets,entries


@nb.njit(parallel=True,fastmath=True,cache=True)
def query_index(tracers,active,centres,radii,box,ncells,offsets,entries):
    survived=np.ones(len(active),dtype=nb.boolean)
    for idx in nb.prange(len(active)):
        pos=tracers[active[idx]]
        x=int(np.floor((float(pos[0])+box/2)*ncells/box))%ncells
        y=int(np.floor((float(pos[1])+box/2)*ncells/box))%ncells
        z=int(np.floor((float(pos[2])+box/2)*ncells/box))%ncells
        cell=(x*ncells+y)*ncells+z
        for entry in range(offsets[cell],offsets[cell+1]):
            b=entries[entry]
            # Preserve float32 subtraction and the original nearest-image
            # operation/order. CSR entries retain original bubble order.
            dx=pos[0]-centres[b,0];dy=pos[1]-centres[b,1];dz=pos[2]-centres[b,2]
            dx-=box*round(dx/box);dy-=box*round(dy/box);dz-=box*round(dz/box)
            radius_sq=radii[b]**2
            if dx*dx+dy*dy+dz*dz<=radius_sq:
                survived[idx]=False;break
    return active[survived]


def accelerated_paper(paper,ncells=16,max_references=8_000_000):
    """Substitute only engulfment; preserve random draws and birth ordering."""
    base=paper['NucleationTracerManager']
    class IndexedTracer(base):
        def update_state(self,centres,radii):
            if self.num_active_tracers==0 or len(centres)==0:return
            offsets,entries=build_index(centres,radii,self.box_size,ncells,max_references)
            if len(offsets)==0:
                super().update_state(centres,radii)
            else:
                self.active_indices=query_index(self.tracers,self.active_indices,centres,radii,
                                                self.box_size,ncells,offsets,entries)
    return dict(paper,NucleationTracerManager=IndexedTracer)
