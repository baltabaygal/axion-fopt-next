"""Unwrapped centroid of one localised pocket on a periodic cube."""
import numpy as np


def periodic_centroid(mask,box):
    n=len(mask)
    if mask.shape!=(n,n,n) or not mask.any():
        raise ValueError('Require a nonempty cubic pocket mask')
    x=(np.arange(n)+.5)*box/n;phase=np.exp(2j*np.pi*x/box)
    centre=[];coherence=[]
    for axis in range(3):
        weights=mask.sum(axis=tuple(a for a in range(3) if a!=axis));total=weights.sum()
        z=weights@phase;coherence.append(float(abs(z)/total))
        if abs(z)/total<1e-8:
            raise ValueError('Periodic pocket centroid is ambiguous')
        circular=(np.angle(z)*box/(2*np.pi))%box
        delta=(x-circular+box/2)%box-box/2
        c=(circular+float(weights@delta/total))%box
        distances=(x-c+box/2)%box-box/2
        if np.any((weights>0)&(np.abs(distances)>box/2-box/n/2)):
            raise ValueError('Pocket extends across its centroid unwrap cut')
        centre.append(c)
    return np.array(centre),coherence
