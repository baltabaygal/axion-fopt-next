"""Geometry diagnostics from saved symmetric second moments."""
import numpy as np


def asphericity(history):
    h=np.atleast_1d(history)
    covariance=np.empty((len(h),3,3))
    covariance[:,0,0]=h['sxx'];covariance[:,1,1]=h['syy'];covariance[:,2,2]=h['szz']
    covariance[:,0,1]=covariance[:,1,0]=h['sxy']
    covariance[:,0,2]=covariance[:,2,0]=h['sxz']
    covariance[:,1,2]=covariance[:,2,1]=h['syz']
    eigenvalues=np.linalg.eigvalsh(covariance)
    result=np.full(len(h),np.nan)
    valid=(h['fv_volume']>0)&(eigenvalues.mean(axis=1)>0)
    result[valid]=(eigenvalues[valid,-1]-eigenvalues[valid,0])/eigenvalues[valid].mean(axis=1)
    return result
