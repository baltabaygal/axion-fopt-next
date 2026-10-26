"""Independent insertion map and shared Poisson-candidate preparation."""
from pathlib import Path
import sys
import numpy as np
from scipy.integrate import quad
from model import PARAMS
from expanding_model import scale,physical_terms

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from axion_fopt_next.core import sample_candidates

def rate(u,h0=.025,beta=.4):
    a=scale(u,h0);t=np.asarray(u)+h0*np.asarray(u)**2/2
    return h0**4*np.exp(beta*t)*a**4

def prepare_candidates(out,box=32.,end=12.,seed=0):
    out=Path(out);u=np.linspace(0,end,1025);q=rate(u)
    integrated=float(np.trapezoid(q,u));truth=float(quad(lambda x:float(rate(x)),0,end,epsabs=1e-13,epsrel=1e-12)[0])
    refined=np.linspace(0,end,2049);refined_integral=float(np.trapezoid(rate(refined),refined))
    candidate=sample_candidates(u,q,box,1.,seed,max_candidates=10000)
    np.savetxt(out/'rate.csv',np.column_stack((u,q)),delimiter=',',header='u,q_Gamma_a4',comments='')
    np.savetxt(out/'candidates.txt',np.column_stack((candidate.times,candidate.positions)),fmt='%.17g',header='u x y z; wall_speed metadata unused; not geometrically excluded')
    counts=np.array([sample_candidates(u,q,box,1.,s,max_candidates=10000).candidate_count for s in range(256)])
    mean=box**3*integrated;mean_se=np.sqrt(mean/len(counts));variance_se=mean*np.sqrt(2/(len(counts)-1))
    checks=dict(rate_quadrature=abs(integrated/truth-1)<.001,rate_refinement=abs(integrated/refined_integral-1)<.001,
        poisson_mean=abs(counts.mean()-mean)<5*mean_se,poisson_variance=abs(counts.var(ddof=1)-mean)<5*variance_se)
    return dict(checks={k:bool(v) for k,v in checks.items()},passed=bool(all(checks.values())),seed=seed,box=box,end=end,H0=.025,beta=.4,
        Gamma0=.025**4,rate_nodes=len(u),candidate_count=candidate.candidate_count,expected_count=mean,
        candidate_mean=float(counts.mean()),candidate_variance=float(counts.var(ddof=1)),ensemble=256,
        rate_quadrature_error=abs(integrated/truth-1),rate_refinement_error=abs(integrated/refined_integral-1),
        numpy=np.__version__,RNG='numpy.default_rng / PCG64; candidate file frozen, no RNG in native evolution')

def insertion(y,box,u,h0,centre,radius=8.,width=2.,tolerance=.05):
    n=len(y);dx=box/n;a=scale(u,h0);x=(np.arange(n)+.5)*dx
    distance=[]
    for c in centre:
        d=abs(x-c);distance.append(np.minimum(d,box-d))
    r=a*np.sqrt(distance[0][:,None,None]**2+distance[1][None,:,None]**2+distance[2][None,None,:]**2)
    core=r<=radius;maxcore=float(abs(y[...,0][core]).max()) if core.any() else 0.
    accepted=bool(core.any() and maxcore<=tolerance)
    if not accepted:return y.copy(),dict(accepted=False,maxcore=maxcore,cells=int(core.sum()))
    start=radius+2*width;stop=radius+4*width
    z=np.clip((r-start)/(stop-start),0,1);taper=1-10*z**3+15*z**4-6*z**5
    b=.5*(1-np.tanh((r-radius)/width))*taper;b[r>=stop]=0
    new=y.copy();new[...,0]=y[...,0]+(1-y[...,0])*b;new[...,2]=(1-b)*y[...,2]
    return new,dict(accepted=True,maxcore=maxcore,cells=int(core.sum()))

def energy_totals(y,box,a):
    terms=physical_terms(y,box/len(y),a);v={k:float(x.sum()*(box/len(y))**3) for k,x in terms.items()}
    scalar=v['K_S']+v['G_S']+v['U'];axion=v['K_A']+v['G_A']+v['V_A']
    return np.array([scalar+axion,scalar,axion])
