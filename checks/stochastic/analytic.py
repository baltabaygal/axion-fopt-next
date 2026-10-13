"""Independent survival references for common-speed bubbles, Rmax<L/2."""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import quad


@dataclass(frozen=True)
class Case:
    name: str
    box: float
    speed: float
    end: float
    constant_q: float | None = None
    h: float = 1.
    beta: float = 1.

    def rate(self, tau):
        tau = np.asarray(tau, dtype=float)
        if self.constant_q is not None:
            return np.full_like(tau, self.constant_q)
        return self.h**4 * np.exp(self.beta*(self.h*tau**2/2-1/(2*self.h))) * (self.h*tau)**4

    def cosmic_time(self, tau):
        return self.h*np.asarray(tau)**2/2 if self.constant_q is None else np.asarray(tau)


CASES = (Case("static",3.,.5,1.,constant_q=4.), Case("radiation",3.,.5,2.5))


def survival(case, times):
    prefactor=4*np.pi/3*case.speed**3
    values=[]
    for tau in times:
        integral=quad(lambda s: float(case.rate(s))*(tau-s)**3,0.,float(tau),epsabs=1e-11,epsrel=1e-11)[0]
        values.append(np.exp(-prefactor*integral))
    return np.array(values)


def cosmic_survival(case, times):
    """Independent cosmic-time evaluation of Gamma*a_birth^3*R_comoving^3."""
    if case.constant_q is not None:
        return np.exp(-np.pi*case.constant_q*case.speed**3*np.asarray(times)**4/3)
    values=[]
    for tau in times:
        t=case.h*tau**2/2
        def integrand(birth):
            a=np.sqrt(2*case.h*birth)
            radius=case.speed*(tau-np.sqrt(2*birth/case.h))
            gamma=case.h**4*np.exp(case.beta*(birth-1/(2*case.h)))
            return gamma*a**3*radius**3
        integral=quad(integrand,0.,t,epsabs=1e-11,epsrel=1e-11)[0]
        values.append(np.exp(-4*np.pi/3*integral))
    return np.array(values)


def tabulated_survival(nodes, q, times, speed):
    """Exact polynomial interval integrals for the sampler's linear rate."""
    answer=[]
    for tau in times:
        total=0.
        for i,left in enumerate(nodes[:-1]):
            if left>=tau: break
            width=min(nodes[i+1],tau)-left
            d=tau-left
            slope=(q[i+1]-q[i])/(nodes[i+1]-left)
            total += q[i]*(d**3*width-1.5*d**2*width**2+d*width**3-width**4/4)
            total += slope*(d**3*width**2/2-d**2*width**3+.75*d*width**4-width**5/5)
        answer.append(np.exp(-4*np.pi/3*speed**3*total))
    return np.array(answer)
