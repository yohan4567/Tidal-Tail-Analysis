"""Globally bounded one-dimensional profile for a mixture of Poisson templates."""
import heapq
import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp,xlogy


def fixed_rate_fit(counts,signals,background,tolerance=1e-8):
    n=np.asarray(counts,float).ravel(); b=np.asarray(background,float).ravel()
    s=np.atleast_2d(np.asarray(signals,float))
    if s.ndim!=2 or s.shape[1]!=n.size or b.shape!=n.shape or not len(s):
        raise ValueError('Mismatched template shapes')
    if any(not np.isfinite(a).all() or (a<0).any() for a in (n,b,s)) or not np.equal(n,np.round(n)).all() or b.sum()<=0:
        raise ValueError('Finite nonnegative templates, integer counts and positive background required')
    if not np.isfinite(tolerance) or tolerance<=0: raise ValueError('Positive profile tolerance required')
    upper=max(float(n.sum()/b.sum()),1e-12)
    valid=~np.any((s==0)&(b[None,:]==0)&(n[None,:]>0),axis=1)
    if not valid.any(): raise ValueError('Templates cannot explain observed counts')
    total=len(s); s=s[valid]
    def component(beta):
        lam=s+np.asarray(beta)[:,None]*b if np.ndim(beta) else s+beta*b
        return np.sum(xlogy(n,lam)-lam,axis=1)
    modes=[]
    for v in s:
        def derivative(beta):
            lam=v+beta*b
            if np.any((lam==0)&(n>0)&(b>0)): return np.inf
            ratio=np.divide(n,lam,out=np.zeros_like(lam),where=lam>0)
            return float(np.sum(b*(ratio-1)))
        if derivative(0)<=0: modes.append(0.)
        elif derivative(upper)>=0: modes.append(upper)
        else: modes.append(brentq(derivative,0,upper,xtol=1e-12))
    modes=np.array(modes)
    def likelihood(beta): return float(logsumexp(component(beta))-np.log(total))
    def bound(lo,hi):
        # Each component is concave in beta. Sum of its separate interval maxima
        # is an upper bound on the mixture likelihood, even for multiple modes.
        separate=float(logsumexp(component(np.clip(modes,lo,hi)))-np.log(total))
        mid=(lo+hi)/2; lam=s+mid*b
        slopes=np.sum(b*(np.divide(n,lam,out=np.zeros_like(lam),where=lam>0)-1),axis=1)
        centre=component(mid)
        # Concavity bounds each component by its tangent. The sum of exponentials
        # of affine tangents is convex, so its maximum is at an interval endpoint.
        tangent=max(float(logsumexp(centre+slopes*(lo-mid))),float(logsumexp(centre+slopes*(hi-mid))))-np.log(total)
        return min(separate,tangent)
    candidates=np.unique(np.r_[0.,upper,modes])
    values=np.array([likelihood(v) for v in candidates]); index=int(np.argmax(values))
    best=float(values[index]); beta=float(candidates[index]); queue=[]
    for lo,hi in zip(candidates[:-1],candidates[1:]):
        heapq.heappush(queue,(-bound(lo,hi),float(lo),float(hi)))
    iterations=0
    while queue:
        neg,lo,hi=heapq.heappop(queue)
        if -neg<=best+tolerance: continue
        mid=(lo+hi)/2
        value=likelihood(mid)
        if value>best: best=value; beta=mid
        if mid==lo or mid==hi: raise RuntimeError('Profile search exhausted floating-point precision')
        for left,right in ((lo,mid),(mid,hi)):
            ceiling=bound(left,right)
            if ceiling>best+tolerance: heapq.heappush(queue,(-ceiling,left,right))
        iterations+=1
        if iterations>100000: raise RuntimeError('Could not certify background profile tolerance')
    return beta,-best
