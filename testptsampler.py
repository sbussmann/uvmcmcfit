import numpy as np
from emcee import PTSampler

# mu1 = [1, 1], mu2 = [-1, -1]
mu1 = np.ones(2)
mu2 = -np.ones(2)

# Width of 0.1 in each dimension
sigma1inv = np.diag([100.0, 100.0])
sigma2inv = np.diag([100.0, 100.0])

def logl(x):
    dx1 = x - mu1
    dx2 = x - mu2
    loglike = np.logaddexp(-np.dot(dx1, np.dot(sigma1inv, dx1))/2.0,
                        -np.dot(dx2, np.dot(sigma2inv, dx2))/2.0)

    return loglike

# Use a flat prior
def logp(x):
    return 0.0

ntemps = 20
nwalkers = 100
ndim = 2

sampler=PTSampler(ntemps, nwalkers, ndim, logl, logp)

p0 = np.random.uniform(low=-1.0, high=1.0, size=(ntemps, nwalkers, ndim))
pexamine = []
for p, lnprob, lnlike in sampler.sample(p0, iterations=1000):
    #print(lnlike[0,0])
    pexamine.extend(lnlike)
    continue

# re-arrange the positions into something more usable
#temperature = 10
#pnew = []
#for i in range(1000):
#    pnew.append(pexamine[i][temperature, 0, 0])
#
#pnew = np.array(pnew).flatten()
from IPython import embed
embed()
import ipdb; ipdb.set_trace()
sampler.reset()

for p, lnprob, lnlike in sampler.sample(p, lnprob0=lnprob,
                                           lnlike0=lnlike,
                                           iterations=5000, thin=10):
    pass

assert sampler.chain.shape == (ntemps, nwalkers, 500, ndim)

# Chain has shape (ntemps, nwalkers, nsteps, ndim)
# Zero temperature mean:
mu0 = np.mean(np.mean(sampler.chain[0,...], axis=0), axis=0)

# Longest autocorrelation length (over any temperature)
max_acl = np.max(sampler.acor)

import ipdb; ipdb.set_trace()
