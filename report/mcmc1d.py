#!/usr/bin/env python3
# Run HMC with a particular choice of potential
import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy.linalg
import argparse
import itertools

np.random.seed(0)



## dq/dt = dH/dp|_{p0, q0}
## dp/dt = -dH/dq|_{p0, q0}
def leapfrog(dhdp, dhdq, q0, p0, dt):
    # kick: half step momentum
    p0 += -dhdq(q0, p0) * 0.5 * dt

    # drift: full step position
    q0 += dhdp(q0, p0) * dt

    # kick: half step momentum
    p0 += -dhdq(q0, p0) * 0.5 * dt
    return (q0, p0)

def euler(dhdp, dhdq, q, p, dt):
   pnew = p + -dhdq(q, p) * dt
   qnew = q + dhdp(q, p) * dt
   return (qnew, pnew)

def mhsimple(x0, prob, prop):
    yield x0; x = x0;
    while True:
        xnext = prop(x); p = prob(x); pnext = prob(xnext)
        r = np.random.uniform() + 1e-5;
        if r < pnext/p: x = xnext
        yield xnext

def mh_uncorr(x0, prob, prop, iters_per_sample):
    yield x0; x = x0;
    while True:
        for i in range(iters_per_sample):
            xnext = prop(x); p = prob(x); pnext = prob(xnext)
            r = np.random.uniform() + 1e-5;
            if np.log(r) < min(0, np.log(pnext) - np.log(p)): x = xnext
        yield xnext

# MH
def exp(x): return np.exp(-x*x)
def expgrad(x): return -2*x
def expprop(x): return np.random.normal(loc=x, scale=1e-1)

### mhsimple ###
NSAMPLES = 1000
COLORTRUTH = "#5C6BC0"; COLORSAMPLES = "#D81B60"
xs = list(itertools.islice(mhsimple(0, exp, expprop), NSAMPLES))
ys = [exp(x) for x in xs] 
fxs = np.arange(np.min(xs)-1e-1, np.max(xs)+1e-1, (np.max(xs)+1e-1 - (np.min(xs) - 1e-1)) / 100.0);
fys = [exp(x) for x in fxs];
fyscum = np.cumsum(fys); fyscum = fyscum / np.max(fyscum)
plt.rcParams.update({'font.size': 10, 'font.family':'monospace'})
fig, ax = plt.subplots(2, 1)
ax[0].plot(fxs, fys, label='prob',
        linewidth=5, color=COLORTRUTH, markersize=4.0, alpha=0.4)
ax[0].plot(xs, ys, 'x', label='prob',
        linewidth=5, color=COLORSAMPLES, markersize=4.0)

legend = plt.legend(frameon=False)
ax[0].set_title("metropolis-hastings-simple simple #samples=%s"% NSAMPLES)
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[0].spines['bottom'].set_visible(False)
ax[0].spines['left'].set_visible(False)

ax[1].set_title("metropolis-hastings-simple (cumulative sum) #samples=%s"% NSAMPLES)
ax[1].hist(xs, bins=80, cumulative=True, density=True, label='prob', linewidth=5, color=COLORSAMPLES)
ax[1].plot(fxs, fyscum, linewidth=5, color=COLORTRUTH)
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['bottom'].set_visible(False)
ax[1].spines['left'].set_visible(False)
fig_size = plt.gcf().get_size_inches() #Get current size
plt.gcf().set_size_inches(2.0 * fig_size) 
plt.savefig("mcmc-mh-simple-1d-exp.png")
plt.show()

## MH UNCORR ##

NSAMPLES = 1000; ITERS_PER_SAMPLE = 20
COLORTRUTH = "#5C6BC0"; COLORSAMPLES = "#D81B60"
xs = list(itertools.islice(mh_uncorr(0, exp, expprop, ITERS_PER_SAMPLE), NSAMPLES))
ys = [exp(x) for x in xs] 
fxs = np.arange(np.min(xs)-1e-1, np.max(xs)+1e-1, (np.max(xs)+1e-1 - (np.min(xs) - 1e-1)) / 100.0);
fys = [exp(x) for x in fxs];
fyscum = np.cumsum(fys); fyscum = fyscum / np.max(fyscum)
plt.rcParams.update({'font.size': 10, 'font.family':'monospace'})
fig, ax = plt.subplots(2, 1)
ax[0].plot(fxs, fys, label='prob',
        linewidth=5, color=COLORTRUTH, markersize=4.0, alpha=0.4)
ax[0].plot(xs, ys, 'x', label='prob',
        linewidth=5, color=COLORSAMPLES, markersize=4.0)

legend = plt.legend(frameon=False)
ax[0].set_title("metropolis-hastings-uncorrelate simple #samples=%s | #iters/sample: %s"% (NSAMPLES, ITERS_PER_SAMPLE))
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[0].spines['bottom'].set_visible(False)
ax[0].spines['left'].set_visible(False)

ax[1].set_title("metropolis-hastings-uncorrelate simple (cumulative sum) #samples=%s | #iters/sample: %s"% (NSAMPLES, ITERS_PER_SAMPLE))
ax[1].hist(xs, bins=80, cumulative=True, density=True, label='prob', linewidth=5, color=COLORSAMPLES)
ax[1].plot(fxs, fyscum, linewidth=5, color=COLORTRUTH)
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['bottom'].set_visible(False)
ax[1].spines['left'].set_visible(False)
fig_size = plt.gcf().get_size_inches() #Get current size
plt.gcf().set_size_inches(2.0 * fig_size) 
plt.savefig("mcmc-mh-uncorr-1d-exp.png")
plt.show()


def hmc(q0, logprob, logprobgrad, iters_per_sample, dt):
    yield q0; q = q0
    while True:
        # number of simulation steps
        p = np.random.normal(0, 1)
        for i in range(iters_per_sample):
            def hdp(p, q): return 0.5*p*p
            def hdq(p, q): return logprobgrad(q)
            (p, q) = leapfrog(hdp, hdq, p, q, dt)
        yield q

# HMC
NSAMPLES = 1000
ITERS_PER_SAMPLE = 10
DT = 1e-2
xs = list(itertools.islice(hmc(3, exp, expgrad, ITERS_PER_SAMPLE, DT), NSAMPLES))
ys = [exp(x) for x in xs]
plt.rcParams.update({'font.size': 10, 'font.family':'monospace'})
fig, ax = plt.subplots(2, 1)
ax[0].plot(xs, ys, 'x', label='prob',
        linewidth=5, color=COLORSAMPLES, markersize=4.0)

legend = plt.legend(frameon=False)
ax[0].set_title("metropolis-hastings: #samples=%s | #iters per sample=%s" % (NSAMPLES, ITERS_PER_SAMPLE))
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[0].spines['bottom'].set_visible(False)
ax[0].spines['left'].set_visible(False)

fxs = np.arange(-3, 3, (3 - (-3)) / 30.0); fys = np.cumsum([exp(x) for x in fxs]); fys = fys / np.max(fys)
ax[1].hist(xs, bins=80, cumulative=True, density=True, label='prob', linewidth=5, color=COLORTRUTH)
ax[1].plot(fxs, fys)
plt.show()
