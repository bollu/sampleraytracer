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
def leapfrog(dhdq, dhdp, q0, p0, dt):
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



def hmc(q0, U, dU, nsteps, dt):
    M = 5.0
    def h(q, p): return U(q) + 0.5 / M  * p*p
    def nextsample(q, p):
        for _ in range(nsteps):
            def hdp(q, p): return 1.0 / M * p
            def hdq(q, p): return dU(q)
            (q, p) = leapfrog(hdq, hdp, q, p, dt)
        return (q, p)

    yield q0; q = q0
    while True:
        # number of simulation steps
        p = np.random.normal(0, M)
        (qnext, pnext) = nextsample(q, p)
        # print("(q:%4.2f, p:%4.2f) -> (q':%4.2f, p':%4.2f)" % (q, p, qnext, pnext))
        r = np.random.uniform();
        # if r < min(1, np.exp(h(q, p) - h(qnext, pnext))): q = qnext
        pnext = -p
        if np.log(r) < h(q, p) - h(qnext, pnext): q = qnext
        yield q

def take_every_nth(n, gen):
    for (i, x) in enumerate(gen):
        if i % n == 0: yield x

# HMC
def exp(x): return np.exp(-x*x)
def expgrad(x): return np.exp(-x*x) * -2*x
def expprop(x): return np.random.normal(loc=x, scale=1e-1)
def logexp(x): return -x*x
def logexpgrad(x): return -2*x

def neglogexp(x): return -1 * logexp(x)
def neglogexpgrad(x): return -1 * logexpgrad(x)

def negexp(x): return -1 * exp(x)
def negexpgrad(x): return -1 * expgrad(x)

COLORTRUTH = "#5C6BC0"; COLORSAMPLES = "#D81B60"
NSAMPLES = 1000; NSTEPS = 2; DT = 0.1; DECORRELATE_STEPS = 10
xs = list(take_every_nth(DECORRELATE_STEPS, itertools.islice(hmc(1, neglogexp, neglogexpgrad, NSTEPS, DT), NSAMPLES*DECORRELATE_STEPS)))
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
ax[0].set_title("HMC: #samples=%s | decorrelation steps: %s | dt: %4.2f | steps inside hmc: %s " % (len(xs), DECORRELATE_STEPS, DT, NSTEPS))
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[0].spines['bottom'].set_visible(False)
ax[0].spines['left'].set_visible(False)

ax[1].hist(xs, bins=200, cumulative=True, density=True, label='prob', linewidth=5, color=COLORSAMPLES)
ax[1].plot(fxs, fyscum, linewidth=5, color=COLORTRUTH, alpha=0.5)

ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['bottom'].set_visible(False)
ax[1].spines['left'].set_visible(False)

fig_size = plt.gcf().get_size_inches() #Get current size
plt.gcf().set_size_inches(2.0 * fig_size) 
plt.savefig("mcmc-hmc-1d-exp.png")
plt.show()


COLORTRUTH = "#5C6BC0"; COLORSAMPLES = "#D81B60"
NSAMPLES = 1000; NSTEPS = 20; DT = 1e-2; DECORRELATE_STEPS = 10
xs = list(take_every_nth(DECORRELATE_STEPS, itertools.islice(hmc(600, neglogexp, neglogexpgrad, NSTEPS, DT), NSAMPLES*DECORRELATE_STEPS)))
ys = [exp(x) for x in xs]
fxs = np.arange(np.min(xs)-1e-1, np.max(xs)+1e-1, (np.max(xs)+1e-1 - (np.min(xs) - 1e-1)) / 1000.0);
fys = [exp(x) for x in fxs];
fyscum = np.cumsum(fys); fyscum = fyscum / np.max(fyscum)
plt.rcParams.update({'font.size': 10, 'font.family':'monospace'})
fig, ax = plt.subplots(3, 1)
ax[0].plot(fxs, fys, label='prob',
        linewidth=5, color=COLORTRUTH, markersize=4.0, alpha=0.4)
ax[0].plot(xs, ys, 'x', label='prob',
        linewidth=5, color=COLORSAMPLES, markersize=4.0)

legend = plt.legend(frameon=False)
ax[0].set_title("HMC: #samples=%s | decorrelation steps: %s | dt: %4.2f | steps inside hmc: %s " % (len(xs), DECORRELATE_STEPS, DT, NSTEPS))
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[0].spines['bottom'].set_visible(False)
ax[0].spines['left'].set_visible(False)

ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['bottom'].set_visible(False)
ax[1].spines['left'].set_visible(False)

ax[1].set_title("HMC cumulative: #samples=%s " % (len(xs)))
ax[1].hist(xs, bins=400, cumulative=True, density=True, label='prob', linewidth=5, color=COLORSAMPLES)
ax[1].plot(fxs, fyscum, linewidth=5, color=COLORTRUTH, alpha=0.5)
fig_size = plt.gcf().get_size_inches() #get current size
plt.gcf().set_size_inches(2.0 * fig_size) 
plt.savefig("mcmc-hmc-1d-exp-startx-600.png")


RCUTOFF = 2
xszoom = [x for x in xs if x < RCUTOFF]
l = np.min(xs)-1e-1; r = min(RCUTOFF, np.max(xszoom)+1e-1)
fxszoom = np.arange(l, r, (r - l) / 1000.0);
fyszoom = [exp(x) for x in fxszoom];
fyszoomcum = np.cumsum(fyszoom); fyszoomcum = fyszoomcum / np.max(fyszoomcum)
ax[2].spines['top'].set_visible(False)
ax[2].spines['right'].set_visible(False)
ax[2].spines['bottom'].set_visible(False)
ax[2].spines['left'].set_visible(False)

ax[2].set_title("HMC: #samples(zoomed: %4.2f < x < %4.2f)=%s | " % (l, r, len(xszoom)))
ax[2].hist(xszoom, bins=400, cumulative=True, density=True, label='prob', linewidth=5, color=COLORSAMPLES)
ax[2].plot(fxszoom, fyszoomcum, linewidth=5, color=COLORTRUTH, alpha=0.5)
fig_size = plt.gcf().get_size_inches() #get current size
plt.gcf().set_size_inches(2.0 * fig_size) 
plt.savefig("mcmc-hmc-1d-exp-startx-600.png")
plt.show()
