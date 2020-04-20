#!/usr/bin/env python3
# Run HMC with a particular choice of potential
import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy.linalg
import argparse



## dq/dt = dH/dp|_{p0, q0}
## dp/dt = -dH/dq|_{p0, q0}
def leapfrog(dhdp, dhdq, q0, p0, dt):
    p0 += -dhdq(q0, p0) * 0.5 * dt

    # full step position
    # q += dt * p
    q0 += dhdp(q0, p0) * dt

    # half step position
    p0 += -dhdq(q0, p0) * 0.5 * dt
    return (q0, p0)

def euler(dhdp, dhdq, q, p, dt):
   pnew = p + -dhdq(q, p) * dt
   qnew = q + dhdp(q, p) * dt
   return (qnew, pnew)

def startpos():
    q = np.array([0.0, 1.0]); p = np.array([-1.0, 0.0]); return (q, p)

def planet(q, p, integrator, n, dt):
    STRENGTH = 0.5

    # minimise potential V(q): q, K(p, q) p^2

    # H = STRENGTH * |q| (potential) + p^2/2 (kinetic)
    def H(qcur, pcur): return STRENGTH * np.linalg.norm(q) + np.dot(p, p) / 2
    def dhdp(qcur, pcur): return p
    def dhdq(qcur, pcur): return STRENGTH * 2 * q / np.linalg.norm(q)

    qs = [q]
    ps = [p]
    for i in range(n):
        print("q: %10s | p: %10s | H: %6.4f" % (q, p, H(q, p)))
        (q, p) = integrator(dhdp, dhdq, q, p, dt)
        qs.append(q.copy())
        ps.append(p.copy())
    return np.asarray(qs), np.asarray(ps)

def plot(NITERS, DT, integrator, integrator_name):
    qs, ps = planet(*startpos(), integrator, NITERS, DT)
    plt.rcParams.update({'font.size': 20, 'font.family':'monospace'})
    fig, ax = plt.subplots()
    ax.plot(qs[:, 0], qs[:, 1], 'x', label='%s (fwd)' % (integrator_name, ),
            linewidth=5, color='#D81B60', markersize=10.0)

    qs, _ = planet(qs[-1], -ps[-1], integrator, NITERS, DT)
    ax.plot(qs[:, 0], qs[:, 1], '+', label='%s (bwd)' % (integrator_name, ),
            linewidth=5, color='#5C6BC0', markersize=10.0)

    legend = plt.legend(frameon=False)
    ax.set_title("%s integrator: NITERS=%s dt=%s" % (integrator_name, NITERS, DT))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.savefig("euler-dt-%s.png" % DT)
    plt.show()
plot(400, 1e-1, euler, 'euler')
# plot(30*100, 1e-2, euler, 'euler')

plot(400, 1e-1, leapfrog, 'leapfrog')
# plot(30*100, 1e-2, leapfrog, 'leapfrog')
