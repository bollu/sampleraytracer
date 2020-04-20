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

def planet(integrator, n, dt):
    STRENGTH = 0.5

    # minimise potential V(q): q, K(p, q) p^2
    q = np.array([0.0, 1.0])
    p = np.array([-1.0, 0.0])

    # H = STRENGTH * |q| (potential) + p^2/2 (kinetic)
    def H(qcur, pcur): return STRENGTH * np.linalg.norm(q) + np.dot(p, p) / 2
    def dhdp(qcur, pcur): return p
    def dhdq(qcur, pcur): return STRENGTH * 2 * q / np.linalg.norm(q)

    qs = []
    for i in range(n):
        print("q: %10s | p: %10s | H: %6.4f" % (q, p, H(q, p)))
        (q, p) = integrator(dhdp, dhdq, q, p, dt)
        qs.append(q.copy())
    return np.asarray(qs)

NITERS = 15

planet_euler = planet(euler, NITERS, dt=1)
plt.rcParams.update({'font.size': 12, 'font.family':'monospace'})
fig, ax = plt.subplots()
ax.plot(planet_euler[:, 0], planet_euler[:, 1], label='euler',
        linewidth=5, color='#D81B60')
legend = plt.legend(frameon=False)
ax.set_title("euler integrator: NITERS=%s dt=%s" % (NITERS, 1))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
plt.savefig("euler-dt-1.png")

planet_euler = planet(euler, NITERS*100, dt=1.0/100)
plt.rcParams.update({'font.size': 12, 'font.family':'monospace'})
fig, ax = plt.subplots()
ax.plot(planet_euler[:, 0], planet_euler[:, 1], label='euler',
        linewidth=5, color='#D81B60')
legend = plt.legend(frameon=False)
ax.set_title("euler integrator: NITERS=%s dt=%s" % (NITERS*100, 1.0/100))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
plt.savefig("euler-dt-1e-2.png")
plt.show()
