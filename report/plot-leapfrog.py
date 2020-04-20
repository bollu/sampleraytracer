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
TIMESTEP = 1

print("planet simulation with leapfrog")
planet_leapfrog = planet(leapfrog, NITERS, TIMESTEP)

plt.rcParams.update({'font.size': 12, 'font.family':'monospace'})
fig, ax = plt.subplots()
print(planet_leapfrog)
ax.plot(planet_leapfrog[:, 0], planet_leapfrog[:, 1], label='leapfrog',
        linewidth=5, color='#00ACC1')
print("planet simulation with euler")
planet_euler = planet(euler, NITERS, TIMESTEP)
ax.plot(planet_euler[:, 0], planet_euler[:, 1], label='euler',
        linewidth=5, color='#D81B60')

legend = plt.legend(frameon=False)
ax.set_title("leapfrog v/s euler: NITERS=%s dt=%s" % (NITERS, TIMESTEP))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
plt.savefig("leapfrog-dt-1.png")


NITERS = 15 * 200
TIMESTEP = 1 * 0.01
fig, ax = plt.subplots()
print("planet simulation with euler")
planet_euler = planet(euler, NITERS, TIMESTEP)
ax.plot(planet_euler[:, 0], planet_euler[:, 1], label='euler',
        linewidth=5, color='#D81B60', alpha=0.5)

planet_leapfrog = planet(leapfrog, NITERS, TIMESTEP)
print(planet_leapfrog)
ax.plot(planet_leapfrog[:, 0], planet_leapfrog[:, 1], label='leapfrog',
        linewidth=5, color='#00ACC1', alpha=0.5)

legend = plt.legend(frameon=False)
ax.set_title("leapfrog v/s euler: NITERS=%s dt=%s" % (NITERS, TIMESTEP))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
plt.savefig("leapfrog-dt-1e-2.png")
