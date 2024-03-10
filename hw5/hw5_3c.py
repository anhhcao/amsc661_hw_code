import numpy as np
import matplotlib.pyplot as plt

z0 = np.array([0, 0.5, 2, 0])

# period and max integration time
T = 9.673596609249161
Tmax = 10 * T

def u_prime(q):
    x, y = q
    r = np.sqrt(x**2 + y**2)**3
    return np.array([x/r, y/r])

def t_prime(p):
    u, v = p
    return np.array([u, v])

def q_step(p, q, h):
    return q + h * t_prime(p)

def p_step(p, q, h):
    return p - h * 0.5 * u_prime(q)

def sv_step(p, q, h):
    p_half_step = p_step(p, q, h)
    q_full_step = q_step(p_half_step, q, h)
    p_full_step = p_step(p_half_step, q_full_step, h)
    u, v = p_full_step
    x, y = q_full_step
    return [u, v, x, y]

def sv(h):

    Nsteps = int(np.ceil(Tmax/h))

    sol = np.zeros((Nsteps+1, 4))
    sol[0,:] = z0

    for j in range(Nsteps):
        p = np.array([sol[j, 0], sol[j, 1]])
        q = np.array([sol[j, 2], sol[j, 3]])
        sol[j+1,:] = sv_step(p, q, h)
    return sol

def plotHamiltonian(steps_per_T, sol, block=False):

    def H(u, v, x, y):
        return 0.5 * u ** 2 + 0.5 * v ** 2 - (1/(np.sqrt(x**2 + y ** 2)))

    Nsteps = 10 * steps_per_T
    h = T/steps_per_T
    ts = np.arange(0, (Nsteps + 1) * h, h)
    Hs = []
    i = 0
    for _ in ts:
        Hs.append(H(sol[i,0], sol[i,1], sol[i,2], sol[i,3]))
        i += 1

    _, ax = plt.subplots()
    ax.set_title(label=f'Stoermer-Verlet, Hamiltonian, h=T/{steps_per_T}')
    ax.set_ylabel('H')
    ax.set_xlabel('t')
    ax.plot(ts, Hs)
    plt.show(block=block)

def plot(steps_per_T, block=False):

    sol = sv(T/steps_per_T)

    plotHamiltonian(steps_per_T, sol)

    _, ax = plt.subplots()
    ax.set_title(label=f'Stoermer-Verlet, xy, h=T/{steps_per_T}')
    ax.set_ylabel('y')
    ax.set_xlabel('x')
    ax.plot(sol[:,2], sol[:,3])
    plt.show(block=block)

plot(100)
plot(1000)
plot(10000, True)