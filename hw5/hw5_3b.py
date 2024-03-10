import numpy as np
import matplotlib.pyplot as plt

z0 = np.array([0, 0.5, 2, 0])

# period and max integration time
T = 9.673596609249161
Tmax = 10 * T

def f(z):
    u, v, x, y = z
    r = np.sqrt(x**2 + y**2)**3
    return np.array([-x/r, -y/r, u, v])

def Df(z):
    _, _, x, y = z
    r = np.sqrt(x**2 + y**2)
    j = np.zeros((4,4))

    j[0, 2] = 3 * x**2 / (r ** 5) - 1/(r**3)
    j[0, 3] = 3 * x * y / (r ** 5)

    j[1, 2] = 3 * x * y / (r ** 5)
    j[1, 3] = 3 * y**2 / (r ** 5) - 1/(r**3)

    j[2, 0] = 1
    
    j[3, 1] = 1

    return j

# parameters for Newton's method
tol = 1.0e-14
itermax = 20

def NewtonIter(z, h, k):
    aux = z + (h/2) * k
    F = k - f(aux)
    DF = np.identity(4) - h * 0.5 * Df(aux)
    return k - np.linalg.solve(DF,F)

def midpointStep(z, h):
    k = np.linalg.solve(np.identity(4) - h * 0.5 * Df(z), f(z))
    
    for _ in range(itermax):
        k = NewtonIter(z, h, k)
        if np.linalg.norm(k - f(z + h * 0.5 * k)) < tol:
            break

    return z + (h * k)

def midpoint(h):
    Nsteps = int(np.ceil(Tmax/h))

    sol = np.zeros((Nsteps+1, 4))
    sol[0,:] = z0

    for j in range(Nsteps):
        sol[j+1,:] = midpointStep(sol[j,:], h)
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
    ax.set_title(label=f'Midpoint, Hamiltonian, h=T/{steps_per_T}')
    ax.set_ylabel('H')
    ax.set_xlabel('t')
    ax.plot(ts, Hs)
    plt.show(block=block)

def plot(steps_per_T, block=False):

    sol = midpoint(T/steps_per_T)

    plotHamiltonian(steps_per_T, sol)

    _, ax = plt.subplots()
    ax.set_title(label=f'Midpoint, xy, h=T/{steps_per_T}')
    ax.set_ylabel('y')
    ax.set_xlabel('x')
    ax.plot(sol[:,2], sol[:,3])
    plt.show(block=block)

plot(100)
plot(1000)
plot(10000, True)