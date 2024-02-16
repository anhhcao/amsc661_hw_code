#! /usr/bin/python

from scipy.integrate import solve_ivp
from time import process_time_ns
from math import log
import numpy as np
import matplotlib.pyplot as plt

def plot_van_der_pol(mu, block=False):

    # time t, ndarray y
    def f(t, y):
        return [y[1], mu * ((1 - y[0] ** 2) * y[1]) - y[0]]

    fig, ax = plt.subplots(2, 4)
    fig.tight_layout()
    fig.subplots_adjust(top = 0.88)
    fig.suptitle(f'Van der Pol Oscillator, μ = {mu}')

    i = 0

    # tolerances
    tols = [1e-6, 1e-9, 1e-12]

    cpu_time = []

    while i < 3:

        t_start = process_time_ns()
    
        sol = solve_ivp(f, np.array([0, 1000]), [2, 0], 'RK45', rtol=tols[i], atol=tols[i])

        t_end = process_time_ns()

        cpu_time.append(log(t_end - t_start))

        ax[0, i].plot(sol.y[0], sol.y[1])
        ax[0, i].set_title(f'RK45, ɛ = {tols[i]}')

        i += 1

    ax[0, 3].set_title(f'RK45, log(ɛ) v.s. log(CPU time)')
    ax[0, 3].plot([log(1e-6), log(1e-9), log(1e-12)], cpu_time)
    
    i = 0

    while i < 3:

        t_start = process_time_ns()

        sol = solve_ivp(f, np.array([0, 1000]), [2, 0], 'LSODA', rtol=tols[i], atol=tols[i])

        t_end = process_time_ns()

        cpu_time[i] = log(t_end - t_start)

        ax[1, i].plot(sol.y[0], sol.y[1])
        ax[1, i].set_title(f'LSODA, ɛ = {tols[i]}')

        i += 1

    ax[1, 3].set_title(f'LSODA, log(ɛ) v.s. log(CPU time)')
    ax[1, 3].plot([log(1e-6), log(1e-9), log(1e-12)], cpu_time)

    plt.show()

#plot_van_der_pol(10)
#plot_van_der_pol(100)
plot_van_der_pol(1000)