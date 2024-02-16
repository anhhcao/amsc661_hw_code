#! /usr/bin/python

from scipy.integrate import solve_ivp
from time import process_time_ns
import matplotlib.pyplot as plt
import numpy as np

def plot_arenstorf(T):
    mu = 0.012277471
    tspan = np.array([0, T])
    initial = [0.994, 0.0, 0.0, -2.00158510637908252240537862224]

    def f(t, y):
        d1 = ((y[0] + mu) ** 2 + y[1] ** 2) ** 1.5
        d2 = ((y[0] - 1 + mu) ** 2 + y[1] ** 2) ** 1.5
        return [y[2], 
                y[3], 
                y[0] + 2 * y[3] - (1 - mu) * (y[0] + mu)/d1 - mu * (y[0] - 1 + mu) / d2,
                y[1] - 2 * y[2] - (1 - mu) * y[1] / d1 - mu * y[1] / d2]

    t_start = process_time_ns()
    sol = solve_ivp(f, tspan, initial, 'RK45', rtol=1e-12, atol=1e-12)
    t_end = process_time_ns()
    print(f'RK45 CPU Time for T={T}: {t_end - t_start}')
    plt.plot(sol.y[0], sol.y[1], 'r')

    t_start = process_time_ns()
    sol = solve_ivp(f, tspan, initial, 'DOP853', rtol=1e-12, atol=1e-12)
    t_end = process_time_ns()
    print(f'RK45 CPU Time for T={T}: {t_end - t_start}')
    plt.plot(sol.y[0], sol.y[1], 'g')

    t_start = process_time_ns()
    sol = solve_ivp(f, tspan, initial, 'Radau', rtol=1e-12, atol=1e-12)
    t_end = process_time_ns()
    print(f'RK45 CPU Time for T={T}: {t_end - t_start}')
    plt.plot(sol.y[0], sol.y[1], 'b')

    plt.title(f'Arenstorf Orbit, T={T}')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend(['RK45', 'DOP853', 'Radau'])
    plt.show()

# plot_arenstorf(17.0652165601579625588917206249)
plot_arenstorf(100)