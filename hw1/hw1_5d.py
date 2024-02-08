import numpy as np
import matplotlib.pyplot as plt

def f(vect):
    div = vect[0]**2 + vect[1]**2
    return np.array([vect[2], vect[3], -vect[0]/div, -vect[1]/div])

def run_5d(N, block=False):

    h = np.pi / N

    u_old = np.array([1, 0, 0, 1])

    # x and y for plotting
    # ignore (x(0), y(0)) since that just interpolates
    # an unecessary line in the plot
    x = np.array([])
    y = np.array([])

    i = 1
    lim = 4 * N
    while i < lim:
        u = u_old + h * f(u_old + (h / 2) * f(u_old))

        # save x and y coordinates
        x = np.append(x, u[1])
        y = np.append(y, u[0])

        u_old = u

        i += 1

    _, ax = plt.subplots()

    ax.plot(x, y)
    ax.set_aspect('equal')

    plt.title(f'N = {N}')
    plt.show(block=block)

# set blocking on the last run call to prevent all the windows from closing
run_5d(20)
run_5d(40)
run_5d(80, True)