import numpy as np
import matplotlib.pyplot as plt

def f(vect):
    div = vect[0]**2 + vect[1]**2
    return np.array([vect[2], vect[3], -vect[0]/div, -vect[1]/div])

def run_5c(N, block=False):
    
    h = np.pi / N

    # cache 
    cache = [np.array([1, 0, 0, 1]), np.array([np.cos(h), np.sin(h), -np.sin(h), np.cos(h)])]

    # x and y for plotting
    x = np.array([cache[0][0], cache[1][0]])
    y = np.array([cache[0][1], cache[1][1]])

    i = 2
    lim = 2 * N
    while i < lim:
        vect = -4 * cache[1] + 5 * cache[0] + h * (4 * f(cache[1]) + 2 * f(cache[0]))
        
        # save x and y coordinates
        x = np.append(x, vect[1])
        y = np.append(y, vect[0])

        cache[0] = cache[1]
        cache[1] = vect
        
        i += 1

    print(f'Norm at t = 4pi, N = {N}: {np.linalg.norm(cache[1])}')

    _, ax = plt.subplots()

    ax.plot(x, y)
    plt.title(f"N = {N}")
    plt.show(block=block)

# set blocking on the last run call to prevent all the windows from closing
run_5c(20)
run_5c(40)
run_5c(80, True)