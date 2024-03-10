import numpy as np
import matplotlib.pyplot as plt

def compute_points(alpha):
    beta_m1 = (4 - alpha)/12
    beta_0 = 2 * (alpha + 2) / 3
    beta_1 = 2 + alpha - beta_m1 - beta_0

    re = []
    im = []

    theta = 0
    while theta < 2 * np.pi:
        f = (np.exp(theta * 2j) + alpha * np.exp(theta * 1j) - (1 + alpha))/(beta_m1 * np.exp(theta * 2j) + beta_0 * np.exp(theta * 1j) + beta_1)
        re.append(f.real)
        im.append(f.imag)
        theta += 0.0001

    return re, im

def plot(alpha_min, alpha_max, block=False):
    _, ax = plt.subplots()

    alpha = alpha_min
    legend = []
    while alpha < alpha_max:
        re, im = compute_points(alpha)
        ax.plot(re, im)
        legend.append(f'α = {round(alpha, 1)}')
        alpha += 0.1

    ax.set_xlabel('Re(hλ)')
    ax.set_ylabel('Im(hλ)')
    ax.set_aspect('equal')
    ax.set_title(f'{alpha_min} ≤ α ≤ {alpha_max}')
    ax.legend(legend)
    plt.show(block=block)

def check_root(alpha, z):

    beta_m1 = (4 - alpha)/12
    beta_0 = 2 * (alpha + 2) / 3
    beta_1 = 2 + alpha - beta_m1 - beta_0
    
    a = 1 - z*beta_m1
    b = alpha - z*beta_0
    c = (- 1 - alpha) - z * beta_1

    return (-b - np.sqrt((b**2) - 4 * a * c)) / (2 * a)

i = -1.8
while i < 0:
    print(f'alpha = {round(i, 1)}, h*lambda = {0.15}: {check_root(i, 0.15)}')
    i += 0.1

plot(-1.8, -1, True)
plot(-1, -0, True)