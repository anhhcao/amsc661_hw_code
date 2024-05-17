import numpy as np
import matplotlib.pyplot as plt

k = 0.0001
h = 0.0001

r3 = np.sqrt(3)
M = 2/(3 * r3)

x0 = -2
t0 = 0
xmax = 2
tmax = 1

def f(u):
    return u * (1 - u*u)

# f is a basic polynomial so it is easy to hardcode the min/max
def F(ul, ur):
    if ul <= ur:
        return min(-M, f(ul), f(ur)) if ul <= -r3 <= ur else min(f(ul), f(ur))
    else:
        return max(M, f(ul), f(ur)) if ur <= r3 <= ul else max(f(ul), f(ur))

def U(n, j):


def godunov()