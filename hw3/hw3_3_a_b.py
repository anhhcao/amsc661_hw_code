#! /usr/bin/python

import numpy as np
import matplotlib.pyplot as plt


def dirk2_err(h):

    pi4 = np.pi/4

    T = 10
    g = 1 - (1 / np.sqrt(2))
    L = 10000
    # y0 = np.sin(pi4)
    y0 = np.sin(pi4) + 10
    
    def y(t):
        return (np.e ** (- L * t)) * 10 + np.sin(t + pi4)
        #return np.sin(t + pi4)

    u = y0
    e = 0

    t = 0
    while t < T:
        k1 = (-L * (u - np.sin(t + pi4)) + np.cos(t + pi4)) / (L * h * g + 1)
        k2 = (-L * (u + h*(1 - g)*k1 - np.sin(t + pi4)) + np.cos(t + pi4)) / (L * h * g + 1)
        u += h*(((1 - g) * k1) + (g * k2))
        e_new = np.abs(u - y(t))
        if e_new > e:
            e = e_new
        t += h
    
    return e

def dirk3_err(h):
    
    pi4 = np.pi/4

    T = 10
    g = 0.5 + np.sqrt(3)/6
    L = 10000
    # y0 = np.sin(pi4)
    y0 = np.sin(pi4) + 10
    
    def y(t):
        return (np.e ** (- L * t)) * 10 + np.sin(t + pi4)
        #return np.sin(t + pi4)

    u = y0
    e = 0

    t = 0
    while t < T:
        k1 = (-L * (u - np.sin(t + pi4)) + np.cos(t + pi4)) / (L * h * g + 1)
        k2 = (-L * (u + h*(1 - 2*g)*k1 - np.sin(t + pi4)) + np.cos(t + pi4)) / (L * h * g + 1)
        u += h*(k1/2 + k2/2)
        e_new = np.abs(u - y(t))
        if e_new > e:
            e = e_new
        t += h
    
    return e

# this is for part a and the beginning of part b
p = 1
d = 5/24
es = []
hs = []
while p <= 6:
    h = 10**(-p)
    hs.append(np.log(h))
    #es.append(np.log(dirk2_err(h)))
    es.append(np.log(dirk3_err(h)))
    # print(p)
    p += d

_, ax = plt.subplots()

#ax.plot(hs, es, label='DIRK2 (log)')
ax.plot(hs, es, label='DIRKo3 (log)')
ax.plot(np.linspace(-14, -8), np.linspace(-2.5, 3), label='slope 2')
ax.plot(np.linspace(-8, -2), np.linspace(3, 3), label='slope 0')
plt.legend()
plt.show()