#! /usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

nx = 100
ny = 160
x = np.linspace(-4,2,nx)
y = np.linspace(-4,4,ny)
xg,yg = np.meshgrid(x,y)

z = xg + 1j*yg 

euler = 1 + z

midpoint = 1 + z + 0.5*z*z

rk3 = 1 + z + 0.5*z*z + (z**3)/6

rk4 = 1 + z + 0.5*z*z + (z**3)/6 + (z**4)/24

a = np.array([[0, 0, 0, 0, 0, 0, 0],
              [1/5, 0, 0, 0, 0, 0, 0],
              [3/40, 9/40, 0, 0, 0, 0, 0],
              [44/45, -56/15, 32/9, 0, 0, 0, 0],
              [19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0, 0],
              [9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0, 0],
              [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0]])


b = np.array([5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40])

one_vect = np.array([1, 1, 1, 1, 1, 1, 1])

dopri5 = 1 + z * (
    np.dot(b, one_vect)
    + z * np.dot(b, np.matmul(a, one_vect))
    + (z**2) * np.dot(b, np.matmul(a**2, one_vect))
    + (z**3) * np.dot(b, np.matmul(a**3, one_vect))
    + (z**4) * np.dot(b, np.matmul(a**4, one_vect))
    + (z**5) * np.dot(b, np.matmul(a**5, one_vect))
    + (z**6) * np.dot(b, np.matmul(a**6, one_vect)))

plt.rcParams.update({'font.size': 22})
fig, ax = plt.subplots(figsize=(8,8))

plt.contour(xg,yg,abs(euler),np.arange(2), colors='red')  
plt.contour(xg,yg,abs(midpoint),np.arange(2), colors='blue')  
plt.contour(xg,yg,abs(rk3),np.arange(2), colors='green')
plt.contour(xg,yg,abs(rk4),np.arange(2), colors='yellow')
plt.contour(xg,yg,abs(dopri5),np.arange(2), colors='purple')

plt.contourf(xg,yg,abs(euler),np.arange(2), alpha=0.5, colors='red')  
plt.contourf(xg,yg,abs(midpoint),np.arange(2), alpha=0.5, colors='blue')  
plt.contourf(xg,yg,abs(rk3),np.arange(2), alpha=0.5, colors='green')
plt.contourf(xg,yg,abs(rk4),np.arange(2), alpha=0.5, colors='yellow')
plt.contourf(xg,yg,abs(dopri5),np.arange(2), alpha=0.5, colors='purple')

plt.title("RAS")
plt.xlabel("Re(z)")
plt.ylabel("Im(z)")
ax.set_aspect(1)
plt.grid(color='k', linestyle='--', linewidth=0.5)
plt.show()