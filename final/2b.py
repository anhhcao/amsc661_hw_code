import numpy as np 
import matplotlib.pyplot as plt 
import random
from scipy.spatial import Delaunay

N = 20
h = 1/N
omega = []

i = -1
while i <= 1:
    j = -1
    while j <= 1:
        k = j # random.uniform(-1, 1)
        if (i - 0.5)**2 + (k - 0.1)**2 >= 0.04 and (i + 0.5)**2 + (k + 0.1)**2 >= 0.04:
            omega.append([i, k])
        j += h
    i += h

omega = np.array(omega)

tri = Delaunay(omega)

plt.triplot(omega[:,0], omega[:,1], tri.simplices, linewidth=0.5)
plt.plot(omega[:,0], omega[:,1], '.')
plt.axis('equal')
plt.show()