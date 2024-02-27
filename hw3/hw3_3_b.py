import numpy as np
import matplotlib.pyplot as plt

def dirk2_err(h, T):

    pi4 = np.pi/4

    g = 1 - (1 / np.sqrt(2))
    L = 10000
    # y0 = np.sin(pi4)
    y0 = np.sin(pi4) + 10
    
    def y(t):
        return (np.e ** (- L * t)) * 10 + np.sin(t + pi4)
        #return np.sin(t + pi4)

    u = y0
    es = [0]
    ts = [0]

    t = 0
    while t < T:
        k1 = (-L * (u - np.sin(t + pi4)) + np.cos(t + pi4)) / (L * h * g + 1)
        k2 = (-L * (u + h*(1 - g)*k1 - np.sin(t + pi4)) + np.cos(t + pi4)) / (L * h * g + 1)
        u += h*(((1 - g) * k1) + (g * k2))
        es.append(np.log(np.abs(u - y(t))))
        t += h
        ts.append(t)
    
    return ts, es

def dirk3_err(h, T):
    
    pi4 = np.pi/4

    g = 0.5 + np.sqrt(3)/6
    L = 10000
    # y0 = np.sin(pi4)
    y0 = np.sin(pi4) + 10
    
    def y(t):
        return (np.e ** (- L * t)) * 10 + np.sin(t + pi4)
        #return np.sin(t + pi4)

    u = y0
    es = [0]
    ts = [0]

    t = 0
    while t < T:
        k1 = (-L * (u - np.sin(t + pi4)) + np.cos(t + pi4)) / (L * h * g + 1)
        k2 = (-L * (u + h*(1 - 2*g)*k1 - np.sin(t + pi4)) + np.cos(t + pi4)) / (L * h * g + 1)
        u += h*(k1/2 + k2/2)
        es.append(np.log(np.abs(u - y(t))))
        t += h
        ts.append(t)
    
    return ts, es

# dirk2, T=10
_, ax = plt.subplots()
ts, es = dirk2_err(10**(-1), 10)
ax.plot(ts, es)
ax.set_title(f'DIRK2, h={1e-1}, T={10}')
plt.xlabel("t")
plt.ylabel("log(|e(t)|)")

_, ax = plt.subplots()
ts, es = dirk2_err(10**(-2), 10)
ax.plot(ts, es)
ax.set_title(f'DIRK2, h={1e-2}, T={10}')
plt.xlabel("t")
plt.ylabel("log(|e(t)|)")

_, ax = plt.subplots()
ts, es = dirk2_err(10**(-3), 10)
ax.plot(ts, es)
ax.set_title(f'DIRK2, h={1e-3}, T={10}')
plt.xlabel("t")
plt.ylabel("log(|e(t)|)")

# dirk2, T=1
_, ax = plt.subplots()
ts, es = dirk2_err(10**(-1), 1)
ax.plot(ts, es)
ax.set_title(f'DIRK2, h={1e-1}, T={1}')
plt.xlabel("t")
plt.ylabel("log(|e(t)|)")

_, ax = plt.subplots()
ts, es = dirk2_err(10**(-2), 1)
ax.plot(ts, es)
ax.set_title(f'DIRK2, h={1e-2}, T={1}')
plt.xlabel("t")
plt.ylabel("log(|e(t)|)")

_, ax = plt.subplots()
ts, es = dirk2_err(10**(-3), 1)
ax.plot(ts, es)
ax.set_title(f'DIRK2, h={1e-3}, T={1}')
plt.xlabel("t")
plt.ylabel("log(|e(t)|)")

# dirko3, T=10
_, ax = plt.subplots()
ts, es = dirk3_err(10**(-1), 10)
ax.plot(ts, es)
ax.set_title(f'DIRKo3, h={1e-1}, T={10}')
plt.xlabel("t")
plt.ylabel("log(|e(t)|)")

_, ax = plt.subplots()
ts, es = dirk3_err(10**(-2), 10)
ax.plot(ts, es)
ax.set_title(f'DIRKo3, h={1e-2}, T={10}')
plt.xlabel("t")
plt.ylabel("log(|e(t)|)")

_, ax = plt.subplots()
ts, es = dirk3_err(10**(-3), 10)
ax.plot(ts, es)
ax.set_title(f'DIRKo3, h={1e-3}, T={10}')
plt.xlabel("t")
plt.ylabel("log(|e(t)|)")

# dirko3, T=1
_, ax = plt.subplots()
ts, es = dirk3_err(10**(-1), 1)
ax.plot(ts, es)
ax.set_title(f'DIRKo3, h={1e-1}, T={1}')
plt.xlabel("t")
plt.ylabel("log(|e(t)|)")

_, ax = plt.subplots()
ts, es = dirk3_err(10**(-2), 1)
ax.plot(ts, es)
ax.set_title(f'DIRKo3, h={1e-2}, T={1}')
plt.xlabel("t")
plt.ylabel("log(|e(t)|)")

_, ax = plt.subplots()
ts, es = dirk3_err(10**(-3), 1)
ax.plot(ts, es)
ax.set_title(f'DIRKo3, h={1e-3}, T={1}')
plt.xlabel("t")
plt.ylabel("log(|e(t)|)")

plt.show()