import matplotlib.pyplot as plt
import numpy as np

ts = np.linspace(0, 1, 100)

x0 = -2.0

def rho(x):
    return (0.1 if x < 0 else 
            0.1 + 0.8 * x if 0 <= x <= 1 else
            0.9)    

while x0 <= 2.0:

    color = (1 if x0 < 0.0 else 0.0,
             1 if 0.0 <= x0 <= 1.0 else 0.0,
             1 if x0 > 1 else 0.0)
    plt.plot((1 - 3 * rho(x0)**2)*ts + x0, ts, color=color)
    x0 += 0.1

plt.xlabel("x")
plt.ylabel("t")
plt.title("Characteristics from x = -2 to x = 2")
plt.show()
