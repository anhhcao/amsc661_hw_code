from hw3_dirk_o3 import DIRKo3
from hw3_dirk2 import DIRK2
from hw3_bdf import BDF
import matplotlib.pyplot as plt
import numpy as np
import time

hs = [10**(-3), 10**(-2), 10**(-1)]
log_hs = [np.log(10**(-3)), np.log(10**(-2)), np.log(10**(-1))]
dirk2_time = []
dirko3_time = []
bdf_time = []

for h in hs:
    
    start = time.time()
    DIRK2(h)
    end = time.time()
    dirk2_time.append(np.log(end - start))
    print(f"method: DIRK2, CPU Time: {end - start}, h:{h}")

    start = time.time()
    DIRKo3(h)
    end = time.time()
    dirko3_time.append(np.log(end - start))
    print(f"method: DIRKo3, CPU Time: {end - start}, h:{h}")

    start = time.time()
    BDF(h)
    end = time.time()
    bdf_time.append(np.log(end - start))
    print(f"method: BDF, CPU Time: {end - start}, h:{h}")

fig, ax = plt.subplots()
ax.title.set_text("DIRK2 CPU Time vs h (log-log)")
plt.plot(log_hs, dirk2_time)

fig, ax = plt.subplots()
ax.title.set_text("DIRKo3 CPU Time vs h (log-log)")
plt.plot(log_hs, dirko3_time)

fig, ax = plt.subplots()
ax.title.set_text("BDF CPU Time vs h (log-log)")
plt.plot(log_hs, bdf_time)

plt.show()