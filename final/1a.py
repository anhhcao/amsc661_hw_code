import numpy as np
import matplotlib.pyplot as plt
import time

mu = 100
T = 200
h = 0.001

Nsteps = int(np.ceil(T/h))
tol = 1.0e-14
itermax = 20

# initial condition
v0 = np.array([2, 0])

# the right-hand side
def func(v): 
    dv = np.zeros(2)
    dv[0] = v[1]
    dv[1] = mu * (1 - v[0] ** 2) * v[1] - v[0]
    return dv

# the Jacobian matrix for the right-hand side
def Jac(v):
    Jac = np.zeros((2, 2))
    Jac[0,0] = 0
    Jac[0,1] = 1
    Jac[1,0] = -2 * mu * v[0] * v[1] - 1
    Jac[1,1] = mu * (1 - v[0]**2)
    return Jac

# DIRK2

def NewtonIterDIRK2(y,h,k,gamma):
    aux = y + h*gamma*k
    F = k - func(aux)
    DF = np.identity(2) - h*gamma*Jac(aux)
    return k - np.linalg.solve(DF,F)

def DIRK2step(y,h):
    gamma = 1.0 - 1.0/np.sqrt(2)
    k1 = func(y)
    # print(y)
    # print(k1)
    # exit()
    for j in range(itermax):
        k1 = NewtonIterDIRK2(y,h,k1,gamma)
        if np.linalg.norm(k1 - func(y + h*gamma*k1)) < tol:
            break
    # print( np.linalg.norm(k1 - func(y + h*gamma*k1)))
    # print(j)
    k2 = k1
    y = y + h*(1-gamma)*k1
    for j in range(itermax):
        k2 = NewtonIterDIRK2(y,h,k2,gamma)
        aux = y + h*gamma*k2
        if np.linalg.norm(k2 - func(aux)) < tol:
            break
    # print( np.linalg.norm(k2 - func(aux)))
    # print(j)
    return aux

def DIRK2(h):
    Nsteps = int(np.ceil(T/h))
    sol = np.zeros((Nsteps+1,2))
    sol[0,:] = v0
    start_time = time.time()

    for j in range(Nsteps): # DIRK2
        sol[j+1,:] = DIRK2step(sol[j,:],h)

    end_time = time.time()
    t_cpu = end_time - start_time
    return sol, t_cpu

t = np.arange(0,(Nsteps+1)*h,h)
method_name = "DIRK2"
sol, t_cpu = DIRK2(h)

print(f'method = {method_name:5}, CPUtime = {t_cpu:.6e}')

# plot the solution
fig, ax = plt.subplots()
ax.set_title("x vs t")
plt.plot(sol[:,0], t)
plt.xlabel("x")
plt.ylabel("t")

fig, ax = plt.subplots()
ax.set_title("y vs t")
plt.plot(sol[:,1], t)
plt.xlabel("y")
plt.ylabel("t")

fig, ax = plt.subplots()
ax.set_title("y vs x")
plt.plot(sol[:,1], sol[:,0])
plt.xlabel("y")
plt.ylabel("x")

plt.show()