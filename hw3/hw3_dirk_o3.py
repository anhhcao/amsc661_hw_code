import numpy as np
import matplotlib.pyplot as plt

# Stiff Robertson's problem from chemical kinetics as in
# https://archimede.uniba.it/~testset/report/rober.pdf
import numpy as np
import matplotlib.pyplot as plt
import time

a = 0.04
b = 1.0e4
c = 3.0e7

# timestep, Tmax, tolearnce for Newton's solver
#h = 1.0e-1
Tmax = 1.0e2 # up to 4.0e10
#Nsteps = int(np.ceil(Tmax/h))
tol = 1.0e-14
itermax = 20

# initial condition
y0 = np.array([1.0,0.0,0.0])

# the right-hand side
def func(y): 
    dy = np.zeros(3)
    byz = b*y[1]*y[2]
    cy2 = c*y[1]*y[1]
    ax = a*y[0]
    dy[0] = -ax + byz
    dy[1] = ax - byz - cy2
    dy[2] = cy2
    return dy

# the Jacobian matrix for the right-hand side
def Jac(y):
    by = b*y[1]
    bz = b*y[2]
    c2y = 2*c*y[1]
    Jac = np.zeros((3,3))
    Jac[0,0] = -a
    Jac[0,1] = bz
    Jac[0,2] = by
    Jac[1,0] = a
    Jac[1,1] = -bz-c2y
    Jac[1,2] = -by
    Jac[2,1] = c2y
    return Jac

# DIRKo3

def NewtonIterDIRKo3(y,h,k,gamma):
    aux = y + h*gamma*k
    F = k - func(aux)
    DF = np.identity(3) - h*gamma*Jac(aux)
    return k - np.linalg.solve(DF,F)

def DIRKo3step(y,h):
    gamma = 0.5 + np.sqrt(3)/6
    k1 = func(y)
    for j in range(itermax):
        k1 = NewtonIterDIRKo3(y,h,k1,gamma)
        if np.linalg.norm(k1 - func(y + h*gamma*k1)) < tol:
            break
    # print( np.linalg.norm(k1 - func(y + h*gamma*k1)))
    # print(j)
    k2 = k1
    y += h*k1/2
    for j in range(itermax):
        k2 = NewtonIterDIRKo3(y,h,k2,gamma)
        if np.linalg.norm(k2 - func(y + h*((1 - 2*gamma)*k1 + gamma*k2))) < tol:
            break
    y += h*k2/2
    # print( np.linalg.norm(k2 - func(aux)))
    # print(j)
    return y

def DIRKo3(h):

    Nsteps = int(np.ceil(Tmax/h))

    sol = np.zeros((Nsteps+1,3))
    sol[0,:] = y0
    start_time = time.time()
    
    for j in range(Nsteps): # DIRK2
        sol[j+1,:] = DIRKo3step(sol[j,:],h)

    end_time = time.time()
    t_cpu = end_time - start_time
    return sol, t_cpu

# sol, t_cpu = DIRKo3(h)

# t = np.arange(0,(Nsteps+1)*h,h)
# method_name = "DIRKo3"
# print(f'method = {method_name:5}, CPUtime = {t_cpu:.6e}')

# # plot the solution
# plt.rcParams.update({'font.size': 22})
# fig, ax = plt.subplots(figsize = (8,2))
# plt.plot(t,sol[:,0],label = "DIRKo3")
# plt.xlabel("x")
# plt.ylabel("t")
# plt.legend()
# fig.savefig('code/dirko3_images/dirko3_x_e-1.png')
# #plt.xscale("log")
# fig, ax = plt.subplots(figsize = (8,2))
# plt.plot(t,sol[:,1],label = "DIRKo3")
# plt.xlabel("y")
# plt.ylabel("t")
# plt.legend()
# fig.savefig('code/dirko3_images/dirko3_y_e-1.png')
# #plt.xscale("log")
# fig, ax = plt.subplots(figsize = (8,2))
# plt.plot(t,sol[:,2],label = "DIRKo3")
# plt.xlabel("z")
# plt.ylabel("t")
# plt.legend()
# fig.savefig('code/dirko3_images/dirko3_z_e-1.png')
# #plt.xscale("log")
# plt.show()