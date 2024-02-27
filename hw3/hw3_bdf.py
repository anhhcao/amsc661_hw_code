from hw3_dirk2 import DIRK2step
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

y0 = np.array([1.0,0.0,0.0])

def func(y): 
    dy = np.zeros(3)
    byz = b*y[1]*y[2]
    cy2 = c*y[1]*y[1]
    ax = a*y[0]
    dy[0] = -ax + byz
    dy[1] = ax - byz - cy2
    dy[2] = cy2
    return dy

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

def NewtonIterBDF2(y, y1, y2, h):
    F = y - (4/3)*y1 + (1/3)*y2 - (2/3)*h*func(y)
    DF = np.identity(3) - (2/3)*h*Jac(y)
    return y - np.linalg.solve(DF,F)

def bdf_step(y1, y2, h):
    y = (4.0/3.0) * y1 - (1.0/3.0)*y2
    for j in range(itermax):
        y = NewtonIterBDF2(y, y1, y2, h)
        if np.linalg.norm(y - (4.0/3.0) * y1 + (1.0/3.0)*y2 - (2/3) * func(y)) < tol:
            break
    return y

def BDF(h):

    Nsteps = int(np.ceil(Tmax/h))

    sol = np.zeros((Nsteps+2,3))
    sol[0,:] = y0
    sol[1,:] = DIRK2step(y0,h)

    start_time = time.time()

    for j in range(Nsteps): # DIRK2
        sol[j+2,:] = bdf_step(sol[j+1,:], sol[j,:], h)

    end_time = time.time()
    t_cpu = end_time - start_time
    return sol, t_cpu

# sol, t_cpu = BDF(h)
# t = np.arange(0,(Nsteps+2)*h,h)
# method_name = "BDF2"
# print(f'method = {method_name:5}, CPUtime = {t_cpu:.6e}')

# # plot the solution
# plt.rcParams.update({'font.size': 22})
# fig, ax = plt.subplots(figsize = (8,2))
# plt.plot(t,sol[:,0],label = "BDF2")
# plt.xlabel("x")
# plt.ylabel("t")
# plt.legend()
# fig.savefig('code/bdf_images/bdf_x_e-1.png')
# #plt.xscale("log")
# fig, ax = plt.subplots(figsize = (8,2))
# plt.plot(t,sol[:,1],label = "BDF2")
# plt.xlabel("y")
# plt.ylabel("t")
# plt.legend()
# fig.savefig('code/bdf_images/bdf_y_e-1.png')
# #plt.xscale("log")
# fig, ax = plt.subplots(figsize = (8,2))
# plt.plot(t,sol[:,2],label = "BDF2")
# plt.xlabel("z")
# plt.ylabel("t")
# plt.legend()
# fig.savefig('code/bdf_images/bdf_z_e-1.png')
# #plt.xscale("log")
# plt.show()