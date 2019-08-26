import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import *

# x --- space

def timezero(x):
    if x <= pi/2 and x >= 0:
        return x
    elif x >= pi/2 and x<= pi:
        return (pi-x)

# k --- nk, time step
# h --- jh, space step
# sigma --- kappa value
# ip --- u0, inital value
# ep --- uM, end value
# M --- bound, space
# N --- bound, time

def fd(k,h,sigma,M,N):
    mlist = np.linspace(0,pi,num=M+1)
    u = np.zeros([N,M])
    alpha = (sigma*k)/(h**2)
    print('Alpha =',alpha)
    for i in range(M):
        u[0,i] = timezero(mlist[i])
    for n in range(1,N):
        for j in range(1,M-1):
            u[n,j] = u[n-1,j]+alpha*(u[n-1,j-1]-2*u[n-1,j]+u[n-1,j+1])
    
    x,y = np.meshgrid(np.arange(0,100),np.arange(0,100))
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_surface(x,y,u)
    plt.show()
    return u

def bd(k,h,sigma,M,N):
    mlist = np.linspace(0,pi,num=M+1)
    u = np.zeros([N,M])
    alpha = (sigma*k)/(h**2)
    print('Alpha =',alpha)
    A = 1/h
    for i in range(M):
        u[0,i] = timezero(mlist[i])
    for n in range(1,N):
        a = (M-1)*[1+2*A]
        b = (M-2)*[-1*A]
        c = b
        d = u[n-1,1:M]
        t = tridiagonal(a,b,c,d)
        u[n,1:M] = t
    x,y = np.meshgrid(np.arange(0,100),np.arange(0,100))
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_surface(x,y,u)
    plt.show()
    return u