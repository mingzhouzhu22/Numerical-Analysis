from math import *
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d

def tridiag(A,b):
    y = [] #y list
    l = [] #L list
    m = [] #m list
    n = A.shape[0]
    m1 = A[0,0] #let m1 = A1, A[0,0] in numpy
    m.append(m1) #m[0] = m1
    for j in range(1,n):
        l.append(A[j,j-1]/m[j-1]) #lj = cj / mj
        m.append(A[j,j]-l[j-1]*A[j-1,j]) #mj+1 = aj+1 - lj * bj
    y.append(b[0]) #y1 <- d1
    for j in range(2,n+1):
        y.append(b[j-1]-l[j-2]*y[j-2])
    x = np.ones((1,n)) #dummy matrix for solutions
    x[0,n-1] = y[n-1] / m[n-1]
    for j in range(n-2,-1,-1):
        x[0,j] = (y[j]-x[0,j+1]*A[j,j+1]) / m[j]
    return x

def u0(x,y): #initial condition
    return sin(pi*x)*sin(pi*y)

def adi(Nx,Nt):
    t = np.linspace(0,1,Nt+1) #time values
    x = np.linspace(0,1,Nx+1) #space values
    y = x #assume h=k
    h = x[1]-x[0]
    k = t[1]-t[0]
    a = k/h**2
    print('Alpha =',a)
    u = np.zeros((Nt+1,Nx+1,Nx+1))
    un = np.zeros((Nt+1,Nx+1,Nx+1))
    u_s = np.zeros((Nt,Nx+1,Nx+1))
    d = []
    b = []
    c = []
    for i in range(Nx-2):
        d.append(-a/2)
        c.append(-a/2)
    for i in range(Nx-1):
        b.append(1+a)
    A = np.diag(d,-1)+np.diag(b,0)+np.diag(c,1)
    for i in range(Nx+1): #at t=0, plug in initial condition formula
        for j in range(Nx+1):
            un[0][i][j] = u0(x[i],y[j]) #push initial values to final matrix
    for n in range(Nt):
        for i in range(1,Nx):
            b = []
            for j in range(1,Nx):
                b.append(a/2*un[n][i-1][j]+(1-a)*un[n][i][j]+a/2*un[n][i+1][j]) #star step
            star = tridiag(A,b)
            star = np.append(star,0)
            star = np.insert(star,0,[0])
            u_s[n][i][:] = star
        for j in range(1,Nx):
            b = []
            for i in range(1,Nx):
                b.append(a/2*u_s[n][i-1][j]+(1-a)*u_s[n][i][j]+a/2*u_s[n][i+1][j])
            ut = tridiag(A,b)
            ut = np.append(ut,0)
            ut = np.insert(ut,0,[0])
            un[n+1][:][j] = ut
    return un