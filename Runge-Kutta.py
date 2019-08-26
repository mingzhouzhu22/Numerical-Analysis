import numpy as np
import matplotlib.pyplot as plt
from decimal import *

def fx(x,y): #derivative of x
    return 1.2*x-x**2-(x*y)/(x+0.2)

def fy(x,y):
    return (1.5*x*y)/(x+0.2)-y

# 4th order Runge-Kutta method
# parameters:
#     size -- k
#     x initial value -- xiv
#     y initial value -- yiv
#     lower bound of t -- lb
#     upper bound of t -- ub
# return: plot of y~x

def rk(k,xiv,yiv,lb,ub):
    xlist = []
    xlist.append(xiv)
    ylist = []
    ylist.append(yiv)
    tlist = []
    t = lb
    tlist.append(lb)
    while t <= ub:
        t = t+k
        tlist.append(t)
    for i in range(len(tlist)):
        k1x = fx(xlist[-1],ylist[-1])
        k1y = fy(xlist[-1],ylist[-1])
        k2x = fx(xlist[-1]+(k*k1x)/2,ylist[-1]+(k*k1y)/2)
        k2y = fy(xlist[-1]+(k*k1x)/2,ylist[-1]+(k*k1y)/2)
        k3x = fx(xlist[-1]+(k*k2x)/2,ylist[-1]+(k*k2y)/2)
        k3y = fy(xlist[-1]+(k*k2x)/2,ylist[-1]+(k*k2y)/2)
        k4x = fx(xlist[-1]+k*k3x,ylist[-1]+k*k3y)
        k4y = fy(xlist[-1]+k*k3x,ylist[-1]+k*k3y)
        xlist.append(xlist[-1]+(k*(k1x+2*k2x+2*k3x+k4x))/6)
        ylist.append(ylist[-1]+(k*(k1y+2*k2y+2*k3y+k4y))/6)
    plt.plot(xlist,ylist)
    plt.show()
    #return xlist,ylist