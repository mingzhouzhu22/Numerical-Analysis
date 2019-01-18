#Randy Zhu
#8554123
#Newton's Method

from math import *
import numpy as np
import matplotlib.pyplot as plt

'''
Compute the divided differences of given x and f(x) values
Input:  f --- the function to compute nodes;
        j --- number of nodes;
        F --- the function to compute f(x) values.
Output: ans --- the list of all divided differences.
Author: Randy Zhu
Date:   7/15/18
'''

def divided_difference(f,j,F):
    xlist = []                                       #list of nodes
    ylist = []                                       #list of corresponding f(x) values
    clist = [[] for _ in range(0,j+1)]               #creat a general divided difference list containing n sublists for every c(n)
    ans   = []                                       #creat a list for the coefficients that are in the uppermost diagonal
    for i in range(0,j+1):
        xlist.append(float(f(i)))                    #fill out xlist with nodes values
    n     = len(xlist)
    for x in xlist:
        ylist.append(float(F(x)))                    #fill out ylist with precise f(x) values

    for i in range(n):
        ans.append(ylist[i])
    for j in range(1,n):
        for i in range(n-1,j-1,-1):
            ans[i] = float(ans[i]-ans[i-1])/float(xlist[i]-xlist[i-j])
    #second way to do this:
    #for y in ylist:                                  #let 0th divided differences be ylist itself, i.e. f(x) values
    #    clist[0].append(y)
    #for i in range(1,j+1):                           #for every ith divided difference, starting from 1, fill the ith order list with computed values
    #    for k in range(1,len(clist[i-1])):
    #        clist[i].append(float(clist[i-1][k]-clist[i-1][k-1])/float(xlist[i+k-1]-xlist[k-1]))
    #for i in range(0,j+1):
    #    ans.append(clist[i][0])
    return ans

'''
Evaluate the corresponding polynomial at a given x using Newton's form
Input:  x --- the value we want to evaluate
        f --- the function to compute nodes;
        j --- number of nodes;
        F --- the function to compute f(x) values.
Output: p --- the approximated value of the corresponding polynomial.
Author: Randy Zhu
Date:   7/16/18
'''

def newton_form(x,f,j,F):
    xlist   = []
    cn      = divided_difference(f,j,F)
    p       = cn[-1]
    for i in range(0,j+1):
        xlist.append(f(i))
    for i in reversed(range(0,j)):
        p=p*(x-xlist[i])+cn[i]
    return p

#j is the order of p and n is the number of nodes we want to approximate
def plot(f,j,F,n):                              #plot the polynomial of a list of given nodes and their approximated values
    xlist  =[]
    ylist  =[]
    for i in range(0,n+1):
        xlist.append(f(i))
    for x in xlist:
        ylist.append(newton_form(x,lambda j:-1+0.2*j,10,lambda x:e**(-x**2)))      #for every nodes, evaluate its approximation
    plt.plot(xlist,ylist)
    plt.show()

def ploterror(f,j,F,n):                         #function to plot the error
    xlist  =[]
    ylist  =[]
    fxlist =[]
    error  =[]                                  #creat an empty list for error, i.e. f(x)-P(x)
    for i in range(0,n+1):
        xlist.append(f(i))
    for x in xlist:
        ylist.append(newton_form(x,lambda j:-1+0.2*j,10,lambda x:e**(-x**2)))
    for x in xlist:
        fxlist.append(F(x))
    for i in range(0,n+1):
        error.append(fxlist[i]-ylist[i])
    plt.plot(xlist,error)
    plt.show()
    #plt.scatter(xlist,error)
    #plt.show()
