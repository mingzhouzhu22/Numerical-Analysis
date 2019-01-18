#Randy Zhu
#8554123
#Natural Spilines

from math import *
import numpy as np
import matplotlib.pyplot as plt

'''
The following four funcitons compute the a,b,c,d coefficients
Input:  hj  --- jth length;
        zjj --- j+1th z;
        zj  --- jth z;
        fjj --- j+1th f(x);
        fj  --- jth f(x).
Output: ans --- the corresponding coefficient.
Author: Randy Zhu
Date:   7/21/18
'''

def a_coe(hj,zjj,zj):
    return (zjj-zj)/(6*(hj))

def b_coe(zj):
    return 0.5*zj

def c_coe(hj,zjj,zj,fjj,fj):
    return ((fjj-fj)/hj)-((hj*(zjj+2*zj))/6)

def d_coe(fj):
    return fj

'''
Compute the natural spline of a given set of points
Input:  datalist --- a list of tuples of points;
Output: alist    --- list of all a coefficients;
        blist    --- list of all b coefficients;
        clist    --- list of all c coefficients;
        dlist    --- list of all d coefficients;
        zlist    --- list of all z values;
Author: Randy Zhu
Date:   7/21/18
'''

def natural_spline(t,datalist):
    alist = []
    blist = []
    clist = []
    dlist = []
    n     = len(datalist)                   #length of the datalist
    h     = []                              #creat a list for length of points
    for i in range(0,n-1):
        h.append(datalist[i+1][0]-datalist[i][0])
    zlist = []
    zlist.append(0)                         #set the first z 0 because of the properties of natural spline
    M = np.zeros(shape=(n-2,n-2))           #creat the matrix of the linear system
    for i in range(0,n-2):                  
        M[i][i]=2*(h[i-1]+h[i])             #fill out the diagonal
    for j in range(0,n-3):
        M[j][j+1]=h[j]                      #fill out the line above the diagonal
    for k in range(1,n-2):
        M[k][k-1]=h[k]                      #fill out the line brlow the diagonal
    D = np.zeros(shape=(n-2,1))             #creat a matrix for d values
    for i in range(0,n-2):
        D[i]=((-6*(datalist[i+1][1]-datalist[i][1]))/h[i])+((6*(datalist[i+2][1]-datalist[i+1][1]))/h[i+1])
    Z = np.linalg.solve(M, D)               #Z matrix is consitituted of z values
    for z in Z:
        zlist.append(z[0])
    zlist.append(0)                         #set the last z 0 because of the properties of natural spline
    for i in range(0,n-1):
        alist.append(a_coe(h[i],zlist[i+1],zlist[i]))
        blist.append(b_coe(zlist[i]))
        clist.append(c_coe(h[i],zlist[i+1],zlist[i],datalist[i+1][1],datalist[i][1]))
        dlist.append(d_coe(datalist[i][1]))
    print("\n","A LIST:",alist,"\n","B LIST:",blist,"\n","C LIST:",clist,"\n","D LIST:",dlist,"\n","Z LIST:",zlist)
    result=0
    for i in range(0,n):
        if (t>=datalist[i][0] and t<datalist[i+1][0]):
            result=alist[i]*((t-datalist[i][0])**3)+blist[i]*((t-datalist[i][0])**2)+clist[i]*(t-datalist[i][0])+dlist[i]
    return result

def plotx():                                #plot the (t,x) graph
    tlist=[]
    xlist=[]
    tlist=np.arange(0.0,3.3,0.01)
    for t in tlist:
        xlist.append(natural_spline(t,[(0,1.5),(0.618,0.9),(0.935,0.6),(1.255,0.35),(1.636,0.2),(1.905,0.1),(2.317,0.5),(2.827,1),(3.33,1.5)]))
    return xlist
    #plt.plot(tlist,xlist)
    #plt.show()

def ploty():                                #plot the (t,y) graph
    tlist=[]
    ylist=[]
    tlist=np.arange(0.0,3.3,0.01)
    for t in tlist:
        ylist.append(natural_spline(t,[(0,0.75),(0.618,0.9),(0.935,1),(1.255,0.8),(1.636,0.45),(1.905,0.2),(2.317,0.1),(2.827,0.2),(3.33,0.25)]))
    return ylist
    #plt.scatter(tlist,ylist)
    #plt.show()

def plotanswer():
    tlist=np.arange(0.0,3.3,0.01)
    a=plotx()
    b=ploty()
    plt.plot(a,b)
    plt.show()