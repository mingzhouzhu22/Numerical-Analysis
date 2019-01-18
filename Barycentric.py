#Randy Zhu
#8554123
#Barycentric Formula for interpolation

from math import *
import numpy as np
import scipy.special
import matplotlib.pyplot as plt

'''
Compute the barycentric weights of at a certain node of a given array
Input:  a --- the array of given nodes and their corresponding values;
        j --- the jth node we want to compute the weights of;
        n --- the highest order of all nodes.
Output: weights --- the barycentric weights of the node we chose.
Author: Randy Zhu
Date:   7/8/18
'''

def barycentric_weights(a,j,n):
    weights=1                                       #set the initial weight to be 1
    for k in range(0,n+1):
        if k != j:                                  #when k=!j, time the reciprocal of the difference between the kth node and the node we chose
            weights=weights*(1/(a[j][0]-a[k][0]))
        else:                                       #when k=j, skip to the next loop, in this case times 1
            weights=weights*1
    return weights

'''
Use Barycentric Formula to calculate the interpolation
Input:  x --- the value we want to approximate at;
        a --- the array of given nodes and their corresponding value;
        n --- the highest order of all nodes.
Output: ans --- the approximated value Barycentric Formula calculated.
Author: Randy Zhu
Date:   7/8/18
'''

def barycentric(x,a,n):
    ans         =0          #set the initial answer to be 0
    numerator   =0          #set the numerator to be 0
    denominator =0          #set the denominator to be 0
    l           =[]
    for i in range(0,n+1):
        l.append(a[i][0])
    if x in l:
        for i in range(0,n+1):
            if x==a[i][0]:
                return a[i][1]
    else:
        for j in range(0,n+1):
            numerator   += ((barycentric_weights(a,j,n))/(x-a[j][0]))*a[j][1]           #numerator of Barycentric Formula
            denominator += (barycentric_weights(a,j,n))/(x-a[j][0])                     #denominator of Barycentric Formula
        ans = numerator/denominator
        return ans

###################equidistributed case##################
def equidistributed_barycentric_weights(n,j):
    return (scipy.special.binom(n,j))*((-1)**j)           #returns barycentric weights when nodes are equidistributed

'''
Use Barycentric Formula to calculate the interpolation when nodes are equidistributed
Input:  x --- the value we want to approximate at;
        a --- the array of given nodes and their corresponding value;
        n --- the highest order of all nodes.
Output: ans --- the approximated value for equidistributed nodes
Author: Randy Zhu
Date:   7/11/18
'''

def barycentric_equal(x,a,n):
    ans         =0
    numerator   =0
    denominator =0
    l           =[]
    for i in range(0,n+1):
        l.append(a[i][0])
    if x in l:
        return 1/(1+x**2)
    else:
        for j in range(0,n+1):
            numerator   += ((equidistributed_barycentric_weights(n,j))/(x-a[j][0]))*(a[j][1])         #call equidistributed barycentric weights
            denominator += (equidistributed_barycentric_weights(n,j))/(x-a[j][0])
        ans = numerator/denominator
        return ans
#########################################################

###################cos case###################
def cos_barycentric_weights(n,j):
    if j==0 or j==n:
        return ((-1)**j)/2
    else:
        return (-1)**j           #returns barycentric weights when nodes have the form cos

'''
Use Barycentric Formula to calculate the interpolation when nodes have the form cos
Input:  x --- the value we want to approximate at;
        a --- the array of given nodes and their corresponding value;
        n --- the highest order of all nodes.
Output: ans --- the approximated value for cos nodes
Author: Randy Zhu
Date:   7/11/18
'''

def barycentric_cos(x,a,n):
    ans         =0
    numerator   =0
    denominator =0
    l           =[]
    for i in range(0,n+1):
        l.append(a[i][0])
    if x in l:
        return 1/(1+x**2)
    else:
        for j in range(0,n+1):
            numerator   += ((cos_barycentric_weights(n,j))/(x-a[j][0]))*a[j][1]         #call cos barycentric weights
            denominator += (cos_barycentric_weights(n,j))/(x-a[j][0])
        ans = numerator/denominator
        return ans

##############################################

'''
Use given precise nodes and their corresponding f(x) values to create an array for interpolating polynomial to use, this is for quidistributed case
Input:  n --- the highest order of all nodes.
Output: a --- the array of given nodes and their precise f(x) values
Author: Randy Zhu
Date:   7/11/18
'''

def getarray(n):
    list    = []
    fxlist  = []
    a       = np.zeros(shape=(n+1,2))
    for j in range(0,n+1):
        list.append(-5+j*(10/n))
    for x in list:
        fxlist.append(e**(-(x**2)/5))
    for i in range(0,n+1):
        a[i]=[list[i],fxlist[i]]
    return a

'''
Use given precise nodes and their corresponding f(x) values to create an array for interpolating polynomial to use, this is for cos case
Input:  n --- the highest order of all nodes.
Output: a --- the array of given nodes and their precise f(x) values
Author: Randy Zhu
Date:   7/11/18
'''

def getcosarray(n):
    list    = []
    fxlist  = []
    a       = np.zeros(shape=(n+1,2))
    for j in range(0,n+1):
        list.append(5*cos((j*pi)/n))
    for x in list:
        fxlist.append(1/(1+x**2))
    for i in range(0,n+1):
        a[i]=[list[i],fxlist[i]]
    return a

'''
Use given nodes and their approximated f(x) values to plot the interpolating polnomial, this is for quidistributed case
Input:  n --- the highest order of all nodes.
        e --- the number of points we want to plot, 5000 in homework
Output: the plot of the interpolating polynomial
Author: Randy Zhu
Date:   7/11/18
'''

def plot(n,e):
    xlist=[]
    ylist=[]
    for i in range(0,e+1):
        xlist.append(-5+i*(10/e))
    for x in xlist:
        #ylist.append(barycentric_equal(x,getarray(n),n))
        ylist.append(1/(1+x**2))
    X=xlist
    Y=ylist
    plt.scatter(X,Y)
    plt.show()

'''
Use given nodes and their approximated f(x) values to plot the interpolating polnomial, this is for cos case
Input:  n --- the highest order of all nodes.
        e --- the number of points we want to plot, 5000 in homework
Output: the plot of the interpolating polynomial
Author: Randy Zhu
Date:   7/11/18
'''

def plotcos(n,e):
    xlist=[]
    ylist=[]
    for i in range(0,e+1):
        xlist.append(5*cos((i*pi)/n))
    for x in xlist:
        #ylist.append(barycentric_cos(x,getcosarray(n),n))
        ylist.append(1/(1+x**2))
    X=xlist
    Y=ylist
    plt.scatter(X,Y)
    plt.show()

def lebesgue_weight(l,j,n):                         #take a list as a parameter
    weights=1                                       #set the initial weight to be 1
    for k in range(0,n+1):
        if k != j:                                  #when k=!j, time the reciprocal of the difference between the kth node and the node we chose
            weights=weights*(1/(l[j]-l[k]))
        else:                                       #when k=j, skip to the next loop, in this case times 1
            weights=weights*1
    return weights


#def lebesgue(x,n):
#    L        =0
#    w        =1
#    nodelist=[]
#    for i in range(0,n+1):
#        nodelist.append(-1+i*(2/n))
#    for j in range(0,n+1):
#        w=w*(x-(-1+j*(2/n)))
#    for j in range(0,n+1):
#        L+=((w*(lebesgue_weight(nodelist,j,n)))/(x-nodelist[j]))
#    return L

#def plotlebesgue(n,e):
#    xlist=[]
#    Llist=[]
#    for i in range(0,e+1):
#        xlist.append(-1+i*(2/n))
#    for x in xlist:
#        Llist.append(lebesgue(x,n))
#    plt.scatter(xlist,ylist)
