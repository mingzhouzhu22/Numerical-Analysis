from math import *
import numpy as np

#test arrays:
A = np.array([[1.0,2.0,0.0,0.0],[1.0,1.0,7.0,0.0],[0.0,3.0,3.0,1.0],[0.0,0.0,4.0,4.0]])
b = np.array([6.0,7.0,8.0,9.0])
b.shape = (4,1)

'''
Tridiagonal solver
Input:  A --- a matrix with tridiagonal structure;
        b --- a column vector of RHS value.
Output: x --- list of solutions to the linear system 
Author: Randy Zhu
Date:   2/23/2019
'''

def tridiagonal(A,b):
    y  = [] #y list
    l  = [] #L list
    m  = [] #m list
    n  = A.shape[0]
    m1 = A[0,0] #let m1 = A1, A[0,0] in numpy
    m.append(m1) #m[0] = m1
    for j in range(1,n):
        l.append(A[j,j-1]/m[j-1]) #lj = cj / mj
        m.append(A[j,j]-l[j-1]*A[j-1,j]) #mj+1 = aj+1 - lj * bj
    y.append(b[0,0]) #y1 <- d1
    for j in range(2,n+1):
        y.append(b[j-1,0]-l[j-2]*y[j-2])
    x = np.ones((1,n)) #dummy matrix for solutions
    x[0,n-1] = y[n-1] / m[n-1]
    for j in range(n-2,-1,-1):
        x[0,j] = (y[j]-x[0,j+1]*A[j,j+1]) / m[j]
    return x