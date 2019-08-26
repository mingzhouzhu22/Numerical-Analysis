from math import *
import numpy as np

A = np.array([[5.0,1.0,0.0,2.0,1.0],[0.0,4.0,0.0,1.0,2.0],[1.0,1.0,4.0,1.0,1.0],[0.0,1.0,2.0,6.0,0.0],[0.0,0.0,1.0,2.0,4.0]])
b = np.array([1.0,2.0,3.0,4.0,5.0])
b.shape = (5,1)

'''
Gaussian Elimination with partial pivoting to solve a linear system
Input:  A --- a m * j matrix;
        b --- a m * 1 vector.
Output: x --- list of solutions to the linear system 
Author: Randy Zhu
Date:   2/12/2019
'''

'''
test pseudo code:
    for j in range(n):
        if abs(A[j,0]) > col_max:
            col_max = A[j,0]
            max_row = j
        if col_max == 0:
            print("Error: matrix A is nonsingular.")
        else:
            A[[0,j]]=A[[j,0]] #switch rows of A
        for i in range(j+1,n):
            m = (A[i,j])/(A[j,j])
            print(m)
            for k in range(j+1,n):
                A[i,k] -= (m * A[j,k])
            A[i,j] = m
    for i in reversed(range(0,n)):
    sigma = 0
    for j in range(i+1,n+1):
        sigma += A[i,j]*x[j]
    x[i] = (A[i,n+1]-sigma)/A[i,i]
    return x


def gaussian(A,b):
    if np.linalg.det(A) == 0: #nonsingular since determinant is 0
        print("Error: nonsingular")
    else:
        dim = A.shape
        n = dim[0] - 1
        for j in range(n):
            m = j #pseudocode
            max_col = abs(A[j,j]) #set A(1,1) to be the default max column value
            for i in range(j+1,n+1):
                if abs(A[i,j]) > max_col:
                    m = i
                    max_col = abs(A[i,j])
            if max_col == 0:
                continue
            for k in range(j,n+1):
                t = A[j,k]
                A[j,k] = A[m,k]
                A[m,k] = t
            t = b[j]
            b[j] = b[m]
            b[m] = t
            for i in range(j+1,n+1):
                m = A[i,j]/A[j,j] #multiplier
                for k in range(j+1,n+1):
                    A[i,k] = A[i,k] - m * A[j,k]
                A[i,j] = 0
                b[i] = b[i] - m*b[j]
        x = np.ones((1,n+1)) #dummy matrix for solutions
        for i in range(n,-1,-1):
            x[0,i] = b[i]
            for j in range(i+1,n+1):
                x[0,i] = x[0,i]-A[i,j]*x[0,j]
            x[0,i] = x[0,i]/A[i,i]
        return x

'''
def gaussian(A,b):
    if np.linalg.det(A) == 0:
        print("Error: nonsingular")
    A   = np.concatenate((A,b),1)
    dim = A.shape #get dimension of A
    n   = dim[0] #number of rows
    c   = dim[1] #number of columns (n+1) since column vector
    x   = np.ones(n) #solutions of x
    col_max = abs(A[0,0]) #set the max pivot A11 (0,0 in numpy)
    max_row = 0 #ith row with the max pivot
    for i in range(n):
        if abs(A[i,0]) > col_max:
            col_max = A[i,0]
            max_row = i
    #after the for loop, max row has the ith row with the greatest value
    A[[0,max_row]] = A[[max_row,0]] #switch rows
    for i in range(n-1):
        for j in range(i+1,n):
            m = A[j,i]/A[i,i] #get multiplier
            A[j,i:] = A[j,i:] - m * A[i,i:] #target row minus ith row times multiplier
    #now A is already Gaussian eliminated
    for i in reversed(range(n)):
        sigma = 0
        for j in range(i+1,n):
            sigma += A[i,j]*x[j]
        x[i] = (A[i,n]-sigma)/A[i,i]
    acc = 0
    for s in x:
        if s == nan or s == inf or s == -inf:
            acc += 1
    if acc == 0:
        return x
    else:
        print('Error: nonsingular.')
