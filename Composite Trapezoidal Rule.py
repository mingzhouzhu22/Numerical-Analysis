#Randy Zhu
#8554123
#Composite Trapezoidal Rule

from math import *

'''
Composite Trapezoidal Rule for a function defined on [a, b]
Input:  a --- the left endpoint of the interval;
        b --- the right endpoint of the interval;
        N --- the number of the partion used in the CTR.
Output: T --- the approximation of the interval given by CTR.
Author: Randy Zhu
Date:   6/26/18
'''

def CTR(f,a,b,N):                     #f has the form lambda x: x, e.g. f=lambda x:x**2
    x=a                               #let x0 equal to a
    h=(b-a)/N                         #set the length
    T=0                               #set the initial value of T
    T+=(f(a)/2)*h
    for i in range(1,N):              #loop through N times
        T+=h*f(x+h)
        x+=h
    T+=(f(b)/2)*h
    return T

#test:
def q():
    return (CTR(lambda x:e**(-x**2),0,1,160)-CTR(lambda x:e**(-x**2),0,1,80))/((CTR(lambda x:e**(-x**2),0,1,320))-(CTR(lambda x:e**(-x**2),0,1,160)))

def E():
    return (4*(CTR(lambda x:e**(-x**2),0,1,160)-CTR(lambda x:e**(-x**2),0,1,80)))/3