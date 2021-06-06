

import numpy as np 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def g_derivate_err(time, distance, f, eps = 0.00001):
    return (f(time + eps) - f(time + eps)) / eps

def g_extract_coef( time, distance, f ,p0=None, bounds=None):
    bounds = (- np.inf, np.inf) if bounds is None else bounds
    popt, pcov = curve_fit(f, time, distance, p0=p0, maxfev=100000, bounds=bounds)
    _range = np.linspace( np.min(time) , np.max(time), 1000 )
    values = f(_range, *popt)
    val_in_given_range = f(time, *popt)
    MSR = np.sqrt(np.sum((np.array(val_in_given_range) - distance)**2))/ len(values)
    # print(MSR)
    return popt, (pcov, (_range, values)) , MSR 

from itertools import combinations
def g_extract_coef_drop( time, distance, f, drop_points, p0=None, bounds=None ):
    _length = len(time) - drop_points
    indices = np.arange( len(time) )

    time, distance = np.array(time), np.array(distance) 

    return min ([ g_extract_coef( np.take(time, _indices, 0), np.take(distance, _indices, 0), f, p0=p0, bounds=bounds )\
         for _indices in  combinations( indices, _length)], key= lambda ret : ret[-1])

def extract_coef( time, distance ):
    def f(x, a, b, c, d):
        return ( b* np.sin(a*x+c  )/(a*x+c ) )**2  
    return g_extract_coef(time, distance, f, p0=[ 10, 2, 0, 0 ])

def extract_coef2( time, distance ):
    def f(x, a, b, c, c2, a2): # c2, d, e, e2):
        return ( b* np.sin(a*x+c)/(a*x+c) *  np.cos(a2*x+c2))**2 
        # return ( b* np.sin(a*np.sin(x) + c)/(a*x + d) *  np.cos(a2*np.sin(x) +c2 ))**2   
    
    return g_extract_coef(time, distance, f, p0=[ 10, 1, 0, 0, 200  ] )  

    
