from random import choice, gammavariate
from copy import deepcopy 
import numpy as np
import scipy.special
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import sys, os
sys.path.append(os.path.realpath(".."))
from lib import fit as ft
from lib.exp import Exp, FittedObj, glabels, plotlabels, gColors, putlabel
from copy import deepcopy
p=0.5
D = { (0,0) : 1 , (1,0) : 1 , (1,1) : p }

from random import random
initiallstate =  { (0,0) : 1 , (1,0) : 1 , (1,1) : p }

def exp(k, p, initx=1):
    x = initx
    t = 0 
    L = [x]
    while x < initx*k:
        
        for _ in range(x):
            x = x+1 if random() < p else x 
        L.append(int(x /initx) )        
    return np.array(L)

def normalize_propb(propb):
    for k in range(propb.shape[0]):
        if np.sum(propb[k]) != 0 :
            propb[k] = propb[k].astype(np.float128) / np.sum(propb[k]) 
    return propb

def transpose_normalize_propb(propb):
    ret = deepcopy(propb)
    ret = ret.T
    for t in range(ret.shape[0]):
        if np.sum(ret[t]) != 0 :
            ret[t] = ret[t].astype(np.float128) / np.sum(ret[t]) 
    return ret


def sim(max_mul_size, max_time, number_of_exp=5000, p =0.6 ):
    empty_propb = np.zeros( shape=(max_mul_size*2, max_time+10) )
    print(empty_propb.shape)
    L = [ ]
    for _ in range(number_of_exp):
        L.append( exp(max_mul_size, p, 1) )

    for vec in L:
        for t,k in enumerate(vec):
            empty_propb[k][t] += 1
    return empty_propb

def f(k , t, p=0.4, _dict =D):
    if (t == 0 and k > 0) or (t< 0) or (k<0) :   
        return 0
    if (k,t) in _dict:
        return _dict[(k,t)]
    else:
        _ret = 0 
        for j in range(int(k/2)+1):
            _ret  += scipy.special.binom( int(k/2) , j) *  p**j *\
                 (1-p)**( int(k/2)  + 1 - j ) * f(k-j,t-1, p=p, _dict=_dict)
    _dict[(k,t)] = _ret
    return _dict[(k,t)]


def estimate_distribute(_shape):
    def _estimate_distribute(empty_propb, p):
        initiallstate_copy = deepcopy(initiallstate)
        ret = deepcopy(empty_propb)

        print ( f" value of p : {p}")        
        # hook
        if ret.shape != _shape:
            ret = np.zeros\
                ( shape=_shape ).astype(np.float128)
                
        for k in range(ret.shape[0]):
            for t in range(ret.shape[1]):
                ret[k][t] = f(k,t, p=p, _dict = initiallstate_copy )
        return normalize_propb(ret).flatten()
    return _estimate_distribute

def model_fit( propb_mes ):
    empty_propb = np.zeros( shape=propb_mes.shape ).astype(np.float128)
    # popt, (pcov, (_range, values)) , MSR  =
    return  ft.g_extract_coef( empty_propb, deepcopy(propb_mes).flatten(), estimate_distribute(empty_propb.shape) , p0=[0.5] )

def plot_propb(propb):
    for i in reversed( range(1,len(propb), 2) ):
        plt.plot(  i + np.arange(propb.shape[-1]) , 0.1*i+  propb[i] , c=next(gColors) )
    # plt.xlabel(r"time ticks")
    # plt.ylabel(r"$ \sum{ \frac{\psi_i}{\psi_0 } } \ge k  $")

if __name__ == "__main__":
    
    prob =  sim(20, 40, number_of_exp=50000)
    plot_propb(normalize_propb(deepcopy( prob) ))
    plt.show()
    
    
    # plt.savefig("./svg/fig2kt.svg")
    # plt.clf()
    # _prob = transpose_normalize_propb(prob)
    # plot_propb(_prob)
    # plt.savefig("./svg/fig3tk.svg")
    
    # X = [f(100, i) for i in range(100)]
    # for i in range(len(X)):
    #     X[i] /= sum(X)
    # for i in range(len(X)):
    #     X[i] /= sum(X)
    # X = np.array(X)
    # plt.plot( X / np.max( X) )
    # plt.show()
    print(prob.shape)
    popt, (pcov, (_range, values)), MSR = model_fit( normalize_propb(deepcopy( prob)))
    print(popt, MSR)

    plot_propb( values.reshape(prob.shape)  )
    plt.show()