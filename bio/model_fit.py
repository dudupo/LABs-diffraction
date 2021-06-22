
from random import choice, gammavariate
from copy import deepcopy
from re import S 
import numpy as np
from numpy.lib import vectorize
import scipy.special
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import sys, os
sys.path.append(os.path.realpath(".."))
from lib import fit as ft
from lib.exp import Exp, FittedObj, gencolor, glabels, plotlabels, gColors, putlabel
from copy import deepcopy

from bio.mc import get_multi_factor
from bio.sim import sim, exp
from bio.utility import histogram_calc, convert_colonys_to_ktlist
import pickle
from datetime import datetime

import matplotlib
matplotlib.rcParams['text.usetex'] = True


p=0.5
D = { (0,0) : 0 , (1,0) : p , (1,1) : 1-p }

from random import random
initiallstate =  { (0,0) : 0 , (1,0) : 1 , (1,1) : 1 }

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

def f(k , t, p=0.125, N =10, z = 10, _dict =D):
  
    if (t == 0 and k > 0) or (t< 0) or (k<0):   
        return 0
    if (k,t) in _dict:
        return _dict[(k,t)]
    else:
        _ret = 0 
        for j in range(int(k)):
            _ret  += scipy.special.binom( k , j) *  p**j *\
                 (1-p)**( k - j ) * f(k-j,t-1, p=p, _dict=_dict)
    _dict[(k,t)] = _ret
    return _dict[(k,t)]


def estimate_distribute(_shape):
    def _estimate_distribute(empty_propb, p, N, z):
        initiallstate_copy = deepcopy(initiallstate)
        ret = deepcopy(empty_propb)

        print ( f" value of p : {p}")        
        hook = False
        
        # hook
        if ret.shape != _shape:
            ret = np.zeros\
                ( shape=_shape ).astype(np.float128)
            hook = True
                
        for k in range(ret.shape[0]):
            for t in range(ret.shape[1]):
                ret[k][t] = f(k,t, p=p, N=N,z=z, _dict = initiallstate_copy )
        return normalize_propb(ret).flatten() 
    return _estimate_distribute

def model_fit( propb_mes ):
    empty_propb = np.zeros( shape=propb_mes.shape ).astype(np.float128)
    return  ft.g_extract_coef( empty_propb, propb_mes.flatten(),
                estimate_distribute(empty_propb.shape) , p0=[0.14, 10, 10] )

def plot_propb(propb, single = False, k=2, scat =False):
    if not single:
        for i in reversed( range(1,len(propb), 2) ):
            if not scat:
                plt.plot(  i + np.arange(propb.shape[-1]) , 0.1*i+  propb[i] , c=next(gColors) )
            else:
                plt.scatter(  i + np.arange(propb.shape[-1]) , 0.1*i+  propb[i] , c=next(gColors), s=0.1 )
    else:
        plt.plot( np.arange(propb.shape[-1]),  propb[k] , c=next(gColors) )

def picklize():
    prob = get_multi_factor( 200, 0, 100 )
    with open(f"colonys-prob_test-{datetime.now()}.pkl", 'wb') as handle:
        pickle.dump(prob, handle, protocol=pickle.HIGHEST_PROTOCOL)


def mean_func_of_time(propb):
    propb = deepcopy(propb)
    At = np.zeros( propb.shape[-1] )
    for j, vec in enumerate( propb ):
        At += vec * (j + 1)
    return At

def plot_mean_func_of_time(propb):
    # propb = propb.T
    At = mean_func_of_time(propb)
    # plt.scatter(  np.arange( len(At[2:] )) , At[2:] )
    plt.plot( At )

def test( gen_methods ):
    def _test( ):
        ret = [ ]
        firat_color = next(gColors)
        for _method in gen_methods:
            prob = _method()
            ret.append( deepcopy ( prob ) )
            plot_propb( prob )
            current_color = next(gColors)
            while current_color != firat_color :
                current_color = next(gColors)
        return ret
    return _test

def load_data_pkl_filter( ):
    # colonies = pickle.load( open("pkl/colonys-prob_test-2021-06-21_14-27-21.872035.pkl", "rb")) #[:120]
    colonies = pickle.load( open("pkl/colonys-prob_test-2021-06-21_16-55-07.727509.pkl", "rb"))
    good_colonies, _ = filter_colonies(colonies)
    pkt =  histogram_calc(  convert_colonys_to_ktlist(good_colonies), 1, 100, 80)
    return pkt[:50,40:]

# fit the model, and returns pkt.
def lam_fit_pkt(prob ):
    def _lam_fit_pkt( ):
        popt, (pcov, (_range, values)), MSR = model_fit( normalize_propb(deepcopy( prob)))    
        return values.reshape( prob.shape )
    return _lam_fit_pkt

simvsdata = test( [
      lambda : normalize_propb
                    (sim( 50, 100 ))[:80*1],
      lambda : normalize_propb
                    (load_data_pkl_filter())[:80*1] 
    ])
datavsmodel = test( [
      lambda : lam_fit_pkt
                (normalize_propb
                    (load_data_pkl_filter())[:80*1]) (),
      lambda : normalize_propb
                (load_data_pkl_filter())[:80*1] 
    ])
simvsmodel =  test( [
      lambda : normalize_propb
                    (sim( 50, 100 ))[:80*1],
      lambda : lam_fit_pkt
                (normalize_propb
                    (sim( 50, 100 ))[:80*1]) (),
    ])

def colony_k(colonies, i, plot=False):
    colony = colonies[i]
    col0 = colony[0].pix_num
    k_arr = []
    for j, col_t in enumerate(colony):
        k = int(col_t.pix_num/col0 + 0.5)
        k_arr.append(k)
    if plot:
        plt.plot(k_arr)
        # plt.show()
    return k_arr

from bio.mc2 import filter_colonies


def plot_model_p( p, emptyshape=(200,500) ):
    pkt = estimate_distribute(emptyshape)( np.zeros((1,1)) , p , 1, 1).reshape(emptyshape)
    plot_mean_func_of_time(pkt[:, :40])

def battle( _method, genplot_pks = lambda : 1, genplot_mean = lambda : 1 , name=str(random()) ):
    pkts = _method()    
    genplot_pks()
    
    plt.savefig(f"./svg/pkt-{name}.svg")    
    plt.clf()
    for pkt in pkts :
        plot_mean_func_of_time( pkt )
    genplot_mean()
    plt.savefig(f"./svg/average-{name}.svg")
    plt.clf()
if __name__ == "__main__":

    # exit(0)
    def plot_axes():
        plt.xlabel(r"time ticks")
        plt.ylabel(r"$ \sum{ \frac{\psi_i}{\psi_0 } } = k  $")

    def plot_title( _title ):
        def lam():
            plt.title(_title)
            plot_axes()
        return lam
    battle( simvsdata, genplot_pks=plot_title(r"sim vs data") , genplot_mean=plot_title(r"average, sim vs data"), name="sim-vs-data")
    battle( simvsmodel,  genplot_pks=plot_title(r"sim vs model"),  genplot_mean=plot_title(r"average, sim vs model"), name="sim-vs-model")
    battle( datavsmodel, genplot_pks=plot_title(r"data vs model"),  genplot_mean=plot_title(r"average, data vs model"), name="data-vs-model") 
    
    plot_model_p(0.108)
    plot_title( r"model for $ p = 0.108 $ ")()
    plt.show()
    exit(0)



