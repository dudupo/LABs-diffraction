
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
from lib.exp import Exp, FittedObj, gencolor, glabels, plotlabels, gColors, putlabel, resetgencolor
from copy import deepcopy

from bio.mc import get_multi_factor
from bio.sim import sim, exp
from bio.utility import histogram_calc, convert_colonys_to_ktlist
import pickle
from datetime import datetime

import matplotlib
matplotlib.rcParams['text.usetex'] = True
matplotlib.rc('font', size=16)


p=0.5
D = { (0,0) : 0 , (1,0) : p }#, (1,1) : 1-p }

from random import random
initiallstate =  { (0,0) : 0 , (1,0) : 1 , (1,1) : 1-p  }

def ginitiallstate_copy(p1):
    return { (0,0) : 0 , (1,0) : 0, (1,1) : 1 }

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
  
    if (k,t) in _dict:
        return _dict[(k,t)]
    
    if (t == 0 and k > 0) or (t< 0) or (k<0):   
        return 0
    
    else:
        _ret = 0 
        for j in range(int(k/2) + 1):
            # scipy.special.factorial(j) *  scipy.special.binom( k -j , j) *
            _ret  +=   scipy.special.binom( k -j , j) * p**j *\
                 (1-p)**( k - 2*j ) * f(k-j,t-1, p=p, _dict=_dict)
    _dict[(k,t)] = _ret
    return _dict[(k,t)]


def estimate_distribute(_shape):
    def _estimate_distribute(empty_propb, p, N, z):
        initiallstate_copy = deepcopy(ginitiallstate_copy(p))
        ret = deepcopy(empty_propb)

        print ( f" value of p : {p}")        
        hook = False
        
        # hook
        if ret.shape != _shape:
            ret = np.zeros\
                ( shape=_shape ).astype(np.float128)
            hook = True
                
        for k in range(1, ret.shape[0]):
            for t in range(ret.shape[1]):
                ret[k][t] = f(k,t, p=p, N=N,z=z, _dict = initiallstate_copy )
        return normalize_propb(ret).flatten() 
    return _estimate_distribute

def model_fit( propb_mes ):
    empty_propb = np.zeros( shape=propb_mes.shape ).astype(np.float128)
    return  ft.g_extract_coef( empty_propb, propb_mes.flatten(),
                estimate_distribute(empty_propb.shape) , p0=[0.14, 10, 10] )

def plot_propb(propb, single = False, k=2, scat =False, hist = 0, c= None, linewidth=1):
    if not single:
        for i in reversed( range(1,len(propb), 2) ):
            if not scat:
                plt.plot(  i + np.arange(propb.shape[-1]) , 0.1*i+  propb[i] , c=next(gColors) )
            else:
                plt.scatter(  i + np.arange(propb.shape[-1]) , 0.1*i+  propb[i] , c=next(gColors), s=0.1 )
    else:
        if c is None:
            c = next(gColors)
        
        def add_zero( _arr , val = 0):
            return [ val ]  + _arr.tolist()
        
        plt.plot(  add_zero( 6*hist + np.arange(  propb.shape[-1])),
         add_zero( 0.01*hist +  propb[k], val = 0.01*hist ),  c=c , linewidth=linewidth)

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
    plt.plot( At, c=next(gColors) )

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

def load_data_pkl_filter( _hist=40):
    # colonies = pickle.load( open("pkl/colonys-prob_test-2021-06-21_14-27-21.872035.pkl", "rb")) #[:120]
    colonies = pickle.load( open("pkl/colonys-prob_test-2021-06-21_16-55-07.727509.pkl", "rb"))
    good_colonies, _ = filter_colonies(colonies)
    pkt =  histogram_calc(  convert_colonys_to_ktlist(good_colonies), 1, 100, 80)
    return pkt[:50,_hist:]

# def model_fit_hist( ):
#     model_fit( )

# fit the model, and returns pkt.
def lam_fit_pkt(prob ):
    def _lam_fit_pkt( ):
        popt, (pcov, (_range, values)), MSR = model_fit( normalize_propb(deepcopy( prob)))    
        return values.reshape( prob.shape )
    return _lam_fit_pkt

simvsdata = test( [
      lambda : normalize_propb
                    (sim( 45, 100, 5000))[:80*1],
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
                    (sim( 45, 100, 5000 ))[:80*1],
      lambda : lam_fit_pkt
                (normalize_propb
                    (sim( 50, 100, 5000 ))[:80*1]) (),
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


def plot_model_p( p, emptyshape=(6,10) ):
    pkt = estimate_distribute(emptyshape)( np.zeros((1,1)) , p , 1, 1).reshape(emptyshape)
    plot_mean_func_of_time(pkt)
    plot_mean_func_of_time(pkt[:, :40])

def battle( _method, genplot_pks = lambda : 1, genplot_mean = lambda : 1 , name=str(random()) ):
    pkts = _method()    
    genplot_pks()
    
    plt.savefig(f"./svg/pkt-{name}.svg", bbox_inches='tight')    
    plt.clf()

    resetgencolor()
    for pkt in pkts :
        plot_mean_func_of_time( pkt )
    genplot_mean()
    plt.savefig(f"./svg/average-{name}.svg", bbox_inches='tight')
    plt.clf()

def demonstrate_sim():

    colors =[ "red" , "blue" , "black" ]
    Legend = [ ]

    for i, p in reversed(list(enumerate([ 0.1, 0.12, 0.14  ]))):
        pkt = sim(50, 200, 200000, p=p)
        pkt = normalize_propb(pkt)
        c   = next(gColors)
        for j, k in reversed(list(enumerate([3 , 7]))):
            plot_propb( pkt[:,:60],single=True,k=k, hist= i, c=c, linewidth=(j+1))
            Legend.append( f"$ \\psi_{k}  $ where $ p = {p} $ "  )
    Legend = [ r""  +  legend for legend in Legend ]
    plt.legend( Legend )
    plt.xlabel(r"time ticks in $\tau$'s ")
    plt.ylabel(r"$ \psi_k \left( t \right) $")
    plt.title( r"simulation cases")
    plt.savefig("simulation_cases.svg", bbox_inches='tight')
    plt.clf()

def pmodel_relative_psim():
    # _file = open("MSR" , "a+")
    ps = np.linspace( 0.05 , 0.3 , 30 )
    pmodel = [ ]
    for p in ps:
        print(p)
        simpkt = sim(50, 200, 5000, p=p)
        simpkt = normalize_propb(simpkt)
        popt, (pcov, (_range, values)), MSR = model_fit( normalize_propb(deepcopy( simpkt[:,:60])))    
        pmodel.append( popt[0] )
        print(MSR)
        # _file.write(f"{MSR},")
    
    
    popt, (pcov, (_range, values)), MSR = ft.g_extract_coef(np.array(pmodel), ps, lambda x,a,b: a*x + b) 

    plt.plot( _range, values, c=next(gColors) )
    plt.errorbar(pmodel, ps , yerr= pcov[0,0] * np.array(pmodel) + pcov[1,1], fmt='o', c= "black" ,alpha=0.2)
    # plt.scatter( pmodel, ps , c=next(gColors) )
    print( f"pcov:{pcov}")
    print( f"popt:{popt}" )
    # plt.xlabel
    # plt.legend( Legend )
    plt.ylabel(r"$p_{sim}$")
    plt.xlabel(r"$p_{model}$")
    plt.title( r"$p_{model}$ relative to $ p_{simulation} $ ")
    fittedline = f"fitted line $ {popt[0]:.3}x + {popt[1]:.3} $"
    plt.legend( [r"" + fittedline, r"$(p_{model}, p_{sim})$"  ] )
    plt.savefig("pmodelvs_psim.svg", bbox_inches='tight')
    plt.clf()



if __name__ == "__main__":

    # pmodel_relative_psim()
    # demonstrate_sim()
    # exit(0)

    def plot_axes():
        plt.xlabel(r"time ticks in $\tau$'s")
        plt.ylabel(r"$ \mathbf{E} [\psi_{k}(t)]$")

    def plot_title( _title ):
        def lam():
            plt.title(_title)
            plot_axes()
            vars  = deepcopy(_title).split()
            print(vars)
            plt.legend(   [ vars[-1] , vars[-3] ]  )
        return lam
    battle( simvsdata, genplot_pks=plot_title(r"sim vs data"), genplot_mean=plot_title(r"average, data vs sim"), name="sim-vs-data")
    # plt.legend( [ "sim" , "model" ]  )
    battle( simvsmodel,  genplot_pks=plot_title(r"sim vs model"),  genplot_mean=plot_title(r"average, model vs sim"), name="sim-vs-model")
    # plt.legend( [ "model" , "data" ]  )
    battle( datavsmodel, genplot_pks=plot_title(r"data vs model"),  genplot_mean=plot_title(r"average, data vs model"), name="data-vs-model") 
    
    # plot_model_p(0.108)
    # plot_title( r"model for $ p = 0.108 $ ")()
    # plt.show()
    exit(0)



