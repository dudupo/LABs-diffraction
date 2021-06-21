
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
# from lib import fit as ft
# from lib.exp import Exp, FittedObj, gencolor, glabels, plotlabels, gColors, putlabel
# from copy import deepcopy
#
# from bio.mc import get_multi_factor
# from bio.sim import sim, exp
# from bio.utility import histogram_calc, convert_colonys_to_ktlist
import pickle
from datetime import datetime

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


    # plt.xlabel(r"time ticks")
    # plt.ylabel(r"$ \sum{ \frac{\psi_i}{\psi_0 } } \ge k  $")


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
    propb = propb.T
    At = np.zeros( propb.shape[-1] )
    for j, vec in enumerate( propb ):
        At += vec * j 
    plt.scatter(  np.arange( len(At[2:] )) , At[2:] )
    plt.show()

def simvsmodel( ):
    simpkt = sim( 60, 100 )[:60]
    plot_propb (  normalize_propb(deepcopy(simpkt)) ) 
    popt, (pcov, (_range, values)), MSR = model_fit( normalize_propb(deepcopy( simpkt)))
    plot_propb( values.reshape( simpkt.shape ) )
    plt.show()

def simvsdata( prob ):


    simpkt = normalize_propb(sim( 60, 100 ))[:80*1]
    
    print(simpkt.shape)
        
    r = next(gColors)
    plot_propb (  simpkt  )
    
    t = next(gColors)
    while t != r :
        t = next(gColors)

    plot_propb ( normalize_propb(prob[:80*1])[:,40:] )
    plt.show()

def datavsmodel( _file ): 
    prob = pickle.load( open(_file, "rb"))[:20, 15:40] 
    plot_propb ( normalize_propb(deepcopy(prob)))
    popt, (pcov, (_range, values)), MSR = model_fit( normalize_propb(deepcopy( prob)))
    plot_propb( values.reshape( prob.shape ) )
    plt.show()

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




if __name__ == "__main__":
    

    def get_multi_factor_pckl(t, col_idx=1):
        k = 1
        colonies = pickle.load( open("colonys-prob_test-2021-06-21_14-27-21.872035.pkl", "rb"))
        pkt = np.zeros((k, t))
        col_to_plot = []
        col_0 = 0
        for i, colony in enumerate(colonies):
            colony_zero = colony[0].pix_num
            start_idx = 0
            found = False
            for j, col_t in enumerate(colony):
                if col_t.pix_num // colony[0].pix_num < 2:
                    continue
                elif not found:
                    colony_zero = col_t.pix_num
                    found = True
                    start_idx = j
                if j < t and found:
                    if i == col_idx:
                        col_to_plot.append(col_t.pix_num/ colony[0].pix_num)
                    cur_k = int(col_t.pix_num // colony_zero)
                    if cur_k <= k:
                        pkt[cur_k - 1][j - start_idx] += 1
                    else:
                        pkt = np.r_[pkt, np.zeros((cur_k - k, t))]
                        k = cur_k
                        pkt[cur_k - 1][j - start_idx] += 1

        if len(col_to_plot) > 0:
            plt.plot(col_to_plot)
            plt.plot(exp(40 , 0.125))
            plt.show()
        else:
            print("single colony plot failed, try different colony")
        return pkt


<<<<<<< HEAD
    colonies = pickle.load( open("pkl/colonys-prob_test-2021-06-21_14-27-21.872035.pkl", "rb"))[:120]
    pkt = histogram_calc(  convert_colonys_to_ktlist(colonies), 1, 100, 80)
    simvsdata(pkt)
=======
    # picklize()
    colonies = pickle.load( open("./pkl/colonys-prob_test-2021-06-21_14-27-21.872035.pkl", "rb"))
    for i in range(10):
        colony_k(colonies, i, True)
    plt.show()
>>>>>>> f19c179e641fd4a7b0afa0f8310e1788d7c4eb81
    exit(0)

    # colonies = pickle.load( open("./pkl/colonys-prob_test-2021-06-19_19-17-09.429712.pkl", "rb"))[:120]
    # pkt = histogram_calc(  convert_colonys_to_ktlist(colonies), 10, 100, 60)
    # simvsdata(pkt)
    # exit(0)

    # pkt = get_multi_factor_pckl( 50 , col_idx=4)

    # simvsmodel( )
    # simvsdata("prob_test-2021-06-12 16:41:23.299589.pkl")
    # datavsmodel("prob_test-2021-06-12 16:41:23.299589.pkl" )
    # plot_mean_func_of_time()
    #     plot_propb ( deepcopy( prob ) )
# 

    # prob = pickle.load( open( "prob_test-2021-06-07 20:34:17.759904.pkl" , "rb"))
    # prob = pickle.load( open("prob_test-2021-06-12 15:35:52.515836.pkl", "rb"))
    # prob = pickle.load( open("prob_test-2021-06-12 14:19:33.846149.pkl", "rb"))
    


    prob = pickle.load( open("./pkl/prob_test-2021-06-12 16:05:21.825747.pkl".replace(":", "-").replace(" ", "_"), "rb"))
    # prob = pickle.load( open("prob_test-2021-06-12 16:41:23.299589.pkl", "rb"))
    # prob = prob[20:,15:40]

    plot_mean_func_of_time(normalize_propb(deepcopy( prob.T )))

    exit(0)
    # simprob = sim( 20, 100, p = 0.125, number_of_exp=5000  )
    # # plot_propb ( deepcopy( prob )) 

    # plot_propb (  simprob , single=True, k=3 )

    plot_propb (  normalize_propb(deepcopy(simprob)) , single=True, k=2 )
    plot_propb (  normalize_propb(deepcopy(prob)) , single=True, k=2 )

    # plot_propb ( normalize_propb(deepcopy( prob ) ))
    plt.show()

    # exit(0)

    # plt.show()
    # def fitted_for_p(p):
    
    #     print(prob.shape)
    #     return popt[0]
    #     # return [popt, MSR] 
    # popt, (pcov, (_range, values)), MSR = model_fit( normalize_propb(deepcopy( prob[0:180])))
    # plot_propb( values.reshape( prob[0:180].shape ) )
    # plt.show()
    plot_mean_func_of_time(normalize_propb(deepcopy( prob[:100].T )))
    exit(0)
    _prob = normalize_propb(prob[:10])
    simprob =  sim(40, 500, p=0.44, number_of_exp=50000)
    # plot_propb ( deepcopy( values.reshape(  prob.shape) ) )


    # _dict = deepcopy( initiallstate )
    # prob_r = np.zeros(  shape=(40,60) )
    # for k in range(40):
    #     for t in range(60):
    #         prob_r[k][t] = f( k,t, p = 0.4, _dict=_dict)
    # normalize_propb(prob_r)
    # plot_propb( _prob )
    plot_propb( normalize_propb(simprob) )
    plt.show()
    # p_axis = np.linspace( 0.3, 1, 10 )
    # res_p_axis = np.vectorize(fitted_for_p ) (p_axis )
    # plt.plot(p_axis, res_p_axis)
    # plt.show()
        
    # plt.savefig("./svg/fig5kt-fitted.svg")
    # plt.clf()
    # plt.savefig("./svg/fig4kt-sim.svg")
    # plt.clf()
    # points =list( map( np.max , _prob ))

    # def find_i(j, points):
    #     eps = 0.1
    #     for _i ,_j in enumerate( points):
    #         if abs(i - _j) < eps :
    #             return _i, _j

    # for i in range(2, 7):
    #     ( _time1, _pop1 ) , ( _time2, _pop2 )  = find_i(i), find_i(i*2)
        
    #     plt.plot(  [_time1, _time2 , ])
        

        
    #     for _j in points:
    #         if abs(i - _j) < eps :
                # print(_j)    
        # j = points.find( i)
        # print(j)
        # points[2*i]

    # plot_propb( _prob )
    # # plt.plot(  np.arange( len(points)) * points )
    # plt.show()