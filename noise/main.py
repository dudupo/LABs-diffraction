
import sys, os
sys.path.append(os.path.realpath(".."))

import matplotlib
matplotlib.rcParams['text.usetex'] = True


from numpy.lib.function_base import append
import pandas as pd 
import numpy as np 
from matplotlib import pyplot as plt 

from datetime import datetime
from scipy.ndimage.filters import convolve  
from scipy.ndimage import gaussian_filter
from scipy.signal import find_peaks
from scipy import signal
import scipy.integrate as integrate

from lib.exp import Exp

mi = 10**-3
ki = 10 ** 3
stage1 =np.array(  [7, 8, 9,  10,  11,  12,  13,  14,  15,  20,  30,  40,  60,  80, 100]) * ki
stage2 = np.array( list(range(0,4)) +  list( _ * mi for _ in range(300,1000,100)) + [ 1 ] + [  100  ] )

def extract_ramps( _array):
    
    _len = len(_array)
    diff = 100* convolve( _array , np.array([1,-1])) 
    indecis = np.abs(diff) < 0.05
    Q = diff.copy()   
    Q[indecis] = 0 
    
    sign = np.sign( np.average( _array ) )             
    peaks2, _ = find_peaks(Q * sign, width=18)
    realpeaks = []
    peaks = np.array( range(_len))
    A = peaks[peaks2]
    
    B = []
    for j in range(len(A)-1):
        B.append( int( (A[j] + A[j+1])/2 ) )
    B.append(int( (A[-1] + _len)/2 ))
    
    realpeaks = [ False ] * len(_array)
    for b in B:   
        realpeaks[b] = True
    
    return np.array(B)[1:], np.array( realpeaks ) 


class NoiseExp(Exp):

    def __init__(self, _filespath):

        def _filter(dataset):
            ret = {}
            ret["pure"] = dataset.values.transpose()[-1].astype(np.float64)         
            ret["gaus"] = gaussian_filter( ret["pure"] , sigma=10)
            peaks, boolarray = extract_ramps(ret["gaus"])
            ret["peax"] = peaks
            ret["peak"] = ret["gaus"][boolarray]
            return ret

        super().__init__(_filespath, _filter=_filter, _metaformat=None)

        self.peaks = []
        for _data in self.data:
            self.peaks = _data["peak"].tolist() + self.peaks 

FREQ =  np.array(stage2.tolist() + stage1.tolist())
Votage =  2

#66.5 * mi

RESISTANCE = 1.003 * 10 **6
RESISTANCE2 = 1.76 * 10 **6   # k <- 8.514540894565726e-24

def main():
    
    vec = np.array([ 0.1, 0.3, 0.5, 0.6, 0.5, 0.3, 0.1 ]) / 2.4
    
    noiseexp = NoiseExp([   "csv/Trace 2021-05-04 12-54-38 0 copy.csv",
                            "csv/100mili-1hz-1khz 2021-05-04 13-31-12 0 copy.csv"])
    

    PEAKS = np.array(noiseexp.peaks)
    print(PEAKS)
    print(len(FREQ) , len(PEAKS))
    PEAKS /= (2 * 10 **-5)
    plt.scatter( FREQ, PEAKS, c="black", s=5)
    plt.plot(FREQ, PEAKS, c="teal")
    plt.title( r"the gain as function of the frequency")
    plt.xlabel(r"frequency $ \nu_{[Hz]} $")
    plt.ylabel(r"$g(\nu)$")
    plt.legend([r"measured points", r"continuously line"])
    plt.savefig("doc/noise3.svg")
    # result = integrate.quad(lambda x: special.jv(2.5,x), 0, 4.5)
    plt.clf()
    print( "GV : {0}".format(PEAKS))
    # plt.show()

    pC = 2* 63 * 10**-12
    
    for T in [296]:
    
        def fun( gv, R, freq ):
            return (gv **2) / (1 + (2*np.pi*freq*pC*R)**2)

        r_arr = np.array([100*(10**3), 200*(10**3), 330*(10**3), 620*(10**3), 820*(10**3), (10**6), 2*(10**6), 3*(10**6),
                    4.7*(10**6), 6.8*(10**6), 8.2*(10**6), 12*(10**6)])

        R_arr = np.array([97.2*(10**3), 198.1*(10**3), 324*(10**3), 619*(10**3), 819.5*(10**3), 1.02*(10**6), 2.09*(10**6),
                    3.2*(10**6), 4.96*(10**6), 6.77*(10**6), 8.48*(10**6), 12.28*(10**6)])

        # r_arr = np.array([100*(10**3), 200*(10**3), 330*(10**3), 620*(10**3), 820*(10**3), 106, 2*10**6, 3*10**6,
                    # 4.7*10**6, 6.8*10**6, 8.2*10**6, 12*10**6])
        v_arr = np.array([0.3, 0.43, 0.52, 0.74, 0.83, 0.91, 1.23, 1.47, 1.64, 1.84, 2.05, 2.1])


        v_arr_2 = [0.89, 1.2, 1.59, 2.18, 2.5, 2.76, 4, 5.1, 6.31, 8.35, 8.7    ]
        # [ 15moh -> 9.3 v]
        # [100*ki, 200*ki, 




        from scipy.optimize import curve_fit

        def extract_coef( time, distance ):
            def f(x, a, b ):
                return a*x + b
            popt, pcov = curve_fit(f, time, distance, p0 = [ 1.38 * 10**-23 , 0] )
            _range = np.linspace( np.min(time) , np.max(time), 100 )
            return popt, (_range, f(_range, *popt)) 


        vals = []
        yerr = []
        # kvals = []
        for r in R_arr: #R_arr:
            G1, G2, G3 = 0, 0, 0
            for j in range(len(FREQ)-1):
                G1 += (FREQ[j+1] - FREQ[j]) * fun( PEAKS[j], r, FREQ[j]) 
                G2 += (FREQ[j+1] - FREQ[j]) * fun(PEAKS[j],r, FREQ[j])
                G3 += (FREQ[j+1] - FREQ[j]) * fun(PEAKS[j] ** 0.5 ,r, FREQ[j])
            print( 4 * G1 * r )
            yerr.append( G3 *T *r * 4 )
            
            G1 = 4 * G1 * r *T
            # G3 = 4 * G3 * r *T

            # yerr.append( G3 )
            # G1 , G2 = 4 * G1 * r *T, 4 * G2 * r *T
            vals.append( G1)
        
        # print(vals)
        # print(vals[0]**0.5)
        # # Votage /= (2 * 10 **-5)
        # print(Votage**2/ vals[-1])
        vals = np.array(vals)

    # Votage**2 / np.pi * vals[-1]
    #     indecis = (v_arr**2 )< 0.8
    #     print(v_arr)
        indecis = [ True ] * len(vals)
    # # indecis[0] = True
    # # indecis[1] = True 
        print( yerr)
        # kvals = v_arr[indecis]**2/ vals
        yerr =  np.array(yerr) 
        # plt.errorbar( list(range(len(kvals))), kvals  , yerr=yerr , fmt='o', c= "black" ,alpha=0.2)
        # plt.show()
        # plt.clf()        
    #     print(vals)
    #     print(vals)
    #     vals = np.array(vals)
        # yerr = yerr*v_arr[indecis]**2 / (vals[indecis] )
        print(yerr)
        plt.title( r"$\langle V^2 \rangle$ as measured vs $4TRG$" )
        plt.xlabel(r"$4TRG_{ [\frac{K}{J}] } $  ")
        plt.ylabel(r"$\langle V^2 \rangle_{ [V^2] }$ ")
        plt.scatter(vals[indecis], v_arr[indecis]**2 , c= "black" , s=2)
        popt, (X, Y) = extract_coef(vals[indecis], v_arr[indecis] **2 )

    # # # plt.plot(vals[indecis], v_arr[indecis] , "ob")
    # # # popt, (X, Y) = extract_coef(vals[indecis], v_arr[indecis] **2)

        plt.plot(X, Y,  c ="teal")
        plt.plot(X, 1.38*10**-23 *X +popt[1],  c ="darkslategray")
        plt.errorbar( vals[indecis] , popt[0] *vals[indecis] + popt[1], yerr= 0.1*(popt[0] *vals[indecis] + popt[1]) , fmt='o', c= "black" ,alpha=0.2)
        
        plt.legend([ "fitted line", "excpected line", "measured points" ])
        
        print("-- bolzman :  -- ")
        print( popt)
        # plt.show()
        plt.savefig("doc/noise1.svg")
        plt.clf()
    # print( (G1 + G2)/2 )


def main2():
    # 10 khz
    datasets = [ pd.read_csv(_file) for _file in [ "csv/Icol 2021-05-20 16-29-23 0.csv"] ] 
    


    ret = None

    for dataset in datasets:  
        V =  dataset.values.transpose()[-1].astype(np.float64)
        U = gaussian_filter( V, sigma=1)
        plt.plot( U )
        # U = gaussian_filter( U, sigma=2)
        Q = U[10:] - U[:-10] #np.abs(U[10:] - U[:-10]) > 0.4
        peaks2, _ = find_peaks(Q , distance=400)
        # plt.plot(U)
        # plt.plot(Q)
        # plt.scatter( Q, U[Q] )
        # peaks2 = peaks2[:-10] 
        plt.plot(np.arange(len(U)-10)[peaks2],  U[10:][peaks2], marker='o')
        
        ret = U[10:][peaks2]
        
        # plt.plot( V)
        plt.show()

    Vin = np.array([20, 40, 60, 80, 120, 200, 240, 280] ) * 10**-3
    
    # plt.plot(Vin,  ret)
    # plt.show()

from Thermocouple import tc
def main3():
    datasets = [ pd.read_csv(_file) for _file in [  "csv/Vol 2021-05-18 12-10-25 0.csv",
     "csv/Temp 2021-05-18 12-09-15 0 2021-05-18 12-09-51 0.csv" ] ] 
    
    final_data = []
    indecis = [ ]
    for _ , dataset in enumerate( datasets):
        V =  dataset.values.transpose()[-1].astype(np.float64)
        _time = dataset.values.transpose()[-2].astype(np.datetime64)
        if _ == 0:
            indecis = V > 0.9
            
            V = V[indecis]
            _time = _time[indecis]

        if ( _ == 1):
            indecis  = V*10**3 > -5.891           
            V = V[indecis]
            _time = _time[indecis]
        
            # plt.plot(V)
            # plt.show()
            
            V =np.array(list(map( lambda x :  tc.Thermocouple.mv_to_typek( x * (10**3) ) + 273.15 , V)))

        if ( _ == 0):
            indecis = np.random.choice(V.shape[0], 1000, replace=False) 
            V = V[indecis]
            _time = _time[indecis]
        
        
        final_data.append((_time, V))
    
    X, Y = [], []
    for (_time_, v) in zip(final_data[0][0],final_data[0][1]) :
        j = np.argmin(  np.array(list(map(lambda x:  abs(x-_time_ ), final_data[1][0]))), 0)   
        X.append( final_data[1][1][j] )
        Y.append(v**2)
    plt.scatter( X,Y, s=0.8, c="black" )

    
    from scipy.optimize import curve_fit

    def extract_coef( time, distance ):
        def f(x, a, b ):
            return a*x + b
        popt, pcov = curve_fit(f, time, distance )
        _range = np.linspace( np.min(time) , np.max(time), 100 )
        return popt, (_range, f(_range, *popt)) 

    popt, (_X,  _Y) = extract_coef( X,Y)
    plt.plot(_X, _Y, c="teal")
    plt.title(r"$\langle V^2 \rangle$ as measured vs $ T $")
    plt.xlabel(r"$T_[K]$")
    plt.ylabel(r"$\langle V^2 \rangle_{[V^2]}$")
    print(popt) 
    plt.savefig("doc/noise2.svg")


if __name__ == "__main__":
    main()
    # main3()
    # main2()
