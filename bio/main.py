
p=0.5

D = { (0,0) : 1 , (1,0) : 1 , (1,1) : p }



from random import choice, gammavariate
import numpy as np
from copy import deepcopy
def exp():
    x = 1
    t = 0 
    L = []
    while x < 100:
        for _ in range(x):
            x = choice([x+1 , x ])
        L.append(deepcopy(x))
        t += 1
    return L




import scipy.special
import scipy.integrate as integrate





def f(k , t, _dict =D):
    if (t == 0 and k > 0) or (t< 0) or (k<0) :   
        return 0
    if (k,t) in _dict:
        return _dict[(k,t)]
    else:
        _ret = 0 
        for j in range(int(k/2)+1):
            _ret  += scipy.special.binom( int(k/2) , j) *  p**j * (1-p)**( int(k/2)  + 1 - j ) * f(k-j,t-1, _dict)
    _dict[(k,t)] = _ret
    return _dict[(k,t)]

import matplotlib.pyplot as plt
import numpy as np



def sim_main():
    X = [f(100, i) for i in range(100)]
    # for i in range(len(X)):
        # X[i] /= sum(X)
    # for i in range(len(X)):
        # X[i] /= sum(X)
    X = np.array(X)
    # plt.plot( X / np.max( X) )
    print( sorted(D.items(),  key= lambda x : x[1]) )

    # # result = integrate.quad( X , 0, 4.5)
    # plt.show()
    plt.plot(exp())
    plt.yscale("log")
    plt.show()



    # X = np.array([ exp( ) for _ in range(2000)]
    
    
    
    # plt.scatter( )
    from numpy import histogram
    T = histogram( X  , bins=20)
    print(T)
    # plt.plot(T[1])
    plt.scatter( T[1][1:] , T[0]/ np.max(T[0]), s=2) 
    plt.show()


from os import walk, system

def contrast():
    _path = "./mgyi/mgyi_1"
    for dirpath, dirname, filenams in walk(_path): 
        for filename in filter( lambda s : "tif" in s , filenams):   
            print(filename)
            DIRname = dirpath.split("/")[-1]
            system( f"echo y | sh ./impr.sh {dirpath}/{filename} ./tif/{DIRname}/{filename}")

import cv2 
def loadf(_path):
    print(_path)
    return cv2.imread(f'{_path}', cv2.IMREAD_GRAYSCALE).astype(np.float128)

def time (_path ):
    return int(_path.split('_')[1])

def norm(X):
    return np.linalg.norm(X, 'fro' )


from scipy import fft
import matplotlib
matplotlib.rcParams['text.usetex'] = True

from copy import deepcopy

from pickle import dumps, dump, load, loads
import pickle


def picklize():
    _path = "./tif/"

    _keyword = "YFP" #"Phasefast"


    _dict = {}
    for dirpath, dirname, filenams in walk(_path): 
        DIRname = dirpath.split("/")[-1]
        _dict[DIRname] = [  ]
        for _filename in filter( lambda s : ("tif" in s ) and ( _keyword in s), filenams):               
            __file = "{0}/{1}".format( dirpath, _filename )
            print(__file)
            _dict[DIRname].append( ( time(__file), norm(loadf( __file ))))
        print(_dict)

    with open(f"{_keyword}.pkl", 'wb') as handle:
        dump(_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


import sys, os

from numpy.lib import vectorize
sys.path.append(os.path.realpath(".."))

from lib.exp import Exp, FittedObj, gencolor, glabels, plotlabels, gColors, putlabel

import matplotlib
matplotlib.rcParams['text.usetex'] = True

from skimage.measure import regionprops







if __name__ == "__main__":
    _keyword = "YFP"

    _dict = None
    with open(f"{_keyword}.pkl", 'rb') as handle:
        _dict = pickle.load(handle)

    propb = np.zeros((50,250))

    for _key, dictvec in _dict.items():
        print("----" + str( _key) + "-----")
        if len(dictvec) > 0:
            vec = (np.array( sorted( dictvec, key=lambda x: x[0]) )).T
            vec[1] /= vec[1][0] 
            # print(vec.T)

            for timetick, mulsize in vec.T:
                propb[int(mulsize)][int(timetick)] += 1
                print(timetick, mulsize, _key)



    

    # for i in reversed( range(2, 15) ) :
    #     plt.plot(  50*i + np.arange(250) , 5*i + propb[i] , c=next(gColors) )
    # plt.xlabel(r"time ticks")
    # plt.ylabel(r"$ \sum{ \frac{\psi_i}{\psi_0 } } \ge k  $")
    # plt.savefig("./svg/fig1.svg")