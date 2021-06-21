
p=0.5

D = { (0,0) : 1 , (1,0) : 1 , (1,1) : p }


import numpy as np
import scipy.special
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np
from os import walk, system
from scipy import fft
import matplotlib
matplotlib.rcParams['text.usetex'] = True
from pickle import dumps, dump, load, loads
import pickle
import cv2 
import sys, os
from numpy.lib import vectorize
sys.path.append(os.path.realpath(".."))
from lib.exp import Exp, FittedObj, gencolor, glabels, plotlabels, gColors, putlabel
import matplotlib
matplotlib.rcParams['text.usetex'] = True

def contrast( _path = "./mgyi/mgyi_1"):
    for dirpath, dirname, filenams in walk(_path): 
        for filename in filter( lambda s : "tif" in s , filenams):   
            print(filename)
            DIRname = dirpath.split("/")[-1]
            system( f"echo y | sh ./scripts/impr.sh {dirpath}/{filename} ./tif3/{DIRname}/{filename}")

def loadf(_path):
    print(_path)
    return cv2.imread(f'{_path}', cv2.IMREAD_GRAYSCALE).astype(np.float128)

def time (_path ):
    return int(_path.split('_')[1])

def norm(X):
    return np.linalg.norm(X, 'fro' )

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


if __name__ == "__main__":


    contrast( _path = "./mgyi_3/" )

    exit(0)

    _keyword = "YFP"

    _dict = None
    with open(f"{_keyword}.pkl", 'rb') as handle:
        _dict = pickle.load(handle)
    # print(_dict)

    propb = np.zeros((50,250))
    # probexcpt = np.zeros((50,250))


    for _key, dictvec in _dict.items():
        print("----" + str( _key) + "-----")
        if len(dictvec) > 0:
            vec = (np.array( sorted( dictvec, key=lambda x: x[0]) )).T
            vec[1] /= vec[1][0] 

            # print(vec.T)

            for timetick, mulsize in vec.T:
                if timetick < 110:
                    propb[int(mulsize)][int(timetick)] += 1
    
    from scipy.signal import find_peaks



    for _key, dictvec in _dict.items():
        print("----" + str( _key) + "-----")
        if len(dictvec) > 0:
            vec = (np.array( sorted( dictvec, key=lambda x: x[0]) )).T
            vec[1] /= vec[1][0] 
            # print(vec.T)

            for timetick, mulsize in vec.T:
                peaks, _ = find_peaks(  propb[int(mulsize)] )
                np.diff(peaks)
                if len(peaks) > 1 :
                    if abs(timetick - peaks[-1]) < 3:
                    # if timetick > 2 *  np.average(propb[int(mulsize)]):
                        print( f" anomaloy :  {_key} , {timetick}, {mulsize}")

    print(propb[7])