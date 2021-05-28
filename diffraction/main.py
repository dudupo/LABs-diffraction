

import sys, os
sys.path.append(os.path.realpath(".."))

import numpy as np 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from lib import fit as ft
from lib.exp import Exp, FittedObj

def merge(f, g, x0, x1):
    ret = np.zeros( len(f) )
    ret[:x0] = f[:x0]
    ret[x0:x1] = g[x0:x1]
    ret[x1:] = f[x1:]
    return ret 

# _files = [
#         "single0.02Bblack.txt",  "single0.04.txt",
#         "single0.02no.txt",      "single0.08.txt",
#         "single0.02B.txt",  "single0.02no2.txt",     "single0.08a.txt" ]

_files = [
    "txt/single0.02exp2.txt",
    "txt/single0.04Good2.txt",
    "txt/single0.04Good.txt",
    "txt/single0.08Good2.txt",
    "txt/single0.08Good.txt",
    "txt/single0.16Good2.txt",
    "txt/single0.16Good.txt",
]

_files2 = [ 
    "txt/double0.5and0.8.txt",
    "txt/double0.5and0.4.txt",
    "txt/double0.25and0.8.txt",
    "txt/double0.5and0.8SYM.txt"
]


    # def plot(self, ) 

# LENGTH = 0.3


def find_wave_length(x_vals, y_vals, drop=2):
    _fited = FittedObj(x_vals,\
         np.pi *x_vals/ y_vals, lambda x,a,b : a*x+b, drop=drop)    
    return _fited


class SingleExp (Exp):
    
    def __init__(self):
        
        def _filter( _matrix ):
            X, Y = _matrix[-1], -_matrix[-2]
            ind = np.abs(_matrix[-1]) < 1 
            X, Y = _matrix[-1][ind], -_matrix[-2][ind]
            return np.array([ X, Y])
        
        def _metaformat(_str):
            def find_width_length(self):
                    return float(_str[10:14])
            return find_width_length(_str)
        super().__init__(_files, _filter=_filter, _metaformat=_metaformat)


    def extract_func(self):
        super().extract_func(ft.extract_coef)

    def find_wave_length(self):
        y_vals = np.array([_parmas[0][0] for _parmas in self.func_parmas])
        x_vals = np.array(self.meta)
        print(x_vals, y_vals)
        _fited = find_wave_length(x_vals, y_vals)
        print(_fited.func_parmas[0], _fited.func_parmas[1][0] )
        _fited.plot()
        plt.show()
        # plt.scatter(x_vals, y_vals, marker='o')
        # plt.plot(_range, values)
        # plt.show()
        # print(MSR)
        # print(np.pi/ (popt[0]))

class doubleSlit(Exp):

    def __init__(self):
        
        def _filter(_matrix):
            X, Y = _matrix[-1], -_matrix[-2]
            ind = np.abs(_matrix[-1]) < 1 
            X, Y = _matrix[-1][ind], -_matrix[-2][ind]
            ind = (Y > 0) & ( np.abs(X) < 0.5)  
            X, Y = X[ind], Y[ind]

            return np.array([ X, Y])
            
        def _metaformat(_str):
            _str = _str[10:]
            findex = _str.index('a')
            x = float(_str[:findex])
            _str = _str[findex+3:]
            findex = _str.index('SYM') if 'SYM' in _str else _str.index('.txt') 
            y = float(_str[:findex])
            return [x , y]
        
        super().__init__(_files2, _filter=_filter, _metaformat=_metaformat)

    def extract_func(self):
        super().extract_func(ft.extract_coef2)
    
    def find_wave_length(self):
        # print( self.meta)
        y_vals = np.array([[_parmas[0][0], _parmas[0][-1]] for _parmas in self.func_parmas]).reshape(1,-1)
        x_vals = np.array(self.meta).reshape(1,-1)
        x_vals, y_vals  = x_vals.transpose(), y_vals.transpose()
        print(x_vals)
        print(y_vals)
        
        _fited = find_wave_length(x_vals[1], y_vals[0], drop=1)
        print(_fited.func_parmas[0], _fited.func_parmas[1][0] )
        _fited.plot()
        plt.show()


def genLegend( _file ):
    return _file[_file.index("."):_file.rindex(".")]

if __name__ == "__main__":

    _singleExp = SingleExp()
    _singleExp.extract_func()
    _singleExp.find_wave_length()
    PPP = doubleSlit()
    PPP.extract_func()
    
    print(PPP.meta)
    PPP.find_wave_length()