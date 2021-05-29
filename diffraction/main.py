
from typing import overload
# from numpy.ma import indices

import sys, os
sys.path.append(os.path.realpath(".."))

import numpy as np 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from lib import fit as ft
from lib.exp import Exp, FittedObj, glabels, plotlabels

import matplotlib
matplotlib.rcParams['text.usetex'] = True


def merge(f, g, x0, x1):
    ret = np.zeros( len(f) )
    ret[:x0] = f[:x0]
    ret[x0:x1] = g[x0:x1]
    ret[x1:] = f[x1:]
    return ret 

_files = [
    "txt/single0.02exp2.txt",
    "txt/single0.04Good2.txt",
    # "txt/single0.04Good.txt",
    "txt/single0.08Good2.txt",
    # "txt/single0.08Good.txt",
    "txt/single0.16Good2.txt",
    # "txt/single0.16Good.txt",
]

_files2 = [ 
    "txt/double0.5and0.08.txt",
    "txt/double0.5and0.04.txt",
    "txt/double0.25and0.04.txt",
    "txt/double0.25and0.08SYM.txt",
    "txt/double0.5and0.08SYM.txt"
]

def find_wave_length(x_vals, y_vals, drop=2):
    _fited = FittedObj(1/x_vals,\
         np.pi *x_vals/ y_vals, lambda x,a,b : a*x, drop=drop)    
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
        _fited = find_wave_length(x_vals, y_vals, drop=None)
        _fited.plot()
        glabels.send(r"single slit, $\lambda = {0}$".format(np.around( _fited.func_parmas[0][0],2)))
        # plt.legend( )
        print(_fited.func_parmas[0], _fited.func_parmas[1][0] )
    
    def fun(self, x, i):
        def f(x, a, b, c, d):
            return ( b* np.sin(a*x+c )/(a*x+c ) )**2
        return f(x, *self.func_parmas[i][0])
    def plot(self):
        _range  = np.linspace( np.min(self.data[0][0]),\
             np.max(self.data[0][0]), 1000 )
        indices = abs(_range) < 0.5        
        _range = _range[indices] 
        for i, _ in enumerate(self.data):
            indices = abs(self.data[i][0]) < 0.5
            plt.scatter( self.data[i][0][indices], self.data[i][1][indices], s=0.5)
            plt.plot( _range, self.fun(_range, i) )
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
            return [y , x]
        
        super().__init__(_files2, _filter=_filter, _metaformat=_metaformat)

    def extract_func(self):
        super().extract_func(ft.extract_coef2)
    
    def find_wave_length(self):
        y_vals = np.array([[_parmas[0][0], _parmas[0][-1]] for _parmas in self.func_parmas])
        x_vals = np.array(self.meta)
        x_vals, y_vals  = x_vals.transpose(), y_vals.transpose()
        print(x_vals)
        print(y_vals)
        for i in range(2):
            if i == 1:            
                _fited = find_wave_length( ( x_vals[1])/2, y_vals[i], drop=None)
                glabels.send(r"double slit, $\lambda(d) = {0}$".format( np.around(_fited.func_parmas[0][0],2)))

            else:
                _fited = find_wave_length(x_vals[i], y_vals[i], drop=None)
                glabels.send(r"double slit, $\lambda(W) = {0}$".format(np.around(_fited.func_parmas[0][0],2)))

            _fited.plot()
            print("res:")
            print(_fited.func_parmas[0], _fited.func_parmas[1][0] )

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
    plt.title(r' $\frac{1}{\alpha} \ - \ W \ linear \ connection $')
    plt.xlabel(r' $ (\frac{1}{W} \ or \ d)_{[mm]} $')
    plt.ylabel(r' $ \frac{1}{\alpha(W)}_{[\cdot]}$ ')
    plotlabels()
    plt.savefig("./graphs/Fig1.svg")
    plt.clf()
    _singleExp.plot()
    plt.show()