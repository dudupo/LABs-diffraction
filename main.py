


import numpy as np 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def g_extract_coef( time, distance, f ,p0=None):
    popt, pcov = curve_fit(f, time, distance, p0=p0)
    _range = np.linspace( np.min(time) , np.max(time), 1000 )
    values = f(_range, *popt)
    val_in_given_range = f(time, *popt)
    MSR = np.sqrt(np.sum((val_in_given_range - distance)**2))/ len(values)
    # print(MSR)
    return popt, (pcov, (_range, values)) , MSR 

from itertools import combinations
def g_extract_coef_drop( time, distance, f, drop_points ):
    _length = len(time) - drop_points
    indices = np.arange( len(time) )

    time, distance = np.array(time), np.array(distance) 

    return min ([ g_extract_coef( np.take(time, _indices, 0), np.take(distance, _indices, 0), f )\
         for _indices in  combinations( indices, _length)], key= lambda ret : ret[-1])

def extract_coef( time, distance ):
    def f(x, a, b, c, d):
        return ( b* np.sin(a*x )/(a*x ) )**2  
    return g_extract_coef(time, distance, f)

def extract_coef2( time, distance ):
    def f(x, a, b, c, a2): # c2, d, e, e2):
        return ( b* np.sin(a*x+c)/(a*x+c) *  np.cos(a2*x))**2 
        # return ( b* np.sin(a*np.sin(x) + c)/(a*x + d) *  np.cos(a2*np.sin(x) +c2 ))**2   
    
    return g_extract_coef(time, distance, f, p0=[ 10, 1, 0.5, 200  ] )  

    


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
    "single0.02exp2.txt",
    "single0.04Good2.txt",
    "single0.04Good.txt",
    "single0.08Good2.txt",
    "single0.08Good.txt",
    "single0.16Good2.txt",
    "single0.16Good.txt",
]

_files2 = [ 
    "double0.5and0.8.txt",
    "double0.5and0.4.txt",
    "double0.25and0.8.txt",
    "double0.5and0.8SYM.txt"
]



class Exp:

    def __init__(self, _filespath, _filter = None, _metaformat= None ):
        self._files = _filespath
        self.data = []
        self.meta = []

        for _file in self._files:
            tempdata = np.array( [  np.fromstring(line, sep=" ")  for line in  open( _file, "r").readlines() ])
            tempdata =  tempdata.transpose()
            if  _filter is not None: 
                tempdata = _filter(tempdata)
            self.data.append( tempdata )

            if _metaformat is not None:
                self.meta.append(  _metaformat( _file ) )
    
    def extract_func(self, extract_coef_func):
        self.func_parmas =  [ [popt, pcov, MSR] for popt , (pcov, points), MSR in\
             map( lambda _matrix :  extract_coef_func(_matrix[0], _matrix[1]) , self.data)] 
    

class SingleExp (Exp):
    
    def __init__(self):
        
        def _filter( _matrix ):
            X, Y = _matrix[-1], -_matrix[-2]
            ind = np.abs(_matrix[-1]) < 1 
            X, Y = _matrix[-1][ind], -_matrix[-2][ind]
            return np.array([ X, Y])
        
        def _metaformat(_str):
            return float(_str[6:10])

        super().__init__(_files, _filter=_filter, _metaformat=_metaformat)


    def extract_func(self):
        super().extract_func(extract_coef)

    def find_width(self):
        def linear_func( x, a, b):
            return a*x + b
        y_vals = [_parmas[0][0] for _parmas in self.func_parmas]
        x_vals = self.meta
        popt, (pcov, (_range, values)) , MSR =  g_extract_coef_drop(x_vals, y_vals, linear_func, 2)
        plt.scatter(x_vals, y_vals, marker='o')
        plt.plot(_range, values)
        plt.show()
        print(MSR)
        print(np.pi/ (popt[0]))

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
            _str = _str[6:]
            findex = _str.index('a')
            x = float(_str[:findex])
            _str = _str[findex+3:]
            findex = _str.index('SYM') if 'SYM' in _str else _str.index('.txt') 
            y = float(_str[:findex])
            return x , y
        
        super().__init__(_files2, _filter=_filter, _metaformat=_metaformat)

    def extract_func(self):
        super().extract_func(extract_coef2)

def genLegend( _file ):
    return _file[_file.index("."):_file.rindex(".")]

if __name__ == "__main__":


    G = []

    _singleExp = SingleExp()
    _singleExp.extract_func()
    # _singleExp.find_width()
    # for data in _singleExp.data:
    #     plt.plot( data[0], data[1])  
    # plt.show()
    # print(QQQ.meta)
    # # print(QQQ.data)
    PPP = doubleSlit()
    PPP.extract_func()
    # print( PPP.meta)
    # exit(0)
    exit(0)

    legends = []
    for i , (__files, g) in enumerate( zip([_files, _files2],  [extract_coef, extract_coef2])):
        for _file in __files:

            A = np.array( [  np.fromstring(line, sep=" ")  for line in  open( _file, "r").readlines() ])
            A = A.transpose()
            
            print(A.shape)
        
            X, Y = A[-1], -A[-2]
            ind = np.abs(A[-1]) < 1 
            X, Y = A[-1][ind], -A[-2][ind]

                
            # Y = Y / np.sqrt( np.sum( Y**2 ) )
            
            if i == 1:
                ind = (Y > 0) & ( np.abs(X) < 0.5)  
                X, Y = X[ind], Y[ind]
            
            plt.plot(X , Y )
            # plt.plot( A[-1], -Y)
            legends.append(genLegend(_file) )
            
            try:
                pass

                __ , (_, ( X0, Y0 )), MSR = g(X, Y)
                print(__)
                plt.plot(X0,Y0)
                legends.append( genLegend( _file ) + " fit"  )
                if i == 0:
                   G.append( __[0])
            except Exception as e:
                print("error {0}".format(e))
                pass
            plt.legend( legends )
            # plt.show()
            plt.savefig( _file +".svg" )
            legends = []
            plt.clf()

    plt.plot( list(range(len(G))) ,G)
    plt.show()