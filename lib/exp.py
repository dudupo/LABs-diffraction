
from . import fit as ft 
import matplotlib.pyplot as plt
import numpy as np

class FittedObj:
    def __init__(self, time, distance, func, drop =None):
        self.time, self.distance, self.func = time, distance, func
        if drop is not None:
            self.func_parmas = ft.g_extract_coef_drop(time, distance, func, drop_points=drop)
        else:
            self.func_parmas = ft.g_extract_coef(time, distance, func) 
        
        self.funclam = lambda  x : func( x, *self.func_parmas[0])
        # self.yerr = g_derivate_err( time, distance, self.funclam)  

    def plot( self ):
        plt.scatter( self.time, self.distance) 
        plt.plot( self.time, self.funclam(self.time) )
        plt.show()


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
    
