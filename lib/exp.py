
from numpy.core.numeric import indices
from . import fit as ft 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


COLORS = [ "teal", "olive", "darkslategray" ]

def gencolor ():
    i = 0
    while True:
        yield COLORS[ i ]
        i += 1
        i %= len(COLORS)


gColors = gencolor()

LABELS = [ ]

def putlabel():
    while True:
        _str = yield
        LABELS.append(_str)

glabels = putlabel()
next(glabels)
def plotlabels():
    global LABELS
    plt.legend(LABELS)
    LABELS = []


class FittedObj:
    def __init__(self, time, distance, func, drop =None):
        self.time, self.distance, self.func = time, distance, func
        if drop is not None:
            self.func_parmas = ft.g_extract_coef_drop(time, distance, func, drop_points=drop)
        else:
            self.func_parmas = ft.g_extract_coef(time, distance, func) 
        
        self.funclam = lambda  x : func( x, *self.func_parmas[0])
        # self.yerr = g_derivate_err( time, distance, self.funclam)  

    def plot( self, max_err= 0.005  ):
        indices = np.abs(self.funclam(self.time) - self.distance) <max_err  
        print(self.func_parmas[-1])
        
        # hook 
        yerr =  np.pi * self.func_parmas[-1] * self.funclam(self.time)[indices]**-1  
        
        # plt.style.use('seaborn-dark-palette')


        plt.errorbar( self.time[indices], self.distance[indices], yerr=yerr , fmt='o', c= "black" ,alpha=0.2) 
        plt.plot( self.time, self.funclam(self.time) , c = next(gColors))
        # plt.show()

def file_type(_str):
    return _str[-3:]

class Exp:

    def __init__(self, _filespath, _filter = None, _metaformat= None ):
        self._files = _filespath
        self.data = []
        self.meta = []

        for _file in self._files:
            tempdata = { 
                "txt" : 
                    lambda : np.array( [  np.fromstring(line, sep=" ")\
                        for line in  open( _file, "r").readlines() ]).transpose(),   
                "csv" :
                    lambda : pd.read_csv(_file)
                }[file_type(_file)]() 

            if  _filter is not None: 
                tempdata = _filter(tempdata)
            self.data.append( tempdata )

            if _metaformat is not None:
                self.meta.append(  _metaformat( _file ) )
    
    def extract_func(self, extract_coef_func):
        self.func_parmas =  [ [popt, pcov, MSR] for popt , (pcov, points), MSR in\
             map( lambda _matrix :  extract_coef_func(_matrix[0], _matrix[1]) , self.data)] 

    
