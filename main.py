


import numpy as np 
import matplotlib.pyplot as plt



from scipy.optimize import curve_fit

def extract_coef( time, distance ):

    def f(x, a, b, c, d):
        return ( b* np.sin(a*x + c)/(a*x + c) )**2  
    popt, pcov = curve_fit(f, time, distance)
    _range = np.linspace( np.min(time) , np.max(time), 1000 )

    values = f(_range, *popt)
    MSR = (values - distance)**2/ len(values)
    return popt, (pcov, (_range, values)) , MSR 


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
    "single0.04Good2.txt",
    "single0.04Good.txt",
    "single0.08Good2.txt",
    "single0.08Good.txt",
    "single0.16Good2.txt",
    "single0.16Good.txt",
]

def genLegend( _file ):
    return _file[_file.index("."):_file.rindex(".")]

if __name__ == "__main__":


    legends = []
    for _file in _files:

        A = np.array( [  np.fromstring(line, sep=" ")  for line in  open( _file, "r").readlines() ])
        A = A.transpose()
        
        print(A.shape)
        
        ind = np.abs(A[-1]) < 1 
        A[-2] = A[-2] / np.sqrt( np.sum( A[-2]**2 ) )
        
        
        plt.plot(A[-1][ind] , -A[-2][ind] )
        # plt.plot( A[-1], -A[-2])
        legends.append(genLegend(_file) )
        
        try:
            pass

            __ , (_, ( X, Y )), MSR = extract_coef(A[-1][ind], -A[-2][ind])
            plt.plot(X,Y)
            print(_)
            legends.append( genLegend( _file ) + " fit"  )
        except:
            print("error {0}".format(_file))
            pass
        plt.legend( legends )
        plt.savefig( _file +".svg" )
        legends = []
        plt.clf()