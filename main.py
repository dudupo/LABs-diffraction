


import numpy as np 
import matplotlib.pyplot as plt



from scipy.optimize import curve_fit

def extract_coef( time, distance ):
    def f(x, a ):
        return np.sin(a*x)/(a*x)   
    popt, pcov = curve_fit(f, time, distance)
    _range = np.linspace( np.min(time) , np.max(time), 100 )
    return popt, (_range, f(_range, *popt)) 

_files = [
        "single0.02Bblack.txt",  "single0.04.txt",
        "single0.02no.txt",      "single0.08.txt",
        "single0.02B.txt",  "single0.02no2.txt",     "single0.08a.txt" ]




if __name__ == "__main__":

    for _file in _files:

        A = np.array( [  np.fromstring(line, sep=" ")  for line in  open( _file, "r").readlines() ])
        A = A.transpose()
        
        print(A.shape)
        
        ind = np.abs(A[-1]) < 1 

        A[-2] = A[-2] / np.sqrt( np.sum( A[-2]**2 ) )
        try:

            _, ( X, Y )= extract_coef(A[-1][ind], -A[-2][ind])
            plt.plot(X,Y)
        except:
            print("error {0}".format(_file))
            pass
        
        plt.plot(A[-1][ind] , -A[-2][ind] )

    plt.legend( [ _file[_file.index("."):_file.rindex(".")] for _file in _files ]   )
    plt.show()