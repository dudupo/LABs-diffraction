import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(os.path.realpath(".."))
from scipy.signal import find_peaks
from lib.exp import Exp, FittedObj, gencolor, glabels, plotlabels, gColors, putlabel
from lib import fit as ft

import matplotlib
matplotlib.rcParams['text.usetex'] = True
from scipy.ndimage.filters import convolve  
from scipy.ndimage import gaussian_filter


def extract_data():
    _files = [ 
		 	
			"35cm_1hz_1khz.txt",
			"42cm_1hz_1khz.txt",
			#"50cm_1hz_1khz.txt",
			"65cm_1hz_1khz.txt",
			"80cm_1hz_1khz.txt",
			"90cm_1hz_1khz.txt",
                        "110cm_1hz_1khz.txt"
            ]
    
    _dis = [ (10 ** -2) * float(_file.split("cm")[0]) for _file in _files ]
    #_files = sorted(_files)
    L = [ ]
    for _file in _files:
        L.append(
		np.array([ float(line) for line in\
			 open('./data/{0}'.format(_file)).readlines() ] ))
    
    #    plt.plot(L[-1])
    #plt.legend( _dis )
    #plt.show()
    return np.array(L), _dis


def extract_gauss_packet():
    _file = "./data/gauss-2.csv"
    data = []
    for line in open(_file).readlines():  
        [x, y] = line.split(",")
        data.append(float(y))
    return np.array(data)

def fft(funcs):
    ret = [ ] 
    for func in funcs:
        ret.append( np.fft.fft(func))
    return np.array(np.abs(ret)) 

def extract_peaks( func ):
    peaks, _ = np.argmax(func), np.max(func)
    print(peaks, _ )
    return peaks, _ 

def plot_distance_Amp_graph():
    allpeaks, allvals = [], []
    funcs, distance = extract_data()
    for peaks, vals in [extract_peaks(func) for func in fft(funcs)]:
        allpeaks, allvals =  allpeaks+ [peaks], allvals + [vals]
    allvals = np.array(allvals)
    allvals /= allvals[0] 
    
    print( distance , allvals )
    
    def f(x, a):
        return a*(x**-2) 
    popt, (pcov, (_range, values)) , MSR = ft.g_extract_coef( np.array(distance), np.array(allvals), f, p0=[1] )
    print(popt)

    plt.plot( _range, values, c=next(gColors))
    distance = np.array(distance)
    yerr= 0.02* distance**-2
    plt.errorbar(distance, allvals, xerr=[0.02] * len(distance), yerr= yerr , fmt="o",  c='black')
    plt.xlabel(r"distance $ r_[m] $ ")
    plt.ylabel(r"$ Amp_r / Amp_{35 cm} $ ")
    plt.title(r"Amplitude as function of the distance from source") 
    plt.legend([r"$ |\psi(t)|^2 $" ,r"$ Measurred $"])
    import tikzplotlib

    tikzplotlib.save("stage1.tex")
    # plt.savefig("stage1.svg" )
    plt.clf()

'''
    f(x)+g(x+t) = x
    f(x)+g(x+t+t') = y
    | f + g| = f + g + f*g 
    |(f + g)**2| - |(f+g)**2|_0  
'''

from scipy.fft import fftfreq

def main():

    plot_distance_Amp_graph()
    # exit(0)
    data = extract_gauss_packet()
    [part1, part2, part3, part4] = [np.fft.fftshift(np.fft.fft(data))\
            [2500 + 700 * i : 2500 + 700 * (i +1)] for i in range(4) ] 
    for part in  [part1, part2, part3, part4]:
        pass
        #plt.plot(gaussian_filter( part, sigma=1))
        #plt.show()
    plt.plot(data)
    # plt.show()
    plt.clf()
    #plt.plot(np.fft.fft(data))
    #plt.plot(gaussian_filter(np.abs(np.fft.fftshift(np.fft.fft(data))), sigma=3))
    print(len(data))
    X = fftfreq( len(data), 1 / 100)
    #X -= X[0]#((bsi
    plt.plot(X, np.abs(np.fft.fft(data)))
    #    plot_distance_Amp_graph()
    plt.title(r"Restore combined signals")
    plt.legend([r"$ |\psi(t)|^2 $"])
    plt.xlabel(r"$ t_{[ 1/kHz]} $ ")
    plt.ylabel(r"$ Amp  $ ")
    # plt.show()
    # plt.title(r"Restoring the amplitude") 

    import tikzplotlib

    tikzplotlib.save("stage2.tex")
    # plt.savefig("stage1.svg" )
    plt.clf()

if __name__ == "__main__":
    main()
