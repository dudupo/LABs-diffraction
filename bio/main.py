
p=0.5

D = { (0,0) : 1 , (1,0) : 1 , (1,1) : p }



from random import choice
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
if __name__ == "__main__":
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