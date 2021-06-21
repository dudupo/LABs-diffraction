from random import random
import numpy as np
import sys, os
sys.path.append(os.path.realpath(".."))

from bio.mc import ColonyRect
from bio.utility import histogram_calc

def exp(k, p, initx=1):
	#initx = initx * int((random()+1) * 1.8)
	x =  initx
	L = [x]
	while x < initx*k:
		for _ in range(x):
			x = x+1 if random() < p else x 
		L.append( x / initx )        
	return np.array(L)

def sim(max_mul_size, max_time, number_of_exp=300, p =0.125 ):
	L = [ ]
	for _ in range(number_of_exp):
		L.append( exp(max_mul_size, p, 1) )
	return histogram_calc( L, 1, max_time, max_mul_size )

def sim_inital_time(colony_vec, p , max_mul_size , max_time ):
	L = [ exp(max_mul_size, p, colony_vec[0].pix_num ) for colony in colony_vec ] 
	return histogram_calc( L, 1, max_time, max_mul_size )
		



