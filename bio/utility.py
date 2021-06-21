import numpy as np
import sys, os
sys.path.append(os.path.realpath(".."))


def convert_colonys_to_ktlist(colonys):
	return [ [_frame.pix_num / colony[0].pix_num \
		for _frame in colony] for colony in colonys ]

def histogram_calc( _list, eps, max_time, max_mul_size):
	empty_propb = np.zeros( (max_mul_size * eps, max_time) )	

	def calc_pos(k):
		return int( k * eps ) 

	for vec in _list:
		for t,k in enumerate(vec):
			if calc_pos(k) < empty_propb.shape[0] and t < max_time:
				empty_propb[ calc_pos(k) ][t] += 1
			# else:
			# 	if calc_pos(k) > empty_propb.shape[0]:
			# 		empty_propb = np.r_[empty_propb,
			# 		 np.zeros((1 + calc_pos(k) - empty_propb.shape[0] , max_time))]
			# 		empty_propb[calc_pos(k)][t] += 1
	print(empty_propb.shape)
	return empty_propb
	
	