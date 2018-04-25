import numpy as np 
import matplotlib.pyplot as plt 

import scipy.optimize as op

from AntennaClasses import ELfeedDir, LWA_like


def main():
	a = ELfeedDir(start_f =  60.0, end_f = 80.0, n_f = 5, alpha = 0, grasp_version = 10.6)
	a.set_number_of_focal_lengths(5)

	names = a.get_parameter_names()
	
	# remove parameters which are altered multiple times (e.g. z_dist)
	# or parameters that are altered once per execution (e.g. n_f)
	for x in ["z_dist", "start_f", "end_f", "n_f", "alpha"]:
		names.remove(x)


	bounds = a.get_bounds()
	# print (bounds["x"])


	####Update parameters
	# 
	# x_new = list(x)
	# for i in range(2000):
	# 	x_new = []
	# 	for idx, k in enumerate(names):
	# 		# print(k)
	# 		x_new.append(np.random.uniform(bounds[k][0], bounds[k][1]))
	# 	print ("loss = ", a.simulate_single_configuration(x_new, names))


	method = 'Nelder-Mead'
	# method = 'Powell'
	# x = [1.05, .77, .16, -.01, 1.06, .6, .3]
	x = [1.307, .087, -.15]
	# 	   #sep,      x,        y,        z,         dirL,    dirW,   dirS
	print(op.minimize(a.simulate_single_configuration, x0=x, args=(names), method=method))#, bounds= bnd, constraints = constr))


if __name__ == '__main__':
	plt.rc('axes', linewidth=1)
	main()