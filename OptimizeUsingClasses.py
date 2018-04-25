import numpy as np 
import matplotlib.pyplot as plt 
import platform

import scipy.optimize as op

from AntennaClasses import ELfeedDir, LWA_like


def main():
	a = LWA_like(start_f =  60.0, end_f = 80.0, n_f = 5, alpha = 0, grasp_version = 10.6)
	a.set_number_of_focal_lengths(5)
	a.init_global_file_log()


	if platform.node() == "Helios":
		print("Executing on Helios")
		a.set_global_directory_name("/mnt/f/Documents/Caltech/LWA/GRASP/")
		a.set_ticra_directory_name("/mnt/f/Program Files/TICRA/")
		a.set_grasp_analysis_extension(".exe")
	elif platform.node() == 'DESKTOP-3UVMJQF':
		print("Executing on G1 Office")
		a.set_global_directory_name("/mnt/c/Users/dcody/Documents/GRASP/")
		a.set_ticra_directory_name("/mnt/c/Program Files/TICRA/")
		a.set_grasp_analysis_extension(".exe")
	else:
		print("Executing on Moore")
		a.set_global_directory_name()
		a.set_ticra_directory_name()
		a.set_grasp_analysis_extension()
	
	a.gen_file_names()
	

	# remove parameters which are altered multiple times (e.g. z_dist)
	# or parameters that are altered once per execution (e.g. n_f)
	names = a.get_optimizable_parameter_names()
	names.remove("alpha")
	
	nelder_mead(a, names)


	# bounds = a.get_bounds()
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

def nelder_mead(a, names):
	method = 'Nelder-Mead'
	# method = 'Powell'
	# x = [1.05, .77, .16, -.01, 1.06, .6, .3]

	x = [.807, .087, -.01]
	# 	   #sep,      x,        y,        z,         dirL,    dirW,   dirS
	print(op.minimize(a.simulate_single_configuration, x0=x, args=(names), method=method))#, bounds= bnd, constraints = constr))


if __name__ == '__main__':
	plt.rc('axes', linewidth=1)
	main()