import numpy as np 
import matplotlib.pyplot as plt 
import platform

import scipy.optimize as op

from AntennaClasses import *


def main():
	## HELIOS
	if platform.node() == "Helios":
		a = ELfeedDirRef(start_f =  60.0, end_f = 80.0, n_f = 5, alpha = 0, grasp_version = 10.3)
		a.set_number_of_focal_lengths(5)
		a.init_global_file_log()

		print("Executing on Helios")
		a.set_global_directory_name("/mnt/f/Documents/Caltech/LWA/GRASP/")
		a.set_ticra_directory_name("/mnt/f/Program Files/TICRA/")
		a.set_grasp_analysis_extension(".exe")

	## G1
	elif platform.node() == 'DESKTOP-3UVMJQF':
		a = ELfeedRef(start_f =  60.0, end_f = 80.0, n_f = 5, alpha = 0, grasp_version = 10.3)
		a.set_number_of_focal_lengths(5)
		a.init_global_file_log()


		print("Executing on G1 Office")
		a.set_global_directory_name("/mnt/c/Users/dcody/Documents/GRASP/")
		a.set_ticra_directory_name("/mnt/c/Program Files/TICRA/")
		a.set_grasp_analysis_extension(".exe")

	## MOORE
	else:
		a = ELfeedDir(start_f =  60.0, end_f = 80.0, n_f = 5, alpha = 0, grasp_version = 10.3)
		a.set_number_of_focal_lengths(5)
		a.init_global_file_log()


		print("Executing on Moore")
		a.set_global_directory_name()
		a.set_ticra_directory_name()#"/cygdrive/c/Program Files/TICRA/")
		a.set_grasp_analysis_extension()
	
	a.gen_file_names()
	

	# remove parameters which are altered multiple times (e.g. z_dist)
	# or parameters that are altered once per execution (e.g. n_f)
	


	# random(a)
	# nelder_mead(a)
	# nelder_mead2(a)
	random(a)
	# nelder_mead(a).



def simulate_single(a):

	a.parameters["sp"] = 1.05
	a.parameters["x"] = .77
	a.parameters["y"] = .16
	a.parameters["z"] = -0.01
	a.parameters["dl"] = 1.06
	a.parameters["dw"] = .60
	a.parameters["dsep"] = .3

	a.parameters["rl"] = 1.06
	a.parameters["rw"] = .60
	a.parameters["rsep"] = -.3
	
	a.parameters["z_dist"] = 16.4

	x=[]

	names = a.get_optimizable_parameter_names()
	names.remove("alpha")
	for nam in names:
		x.append(a.parameters[nam])
	a.simulate_single_configuration(x, names, plot_feed = True)	

def setup(a):

	a.parameters["sp"] = 1.03
	a.parameters["x"] = .77
	a.parameters["y"] = .16
	a.parameters["z"] = -.5
	a.parameters["dl"] = 1.46
	a.parameters["dw"] = .63
	a.parameters["dsep"] = 1.37

	a.parameters["rl"] = 1.06
	a.parameters["rw"] = .60
	a.parameters["rsep"] = -.3
	a.parameters["z_dist"] = 16.625

	print(a.parameters["x"])

	a.edit_msh()
	a.edit_tor()

def nelder_mead(a):
	method = 'Nelder-Mead'
	names = a.get_optimizable_parameter_names()
	names.remove("alpha")
	# method = 'Powell'

	# x = [1.0761, .7498, 0.4728]  						#LWA_LIKE
	x = [1.0761, .7498, -3, 1.43, .8]  						#LWA_DIR
	# x = [1.05, .25, -.16, 1.5, 1.06, .6, .3]			#11 DIR
	for name, val in zip(names, x):
		print (name,":=  ", val)

	
	# 	   #sep,      x,        y,        z,         dirL,    dirW,   dirS
	print(op.minimize(a.simulate_single_configuration, x0=x, args=(names, True), method=method))#, bounds= bnd, constraints = constr))

def nelder_mead2(a):
	method = 'Nelder-Mead'
	names = a.get_optimizable_parameter_names()
	names.remove("alpha")
	# method = 'Powell'

	# x = [1.0761, .7498, 0.4728]  						#LWA_LIKE
	x = [1.0761, .7498, 0, 1.43, .8]  						#LWA_DIR
	# x = [1.05, .25, -.16, 1.5, 1.06, .6, .3]			#11 DIR
	for name, val in zip(names, x):
		print (name,":=  ", val)
	
	# 	   #sep,      x,        y,        z,         dirL,    dirW,   dirS
	print(op.minimize(a.simulate_single_configuration, x0=x, args=(names, True), method=method))#, bounds= bnd, constraints = constr))


def random(a):
	bounds = a.get_bounds()


	names = a.get_optimizable_parameter_names()
	names.remove("alpha")

	for i in range(2000):
		x_new = []
		for idx, k in enumerate(names):
			# print(k)
			x_new.append(np.random.uniform(bounds[k][0], bounds[k][1]))
		print ("loss = ", a.simulate_single_configuration(x_new, names, plot_feed = True))


if __name__ == '__main__':
	plt.rc('axes', linewidth=1)
	main()