import numpy as np 
import matplotlib.pyplot as plt 
import platform

import scipy.optimize as op

from AntennaClasses import *


def main():

	## HELIOS
	if platform.node() == "Helios":
		a = ELfeedExt(start_f =  60.0, end_f = 80.0, n_f = 5, alpha = 0, grasp_version = 10.3)
		a.set_number_of_focal_lengths(5)

		print("Executing on Helios")
		print("%s"%a)
		a.set_global_directory_name("/mnt/f/Documents/Caltech/LWA/GRASP/")
		a.set_ticra_directory_name("/mnt/f/Program Files/TICRA/")
		a.set_grasp_analysis_extension(".exe")

	## G1
	elif platform.node() == 'DESKTOP-3UVMJQF':
		a = ELfeedExt(start_f =  60.0, end_f = 80.0, n_f = 5, alpha = 0, grasp_version = 10.3)
		a.set_number_of_focal_lengths(5)

		print("Executing on G1 Office")
		print("%s"%a)
		a.set_global_directory_name("/mnt/c/Users/dcody/Documents/GRASP/")
		a.set_ticra_directory_name("/mnt/c/Program Files/TICRA/")
		a.set_grasp_analysis_extension(".exe")

	## MOORE
	else:
		a = ELfeedExt(start_f =  60.0, end_f = 80.0, n_f = 5, alpha = 0, grasp_version = 10.3)
		a.set_number_of_focal_lengths(5)

		print("Executing on Moore")
		print("%s"%a)
		a.set_global_directory_name()
		a.set_ticra_directory_name()#"/cygdrive/c/Program Files/TICRA/")
		a.set_grasp_analysis_extension()
	
	a.set_method_name("general")
	
	

	# remove parameters which are altered multiple times (e.g. z_dist)
	# or parameters that are altered once per execution (e.g. n_f)
	
	x = [0.8756, 0.083, -0.6018, 1.0773, 1.1448, 0.6672] 

	# random(a)
	nelder_mead(a, x)
	# nelder_mead2(a)
	# random(a)
	# nelder_mead(a)
	# setup_configuration(a)
	# simulate_single(a)

def setup_simulation_files(a, method_name):
	a.set_method_name(method_name)
	a.init_global_file_log()
	a.gen_file_names()


def simulate_single(a):
	setup_simulation_files(a, "sing")
	a.parameters["x"] = 		0.7835
	a.parameters["y"] = 		0.048
	a.parameters["z"] = 		0.5627
	a.parameters["sp"] =		1.0254
	a.parameters["dl"] =		1.7321
	a.parameters["dw"] =		0.3172
	a.parameters["dsep"] =		-1.2299

	# a.parameters["rl"] = 1.06
	# a.parameters["rw"] = .60
	# a.parameters["rsep"] = -.3
	
	a.parameters["z_dist"] = 16.4

	x=[]

	names = a.get_optimizable_parameter_names()
	names.remove("alpha")
	for nam in names:
		x.append(a.parameters[nam])
	a.simulate_single_configuration(x, names, plot_feed = True)	

def setup_configuration(a):
	setup_simulation_files(a, "setup")

	if a.model_name == "40mQuadDipole_High_Freq":
		scale = 0.1
	else:
		scale = 1


	a.parameters["x"] = 		0.6332
	a.parameters["y"] = 		0.1347
	a.parameters["z"] = 		0.6752
	a.parameters["sp"] =		0.9175
	# a.parameters["dl"] =		1.7321
	# a.parameters["dw"] =		0.3172
	# a.parameters["dsep"] =		-1.2299

	# a.parameters["rl"] = 1.06
	# a.parameters["rw"] = .60
	# a.parameters["rsep"] = -.3


	a.parameters["z_dist"] = 16.625

	print(a.parameters["x"])

	a.edit_msh()
	a.edit_tor()

def nelder_mead(a, x):
	method = 'NM'

	setup_simulation_files(a, method)
	
	names = a.get_optimizable_parameter_names()
	names.remove("alpha")
	# method = 'Powell'

	# # x = [1.0761, .7498, 0.4728]  						#LWA_LIKE
	# # x = [1.0761, .7498, -3, 1.43, .8]  						#LWA_DIR
	# x =  				#11
	# # x = [1.05, .25, -.16, 1.5, 1.06, .6, .3]			#11 DIR
	# # x = [0.6604, 0.356, -0.3087, 1.1227, 0.1509, 0.0748, -2.02]	#11 REF

	# # x = [0.8161,0.0318,-0.3301,0.9399,0.7703,0.0045, 0.1227]  #11 DIR

	for name, val in zip(names, x):
		print (name,":=  ", val)

	
	# 	   #sep,      x,        y,        z,         dirL,    dirW,   dirS
	print(op.minimize(a.simulate_single_configuration, x0=x, args=(names, True), method=method))#, bounds= bnd, constraints = constr))


def random(a):
	setup_simulation_files(a, "Rand")
	bounds = a.get_bounds()


	names = a.get_optimizable_parameter_names()
	names.remove("alpha")
	# print names

	for i in range(20000):
		x_new = []
		for idx, k in enumerate(names):
			# print(k)
			x_new.append(np.random.uniform(bounds[k][0], bounds[k][1]))
		print ("loss = ", a.simulate_single_configuration(x_new, names, plot_feed = True))


if __name__ == '__main__':
	plt.rc('axes', linewidth=1)
	main()