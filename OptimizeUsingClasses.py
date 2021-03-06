import numpy as np 
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt 
import platform
import os
import shutil

import scipy.optimize as op
try:
	from simanneal import Annealer
except:
	pass

from AntennaClasses import *


def main():

	## HELIOS
	if platform.node() == "Helios":
		a = ELfeedExt(start_f =  40.0, end_f = 100.0, n_f = 50, alpha = 0, grasp_version = 10.3)
		a.set_number_of_focal_lengths(1)

		print("Executing on Helios")
		print("%s"%a)
		a.set_global_directory_name("/mnt/f/Documents/Caltech/LWA/GRASP/")
		a.set_ticra_directory_name("/mnt/f/Program Files/TICRA/")
		a.set_grasp_analysis_extension(".exe")

	## G1 (ASTROS)
	elif platform.node() == 'DESKTOP-3UVMJQF' or platform.node() == 'ASTROS':
		# a = ELfeedExt(start_f =  60.0, end_f = 85.0, n_f = 10, alpha = 0, grasp_version = 10.3)
		# a.set_number_of_focal_lengths(1)
		a = QRFH(freq = 60, grasp_version = 10.6, z_phase = 340, off_axis = 0.67)
		a.set_number_of_focal_lengths(1)

		print("Executing on G1 Office")
		print("%s"%a)
		a.set_global_directory_name("/mnt/c/Users/dcody/Documents/GRASP/")
		a.set_ticra_directory_name("/mnt/c/Program Files/TICRA/")
		a.set_grasp_analysis_extension(".exe")

		cst_dir = "/mnt/c/Users/dcody/Documents/GRASP/40mQRFHsim106/DSAfeedPats/"
		# cst_dir = "/mnt/c/Users/dcody/Documents/GRASP/40mQRFHsim106/QRFHfarfields/"

	## AWS
	elif platform.node() == 'ip-172-31-33-156':
		a = Elfeed(start_f =  60.0, end_f = 85.0, n_f = 10, alpha = 0, grasp_version = 10.3)
		a.set_number_of_focal_lengths(1)

		print("Executing on AWS")
		print("%s"%a)
		a.set_global_directory_name("/home/ubuntu/GRASP/")
		a.set_ticra_directory_name("/home/ubuntu/TICRA/")
		a.set_grasp_analysis_extension()

	## MOORE
	else:

		#a = ELfeed(start_f =  60.0, end_f = 85.0, n_f = 10, alpha = 0, grasp_version = 10.3)
		a = QRFH(freq = 60, grasp_version = 10.6, z_phase = -160)#340)
		a.set_number_of_focal_lengths(5)

		print("Executing on Moore")
		print("%s"%a)
		a.set_global_directory_name()
		a.set_ticra_directory_name()#"/cygdrive/c/Program Files/TICRA/")
		a.set_grasp_analysis_extension()

		cst_dir = "F:\\Devin\\CST\\QRFH\\qrfh_v0_aper_circ_HF_donutnewnew_DC_COPY_noscale\\Result\\"
		cst_dir = "F:\\Devin\\CST\\Midband Antenna\\spdipole_midband\\Result\\"
		#cst_dir = "F:\\Devin\\CST\\DSAfeed\\Antenna DSA\\Result\\"



	if 0:
		a.parameters["x"] = 		0.739
		a.parameters["y"] = 		0.115
		a.parameters["z"] = 		-0.051
		a.parameters["alpha"] = 	0.000
		a.parameters["sp"] =		0.946
		a.parameters["el"] =		1.966
		a.parameters["ew"] =		0.590
		a.parameters["ed"] =		0.000
		
		a.parameters["z_dist"] = 16.50
		a.bounds.update({"z_dist":[16.50, 17.5]})
		random(a)
		# setup_configuration(a)
		# simulate_single(a, override_frequency = False, plot_feed = True)

	if 1:
		a.parameters["z_dist"] = 15.86
		a.bounds.update({"z_dist":[15.86, 17.5]})
		iterate_over_cut_files(a, cst_dir, frequency_scale = 1, freq_scale_to_MHz=1)#5.75)
		#if frquencies in file are in MHz, then freq_scale_to_MHz = 1
		#if frquencies in file are in GHz, then freq_scale_to_MHz = 1000


	if 0:
		a.parameters["x"] = 	0.965
		a.parameters["y"] = 	0.011
		a.parameters["z"] = 	0.218
		# a.parameters["sp"] =	1.068
		# a.parameters["el"] =	1.967
		# a.parameters["ew"] =	0.708
		# a.parameters["ed"] =	0.000
		a.bounds.update({"z_dist":[16.5,17.5]})
		random(a)

	# nelder_mead(a, x)
	if 0:
		x=[0.970, 0.801, 0.347]
		a.bounds.update({"z_dist":[16.5,17.5]})
		nelder_mead(a, x)


	# nelder_mead2(a)
	if 0:
		a.bounds.update({"z_dist":[16.5,17.5]})
		random(a)
	# setup_configuration(a)


	# grid(a)

	# simulate_single(a, override_frequency = False, plot_feed = True)

	# simulate_single(a, plot_feed = True, override_frequency = False)
	# a.parameters["x"] = 		0.990
	# a.parameters["y"] = 		0.737
	# a.parameters["z"] = 		0.501
	# simulate_single(a, plot_feed = True, override_frequency = False)



	# iterate_over_cut_files(a, cst_dir)





def anneal(a):
	class AntennaProblem(Annealer):
		def __init__(self, state, antenna):
			self.a = antenna
			self.names = self.a.get_optimizable_parameter_names()
			self.nxt_energy = self.a.simulate_single_configuration([], [], plot_feed = True)

			super(AntennaProblem, self).__init__(state)

		def move(self):
			e = 1										# Define energy which is returned
			while e > 0:								# IF invalid configuration, then
				for idx, k in enumerate(self.names):
					if (k != "alpha"):
						bnd = self.a.bounds[k][1] - self.a.bounds[k][0]
						self.a.parameters[k] = self.state[k] + np.random.uniform(-bnd/10, bnd/10)	# Update antenna class
					else:
						if (np.random.rand() > .5):
							self.a.parameters["alpha"] = 45
						else:
							self.a.parameters["alpha"] = 0
				e = self.a.simulate_single_configuration([], [], plot_feed = True)	# Calc returnable energy

			for idx, k in enumerate(self.names):				# Should have a vaild configuration now...
				self.state.update({k:self.a.parameters[k]})		# update the state with all the antenna class parameters
			self.nxt_energy = e

		def energy(self):
			return self.nxt_energy

	setup_simulation_files(a, "anneal")
	ap = AntennaProblem(a.parameters, a)
	ap.steps = 300
	ap.tmax = 2000
	ap.anneal()

def iterate_over_cut_files(a, cst_dir, frequency_scale = 5.75, freq_scale_to_MHz=1):
		#if frquencies in file are in MHz, then freq_scale_to_MHz = 1
		#if frquencies in file are in GHz, then freq_scale_to_MHz = 1000

	setup_simulation_files(a, "po_tabs")
	a.include_freq_in_title = False
	#move new file to working directory overwriting last
	#a.GRASP_working_file
	files = os.listdir(cst_dir)
	c = 0
	for file in files:
		if ("[1]_theta-phi).cut" in file[-20:]):
			if (c%3 == 0):
				print (file)
				freq = float(file.split()[1][3:-1]) #find the freq
				print  (freq*freq_scale_to_MHz)
				a.parameters["freq"] = freq_scale_to_MHz*freq/frequency_scale
				shutil.copy2(cst_dir + file, a.GRASP_working_file + "pat.cut")

				a.simulate_single_configuration([],[], plot_feed = False, override_frequency = True, off_axis = True)
			c += 1


def setup_simulation_files(a, method_name):
	a.set_method_name(method_name)
	a.init_global_file_log()
	a.gen_file_names()


def simulate_single(a, plot_feed = True, override_frequency = False):
	setup_simulation_files(a, "sing")
	a.simulate_single_configuration([],[], plot_feed = plot_feed, override_frequency = override_frequency)	

def grid(a):
	setup_simulation_files(a, "grid")
	bounds = a.get_bounds()


	names = a.get_optimizable_parameter_names()
	try:
			names.remove("alpha")
	except:
		pass
	# print(names, bounds)
	number_of_points = 12

	for i in range(number_of_points):
		nm = "taper"
		x_new = []
		x_new.append(bounds[nm][0] + (bounds[nm][1] - bounds[nm][0])/(number_of_points-1)*i)
		for j in range(number_of_points):
			nm = "angle"
			x_new.append(bounds[nm][0] + (bounds[nm][1] - bounds[nm][0])/(number_of_points-1)*j)
			print ("loss = ", a.simulate_single_configuration(x_new, names, plot_feed = False, override_frequency = True))
			x_new = x_new[:-1]

def setup_configuration(a):
	setup_simulation_files(a, "setup")

	if a.model_name == "40mQuadDipole_High_Freq":
		scale = 0.1
	else:
		scale = 1

	a.edit_msh()
	a.edit_tor()

def nelder_mead(a, x):
	method = 'Nelder-Mead'

	setup_simulation_files(a, "NM")
	
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
		if (np.random.uniform(0,1) < 0.7):
			a.parameters["alpha"] = 0
		else:
			a.parameters["alpha"] = 45
			
		print ("loss = ", a.simulate_single_configuration(x_new, names, plot_feed = True))


if __name__ == '__main__':
	plt.rc('axes', linewidth=1)
	main()