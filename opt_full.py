import numpy as np 
import matplotlib.pyplot as plt 
import process_grasp
import subprocess
import time
from datetime import datetime
import pandas as pd
import os
import scipy.optimize as op

def main():
	results_path = gen_results_folder()
	path = tuple([results_path])
	# print (run_GRASP([1.2, .7, .4, -1.0], path))

	x_init = [1.06, .8, .19, -1.34, 1.23, .482, .232]
			 #sep,    x,   y,     z, dirL, dirW, dirS
	method = 'Nelder-Mead'
	#bnd = ((0, 1.5), (0, 1.5), (0, .75), (-1.5, 0), (0,3.5), (0, 2), (0,1.5))
		   #sep,      x,        y,        z,         dirL,    dirW,   dirS
	#constr = ({'type':'ineq', 'fun': lambda x: x[0] - x[1] - x[2] - 0.08 - 0.012}) # 
	print(op.minimize(run_GRASP, x0=x_init, args=(results_path), method=method))#, bounds= bnd, constraints = constr))
	


def edit_tor(template, out_file, change_list, modified_line): #template 
	print("Opening Files")
	f = open(template,'r')
	g = open(out_file,'w+')

	for i, line in enumerate(f): #0-indexed line number
		if i in change_list:
			g.write(modified_line[change_list.index(i)])
		else:
			g.write(line)

	g.close()
	f.close()
	print("Done writing TOR")

def edit_msh(x, y, z, GRASP_working_file):
	msh_template = "patch_leaf_inverted_V.msh"
	msh_out = GRASP_working_file + msh_template

	modified_line = []
	change_list = range(28,39,2)
	origin = np.array([[.08,-.012, 1.5],[.08,    0, 1.5],[.08, .012, 1.5]])
	point_data = [[2,  1],	#10
				  [0,  1],  #11
				  [2, .5], 	#12
				  [1,  1],  #13
				  [0, .5],	#14
				  [1, .5]]	#15

	for pt in point_data:
		diff = np.array([x, y, z])
		if pt[0] == 0:
			diff[1] *= -1
		elif pt[0] == 1:
			diff[1] = 0
		diff *= pt[1]
		modified_line.append("%6.3f  %6.3f  %6.3f\n" % tuple(diff+origin[pt[0]]))

	f = open(msh_template,'r')
	g = open(msh_out,'w+')

	for i, line in enumerate(f): #0-indexed line number
		if i in change_list:
			g.write(modified_line[change_list.index(i)])
		else:
			g.write(line)

	g.close()
	f.close()
	print("Done writing MSH")


def create_file_w_timestamp(location = os.getcwd()):
	file_name = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
	mydir = os.path.join(
		location, 
		file_name)
	try:
		os.makedirs(mydir)
	except OSError as e:
		if e.errno != errno.EEXIST:
			raise Exception('File already exists!')
	return file_name


def gen_results_folder():
	results_path = "../Results/"
	results_path += create_file_w_timestamp(results_path)
	return results_path
	

def check_ant_intersection(parameters):
	return parameters[0] - parameters[1] - parameters[2] - 0.08 - 0.012

def run_GRASP(parameters, results_path):
	# results_path = results_path[0]
	print(results_path)
	plt.rc('axes', linewidth=1)
	GRASP_working_file = "F:/Devin/Grasp/LWASandbox/40mQuadDipoleWDir/working/"
	model_file = "40mQuadDipoleWDir.tor"
	out_file = GRASP_working_file + "40mQuadDipoleWDir.tor"


	separation 	= parameters[0]
	x 	= parameters[1]
	y 	= parameters[2]
	z	= parameters[3]
	dir_len = parameters[4]
	dir_wid = parameters[5]
	dir_sep = parameters[6]

	

	#########  Generate folders  ###########
	results_path += "/sp=%4.2f_x=%4.2f_y=%4.2f_z=%4.2f_dl=%4.2f_dw=%4.2f_ds=%4.2f" % tuple(parameters)
	if not os.path.exists(results_path):
		os.mkdir(results_path)
	print(results_path)

	# Check out of Bounds
	bnd = [[.75, 2.0], [0, 1.5], [0, .75], [-2.5, 0], [0,3.5], [0, 2], [0,1.5]]
	error = 0.0
	for i, par in enumerate(parameters):
		if par < bnd[i][0]:
			print("Out of bounds, parameter: %d", i)
			error += (par-bnd[i][0])**2
		if par > bnd[i][1]:
			print("Out of bounds, parameter: %d", i)
			error += (par-bnd[i][1])**2
	if 0 < check_ant_intersection(parameters):
		print("Antenna Intersection")
		error += check_ant_intersection(parameters)**2
	if error > 0.000000001:
		try:
			os.rename(results_path, results_path + "_minLoss = %4.3f" % error)
		except:
			print("Error: Cannot rename file %f" %results_path)
		return error

	# Generate Sub Folders
	for direct in ["/plots/", "/plots/Efficiency/", "/plots/SEFD/", "/plots/Patterns/", "/data/"]:
		if not os.path.exists(results_path + direct):
			os.mkdir(results_path + direct)

	#Edit the Antenna Mesh File
	edit_msh(x, y, z, GRASP_working_file)

	# Optimize over phase center
	losses = []
	for AntPos in np.linspace(14, 16.7, 10):

		change_list = [489, 333, 502, 507, 512]
		       #dipole_sep, sep, dirL,dirW,dirS
		mod_line 	= ["  value            : %6.3f\n" % AntPos, "  value            : %6.3f\n" % separation,
					   "  value            : %6.3f\n" % dir_len, "  value            : %6.3f\n" % dir_wid,
					   "  value            : %6.3f\n" % dir_sep]
		edit_tor(model_file, out_file, change_list, mod_line) # <====== edit all tor files here


		command = "grasp-analysis batch.gxp out.out out.log"
		process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, cwd = GRASP_working_file)
		output, error = process.communicate()
		# print("Done executing grasp")
		#print(output, error)

		#########  Process data files  ###########
		# print("preparing to Plot")
		freq, s11 = process_grasp.process_par(GRASP_working_file + "S_parameters.par")
		dmax, cut = process_grasp.process_cut(GRASP_working_file + "Field_Data.cut", freq)
		mis = process_grasp.calc_mismatch(s11)
		aperture = process_grasp.calc_app_eff(freq, dmax)
		tsys = process_grasp.Tsys(freq)
		SEFD = process_grasp.SEFD(freq, dmax)
		loss1 = loss(freq, aperture)
		losses.append(loss1)

		plots_directory = results_path + "/plots/"
		data_directory = results_path + "/data/"
		pattern_directory = plots_directory + "Patterns/AntPos=%4.2f_loss=%4.2f/"% (AntPos, loss1)

		if not os.path.exists(pattern_directory):
			os.mkdir(pattern_directory)

		process_grasp.plot_pair_efficiencies(freq, s11, dmax, plots_directory + "Efficiency/Efficiencies Antenna Position = %4.2f Loss = %4.3f.png" % (AntPos, loss1), AntPos)
		process_grasp.plot_SEFD(freq, dmax, plots_directory + "SEFD/SEFD Antenna Position = %4.2f Loss = %4.3f.png" % (AntPos, loss1), AntPos)
		for frequency in freq:
			process_grasp.plot_cut(frequency, cut, AntPos, pattern_directory +"Radiation Pattern Position = %4.2f Freq = %4.2f Loss = %4.3f.png" %(AntPos,frequency, loss1))


		x = open(data_directory+"Results AntPos=%4.2f Loss=%4.3f.csv" %(AntPos,loss1), 'w+')
		x.write(time.strftime("%c"))
		x.write("\n")
		x.write("Frequency [MHz], S11 [dB], Dmax [dBi], Mismatch Efficiency, Aperture Efficiency, Tsys [1E3 K], SEFD [1E6 Jy]\n")
		for ii in range(len(freq)):
			x.write("%6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f\n" % (freq[ii], s11[ii], dmax[ii], mis[ii], aperture[ii], tsys[ii]/1E3, SEFD[ii]/1E6))
		x.close()

	minLoss = np.min(losses)
	try:
		os.rename(results_path, results_path + "_minLoss = %4.3f" % minLoss)
	except:
		print("Error: Cannot rename file %f" %results_path)
	return minLoss

def loss(freq, effic):
	#calculate the loss funciton
	ans = 0
	n = 0
	for i, f in enumerate(freq):
		if (f >= 60 and f <= 80):
			ans -= effic[i]
			n+= 1
	return ans/n




if __name__ == '__main__':
	main()