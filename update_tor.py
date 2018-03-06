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
	print(op.minimize(run_GRASP, x0=[1.2, .7, .4, -1.0], args=(results_path), method='Nelder-Mead', bounds=((0, 1.5), (0, 1.5), (0, .75), (-1.2, 0)) ))
	


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
	
def run_GRASP(parameters, results_path):
	# results_path = results_path[0]
	print(results_path)
	plt.rc('axes', linewidth=1)
	GRASP_working_file = "F:/Devin/Grasp/LWASandbox/40mQuadDipole/working/"
	model_file = "40mQuadDipole.tor"
	out_file = GRASP_working_file + "40mQuadDipole.tor"


	separation 	= parameters[0]
	x 	= parameters[1]
	y 	= parameters[2]
	z	= parameters[3]

	edit_msh(x, y, z, GRASP_working_file)

	#########  Generate folders  ###########
	results_path += "/separation=%4.2f_x=%4.2f_y=%4.2f_z=%4.2f" % tuple(parameters)
	# results_path += "/separation=%4.2f" % tuple(parameters)
	if not os.path.exists(results_path):
		os.mkdir(results_path)
	print(results_path)

	for direct in ["/plots/", "/plots/Efficiency/", "/plots/SEFD/", "/plots/Patterns/", "/data/"]:
		if not os.path.exists(results_path + direct):
			os.mkdir(results_path + direct)

	losses = []
	for AntPos in np.linspace(13,17.5, 15):

		change_list = [489, 333]
		mod_line 	= ["  value            : %4.2f\n" % AntPos, "  value            : %4.2f\n" % separation]
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