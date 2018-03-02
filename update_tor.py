import numpy as np 
import matplotlib.pyplot as plt 
import process_grasp
import subprocess
import time
from datetime import datetime
import pandas as pd
import os

def main():
	results_path = gen_results_folder()
	print (run_GRASP([1,2,3], results_path))
	print(run_GRASP([4.1,2,3], results_path))

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

	plt.rc('axes', linewidth=1)
	GRASP_working_file = "F:/Devin/Grasp/LWASandbox/40mQuadDipole/working/"
	model_file = "40mQuadDipole.tor"
	out_file = GRASP_working_file + "40mQuadDipole.tor"


	x 	= parameters[0]
	y 	= parameters[1]
	z 	= parameters[2]
	#	= parameters[3]

	#########  Generate folders  ###########

	results_path += "/x=%4.2f_y=%4.2f_z=%4.2f" % tuple(parameters)
	if not os.path.exists(results_path):
		os.mkdir(results_path)
	print(results_path)

	for direct in ["/plots/", "/plots/Efficiency/", "/plots/SEFD/", "/plots/Patterns/", "/data/"]:
		if not os.path.exists(results_path + direct):
			os.mkdir(results_path + direct)

	losses = []
	for AntPos in np.linspace(13,17.5, 2):

		change_list = [489]
		mod_line 	= ["  value            : %4.2f\n" % AntPos]
		edit_tor(model_file, out_file, change_list, mod_line) # <====== edit all tor files here


		command = "grasp-analysis batch.gxp out.out out.log"
		process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, cwd = GRASP_working_file)
		output, error = process.communicate()
		print("Done executing grasp")
		#print(output, error)

		#########  Process data files  ###########
		print("preparing to Plot")
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
	os.rename(results_path, results_path + "_minLoss = %4.2f" % minLoss)
	return minLoss

def loss(freq, effic):
	#calculate the loss funciton
	ans = 0
	n = 0
	for i, f in enumerate(freq):
		if (f > 30 and f < 90):
			ans -= effic[i]
			n+= 1
	return ans/n




if __name__ == '__main__':
	main()