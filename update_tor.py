import numpy as np 
import matplotlib.pyplot as plt 
import process_grasp
import subprocess
import time
from datetime import datetime
import pandas as pd
import os

GRASP_working_file = "F:/Devin/Grasp/LWASandbox/40mLWA/working/"
model_file = "model_tor.txt"
out_file = GRASP_working_file + "40mLWA.tor"

def edit_tor(model, out, z):
	print("Opening Files")
	f = open(model,'r')
	g = open(out,'w+')

	for i in range(8):
		g.write(f.readline())

	g.write("  origin           : struct(x: 0.0 m, y: 0.0 m, z: %4.2f m),\n" % z)
	f.readline()

	for i in f:
		g.write(i)

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

#########  Generate folders  ###########
results_path = "../Results/"
results_path += create_file_w_timestamp(results_path)
print(results_path)

for direct in ["/plots/", "/plots/Efficiency/", "/plots/SEFD/", "/plots/patterns/"]:
	if not os.path.exists(results_path + direct):
		os.mkdir(results_path + direct)

for z in [16]:#np.linspace(11,17.5, 50):
	direct = results_path + "/plots/patterns/z-%4.2f/"%z
	if not os.path.exists(direct):
		os.mkdir(direct)

	edit_tor(model_file, out_file, z)

	command = "grasp-analysis batch.gxp out.out out.log"
	process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, cwd = GRASP_working_file)
	output, error = process.communicate()
	print("Done executing grasp, Output:\n")
	#print(output, error)


	#########  Process data files  ###########
	print("\npreparing to Plot")
	freq, s11 = process_grasp.process_par(GRASP_working_file + "S_parameters.par")
	dmax, cut = process_grasp.process_cut(GRASP_working_file + "Field_Data.cut", freq)
	mis = process_grasp.calc_mismatch(s11)
	aperture = process_grasp.calc_app_eff(freq, dmax)
	tsys = process_grasp.Tsys(freq)
	SEFD = process_grasp.SEFD(freq, dmax)
	print("Dmax = ", dmax)
	process_grasp.plot_pair_efficiencies(freq, s11, dmax, "Efficiencies Antenna Position = %4.2f.png" % z, z)
	process_grasp.plot_SEFD(freq, dmax, "SEFD Antenna Position = %4.2f.png" % z, z)

	x = open("../data/Results z=%4.2f.csv" %z, 'w+')
	x.write(time.strftime("%c"))
	x.write("\n")
	x.write("Frequency [MHz], S11 [dB], Dmax [dBi], Mismatch Efficiency, Aperture Efficiency, Tsys [1E3K], SEFD [1E6Jy]\n")
	for ii in range(len(freq)):
		x.write("%6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f\n" % (freq[ii], s11[ii], dmax[ii], mis[ii], aperture[ii], tsys[ii]/1E3, SEFD[ii]/1E6))
	x.close()