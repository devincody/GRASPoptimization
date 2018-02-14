import numpy as np 
import matplotlib.pyplot as plt 
import process_grasp
import subprocess

GRASP_working_file = "F:/Devin/Grasp/LWASandbox/40mLWA/working/"
model_file = "model_tor.txt"
out_file = GRASP_working_file + "40mLWA.tor"

for z in [4, 6, 8, 10, 12, 14, 14.5, 15, 15.5, 16]:

	def edit_tor(model, out, z):
		print("Opening Files")
		f = open(model,'r')
		g = open(out,'w+')


		for i in range(8):
			g.write(f.readline())

		g.write("  origin           : struct(x: 0.0 m, y: 0.0 m, z: %3.1f m),\n" % z)
		f.readline()

		for i in f:
			g.write(i)

		g.close()
		f.close()
		print("Done writing TOR")


	edit_tor(model_file, out_file, z)

	command = "grasp-analysis batch.gxp out.out out.log"
	process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, cwd = GRASP_working_file)
	output, error = process.communicate()
	print("Done executing grasp, Output:\n")
	print(output, error)


	print("\npreparing to Plot")
	freq, s11 = process_grasp.process_par(GRASP_working_file + "S_parameters.par")

	dmax = process_grasp.process_cut(GRASP_working_file + "Field_Data.cut", freq)
	mis = process_grasp.calc_mismatch(s11)
	aperture = process_grasp.calc_app_eff(freq, dmax)
	tsys = process_grasp.Tsys(freq)
	SEFD = process_grasp.SEFD(freq, dmax)
	print("Dmax = ", dmax)
	process_grasp.plot_pair_efficiencies(freq, s11, dmax, "Efficiencies Antenna Position = %3.1f.png" % z, z)

	x = open("../data/Results z=%3.1f.csv" %z, 'w+')
	x.write("frequency [MHz], S11 [dB], Dmax [dBi], Mismatch Efficiency [%], Aperture Efficiency [%], Tsys [K], SEFD [Jy]\n")
	for ii in range(len(freq)):
		x.write("%f, %f, %f, %f, %f, %f, %f\n" % (freq[ii], s11[ii], dmax[ii], mis[ii], aperture[ii], tsys[ii], SEFD[ii]))
	x.close()