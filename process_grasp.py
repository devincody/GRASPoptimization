import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd
import time

def main():
	PA = "F:/Devin/Grasp/LWASandbox/40mLWA/Job_26/"
	freq, s11 = process_par(PA + "Sparameters.par")
	dmax, cut = process_cut(PA + "FieldData.cut", freq)
	plot_pair_efficiencies(freq, s11, dmax, "Efficiencies.png", 16)


def process_par(f_name): #Process S parameters document return  1D numpy arrays of frequencies and s11s
	f = open(f_name)
	f.readline() #ignore header
	freq = []
	s11 = []
	for line in f:
		line = line.split()
		freq.append(float(line[0])*1E3) #grab freq, convert to MHz
		s11.append(float(line[1]))
	f.close()
	# print("freq = ", freq)
	return np.array(freq), np.array(s11)


def process_cut(f_name, freq):
	f = open(f_name)
	line = f.readline()
	dmax = []
	i = -1 # Frequency Index
	cut = pd.DataFrame()
	while ("Field data" in line):
		line = f.readline()
		line = line.split()
		line = [float(x) for x in line]

		'''
		line: Headder defined as per the grasp manual:
		0: starting angle
		1: angle delta
		2: number of angles
		3: phi cut
		4-6: Various parameters related to polarization and radiation pattern data representation
		'''

		# print ("Headder = ", line)
		phi = line[3] #should be 3 possible values: 0, 45, 90
		if (phi == 0):
			i += 1
		# print(line, i, freq, f_name)

		frequency = freq[i]
		series_name_co = "f%4.2f:p%4.2f:co" %(frequency, phi)
		series_name_cx = "f%4.2f:p%4.2f:cx" %(frequency, phi)

		angles = []
		dbi_co = []
		dbi_cx = []

		for ii in range(int(line[2])):
			angle= -180. + ii*line[1] ####- 180. #180 used because antenna is technically upside down
			angles.append(angle)
			fields = [float(x) for x in f.readline().split()]
			cx = 10*np.log10(fields[0]**2 + fields[1]**2)
			co = 10*np.log10(fields[2]**2 + fields[3]**2)

			dbi_co.append(co)
			dbi_cx.append(cx)
			if ii == 100:
				dmax.append(np.max([co, cx]))

		cut[series_name_co] = dbi_co
		cut[series_name_cx] = dbi_cx
		line = f.readline() #should be Field data... if more data
	cut["angles"] = angles

	dmax_f = np.zeros(len(freq))
	for i in range(len(freq)):
		dmax_f[i] = np.max(dmax[i*3:(i+1)*3])

	return dmax_f, cut

def calc_mismatch(s11): 
	#s11 in dB
	#returns %
	return (1-10**(s11/10))

def calc_app_eff(freq, dmax): 
	#freq in MHz, dmax in dB
	#returns %
	c = 2.997E8#meters/sec
	r = 20.0#meters
	wave_len = c/(freq*1E6)
	aphy = np.pi*r**2
	dmax_lin = 10**(dmax/10)
	return dmax_lin*wave_len**2 /(4*np.pi*aphy)


def plot_pair_efficiencies(freq, s11, dmax, location, z):
	plt.figure()
	plt.plot(freq,calc_mismatch(s11), 'b', label= "Mismatch Efficiency")
	plt.title("Mismatch and Aperture Efficiencies z = %4.2f" % z)
	plt.xlabel("Frequency [MHz]")
	plt.ylabel("Efficiency")
	plt.ylim([0, 1])
	plt.yticks(np.arange(0, 1.01,.1))

	plt.plot(freq,calc_app_eff(freq, dmax), 'r', label = "Aperture Efficiency")
	# plt.set_ylabel("Aperture Efficiency [%]")
	# lns = ms+ap
	# labs = [x.get_label() for x in lns]
	# plt.legend(lns, labs)
	plt.legend()
	plt.savefig(location)

def plot_SEFD(freq, dmax, location, z):
	plt.figure()
	plt.semilogy(freq,SEFD(freq, dmax), 'b', label= "SEFD")
	plt.title("SEFD at z = %4.2f" % z)
	plt.xlabel("Frequency [MHz]")
	plt.ylabel("SEFD")
	plt.ylim([1, 5E5])

	plt.legend()
	# print(location)
	plt.savefig(location)


def plot_cut(frequency, cut, z, title, feed_pattern = False):
	#freq is which frequency to use
	#
	plt.rc('axes', linewidth=2)
	fig, ax = plt.subplots(1,3, figsize=(30,10)) #each of the phis

	for i, phi in enumerate([0, 45, 90]):
		axi = ax[i]
		series_name_co = "f%4.2f:p%4.2f:co" %(frequency, phi)
		series_name_cx = "f%4.2f:p%4.2f:cx" %(frequency, phi)


		data = np.copy(cut[series_name_co])


		if feed_pattern:
			# For Feed Patterns
			half = int(len(data)/2)
			temp = np.copy(data[half:])
			temp1 = np.copy(data[:half])
			data[len(temp):] = temp1
			data[:len(temp)] = temp

		# print(cut["angles"])
		axi.plot(cut["angles"] ,data, 'b', label= "co")


		# axi.plot(cut["angles"]-180,cut[series_name_cx], 'r', label= "cx")
		axi.set_title("$\phi$ = %4.2f" %phi, fontsize = 20)
		axi.legend()
		axi.set_ylim([-20,30])
		axi.set_xlim([-180,180])
		axi.set_xticks(range(-180, 181, 30))
		axi.grid(linewidth = 2, linestyle = '--')
		axi.set_yticks(range(-20,31,5))

		fontsize = 14
		for tick in axi.xaxis.get_major_ticks():
		    tick.label1.set_fontsize(fontsize)
		    # tick.label1.set_fontweight('bold')
		for tick in axi.yaxis.get_major_ticks():
		    tick.label1.set_fontsize(fontsize)
		    # tick.label1.set_fontweight('bold')
		if i == 0:
			axi.set_ylabel("Amplitude [dBi]", fontsize = 16)

		if i == 1:
			axi.set_xlabel("Angle [Degrees]", fontsize = 16)
			#axi.set_title()
			fig.suptitle(r"Radiation Pattern at z = %4.2f, $\nu$ = %4.2f" % (z,frequency), fontsize=25)

	# 
	# fig.show()
	# time.sleep(1)
	# print(title)
	fig.savefig(title)
	plt.close('all')
	plt.rc('axes', linewidth=1)


def Tsys(freq):
	return 300*(150/freq)**2.5 + 600


def SEFD(freq, dmax):
	k = 1.38E3 #jy m^2 s k^-1
	r = 20#m
	return 2*k*Tsys(freq)/(calc_app_eff(freq,dmax)*np.pi*r**2)




if __name__ == '__main__':
	main()





# print(calc_missmatch(-.46))
# print(calc_app_eff(50, 21.32))