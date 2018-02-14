import numpy as np 
import matplotlib.pyplot as plt 


def main():
	PA = "F:/Devin/Grasp/LWA Sandbox/40mLWA/Job_26/"
	freq, s11 = process_par(PA + "Sparameters.par")
	dmax = process_cut(PA + "FieldData.cut", freq)
	plot_pair_efficiencies(freq, s11, dmax, "Efficiencies.png")


def process_par(f_name):
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
	while ("Field data" in line):
		line = f.readline()
		line = line.split()
		line = [float(x) for x in line]
		print ("Headder = ", line)
		for ii in range(int(line[2])):
			angle=line[0] + ii*line[1]
			fields = [float(x) for x in f.readline().split()]
			co = 10*np.log10(fields[0]**2+fields[1]**2)
			cx = 10*np.log10(fields[2]**2+fields[3]**2)
			if ii == 100:
				dmax.append(np.max([co, cx]))
		line = f.readline() #should be Field data... if more data

	dmax_f = np.zeros(len(freq))
	for i in range(len(freq)):
		dmax_f[i] = np.max(dmax[i*3:(i+1)*3])

	return dmax_f

def calc_mismatch(s11): 
	#s11 in dB
	#returns %
	return (1-10**(s11/10))*100

def calc_app_eff(freq, dmax): 
	#freq in MHz, dmax in dB
	#returns %
	c = 3E8#meters/sec
	r = 20.0#meters
	wave_len = c/(freq*1E6)
	aphy = np.pi*r**2
	dmax_lin = 10**(dmax/10)
	return dmax_lin*wave_len**2 /(4*np.pi*aphy)*100


def plot_pair_efficiencies(freq, s11, dmax, title):
	fig, ax1 = plt.subplots()
	ax2 = ax1.twinx()
	ms = ax1.plot(freq,calc_mismatch(s11), 'b', label= "Mismatch Efficiency")
	ax1.set_title("Mismatch and Aperture Efficiencies vs Frequency")
	ax1.set_xlabel("Frequency [MHz]")
	ax1.set_ylabel("Mismatch Efficiency[%]")

	ap = ax2.plot(freq,calc_app_eff(freq, dmax), 'r', label = "Aperture Efficiency")
	ax2.set_ylabel("Aperture Efficiency [%]")

	lns = ms+ap
	labs = [x.get_label() for x in lns]
	ax1.legend(lns, labs)

	# ax1.legend()
	plt.savefig("../plots/" + title)



if __name__ == '__main__':
	main()





# print(calc_missmatch(-.46))
# print(calc_app_eff(50, 21.32))