import numpy as np
from copy import deepcopy
import os
import sys
import time
from datetime import datetime
import process_grasp
import subprocess
import scipy.optimize as op


'''
Class based python architecture to keep track of antenna dependencies and parameters

Todo:
x Release log after writing so it can be read during simulations
* Make setting of frequencies, zdists, n_freq easier

'''


class antenna(object):
	def __init__(self,
				 model_name = "",			 
				 parameter_names = ["z_dist"],
				 parameters = {"z_dist":16.6},
				 bounds = {"z_dist":[15.5,17]},
				 grasp_version = 10.3
				 ):
		self.model_name = model_name
		
		self.parameter_names = parameter_names
		self.parameters = parameters
		self.bounds = bounds
		self.grasp_version = grasp_version

		self.bool_gen_results_path_yet = False


	def __str__(self):
		return "Antenna"

	def set_global_directory_name(self, name =  "F:/Devin/Grasp/LWASandbox/"):
		self.global_directory_name = name

	def set_ticra_directory_name(self, name = "C:/Program Files/TICRA/"):
		self.ticra_directory_name = name

	def set_grasp_analysis_extension(self, name = ""):
		self.grasp_analysis_extension = name

	def gen_file_names(self):
		self.GRASP_working_file = self.global_directory_name + self.model_name + "/working/"
		self.in_tor_file = "bin/" + self.model_name + ".tor"
		self.out_tor_file = self.GRASP_working_file + self.model_name + ".tor"		

	def gen_specific_result_folder(self): #specific as in for a particular antenna configuration
		self.specific_result_folder_path = self.get_results_path() + "/" + self.get_datapoint_string()
		if not os.path.exists(self.specific_result_folder_path):
			os.mkdir(self.specific_result_folder_path)
			print ("Generated: " + self.get_datapoint_string(format_str = "%6.4f"))

	def get_specific_results_path(self):
		return self.get_results_path() + "/" + self.get_datapoint_string()

	def get_model_name(self):
		return self.model_name

	def get_bounds(self):
		return deepcopy(self.bounds)

	def get_results_path(self):
		if self.bool_gen_results_path_yet:
			return self.results_path
		else:
			self.gen_results_folder()
			return self.get_results_path()

	# def get_tor_line_numbers(self):
	# 	# print(self.tor_line_numbers)
	# 	return deepcopy(self.tor_line_numbers)

	def get_optimizable_parameter_names(self):
		names = self.get_parameter_names()
		for x in ["z_dist", "start_f", "end_f", "n_f"]:
			if x in names:
				names.remove(x)
		# print("OPT NAMES: ", names)
		return names

	def get_datapoint_string(self, format_str = "%4.3f"):
		'''
		Function which returns a string that describes the current antenna configuration.
		Often used for file generation or for printing to the terminal 
		'''
		#does NOT include / before string
		ans = ""

		# Dont want any of these names in the descriptor
		names = self.get_optimizable_parameter_names()

		# Assemble "local" datapoint descriptor
		for name in names:
			ans += name
			ans += "=" + format_str%self.parameters[name] + "_"

		if ans[-1] == '_':
			ans = ans[:-1] #remove last char if it is "_"
		return ans

	def get_parameters(self):
		return deepcopy(self.parameters)

	def get_parameter_names(self):
		return deepcopy(self.parameter_names)

	def set_parameters(self, new_parameters):
		self.parameters.update(new_parameters)

	def set_number_of_focal_lengths(self, number):
		self.number_of_z_dists = number

	def get_error_intersection(self):
		error = 0
		for name in self.bounds:
			if self.parameters[name] < self.bounds[name][0]:
				error += (self.parameters[name] - self.bounds[name][0])**2
			elif self.parameters[name] > self.bounds[name][1]:
				error += (self.parameters[name] - self.bounds[name][1])**2
		# print("error bounds: ", error)
		return error

	def _create_file_w_timestamp(self, location = os.getcwd()):
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

	def gen_results_folder(self):
		if self.bool_gen_results_path_yet == False:
			results_path = "../Results/"
			results_path += self._create_file_w_timestamp(results_path)
			self.results_path = results_path
			self.bool_gen_results_path_yet = True


	def gen_sub_folders(self, plot_feed = False):
		for directory in ["/plots/", "/plots/Efficiency/", "/plots/SEFD/", "/plots/Patterns/", "/data/"]:
			if not os.path.exists(self.get_specific_results_path() + directory):
				os.mkdir(self.get_specific_results_path() + directory)
		if plot_feed:
			feed_dir = "/plots/FeedPatterns/"
			if not os.path.exists(self.get_specific_results_path() + feed_dir):
				os.mkdir(self.get_specific_results_path() + feed_dir)

	def init_global_file_log(self):
		self.log_name = self.get_results_path() +"/log.csv"
		self.log = open(self.log_name, 'w+')
		self.log.write(self.__str__() + "\n") #write name of antenna here
		self.log.write("GRASP Version: " + str(self.grasp_version) + "\n")
		self.log.close()
		self._gen_global_log_headers()
		

	# def _gen_global_log_headers(self):
	# 	self.log = open(self.log_name, 'a')
	# 	self.log.write("Efficiency, Loss\n")
	# 	self.log.close()

	def _gen_global_log_headers(self):
		self.log = open(self.log_name, 'a')
		names = self.get_optimizable_parameter_names()
		# print("gen log head", names)
		self.log.write("z_dist,")
		for x in names:
			self.log.write(x +",")
		self.log.write("Efficiency, Loss\n")
		self.log.close()

	def edit_msh(self):
		return 0

	def edit_tor(self): #template 
		'''
		Edits the GRASP TOR file
		(i.e. everything but the antenna msh file)

		Updated version which does not rely on manually entering line numbers
		'''
		print("Opening TOR Files")
		print(self.in_tor_file, self.out_tor_file)
		f = open(self.in_tor_file,'r')
		g = open(self.out_tor_file,'w+')


		# Generate list of things that 
		change_list = self.get_parameter_names()
		msh_items = ["x", "y", "z"]
		for x in msh_items:
			change_list.remove(x)
		# print(change_list)
		# modified_line = self.get_tor_line_replacement()

		line = f.readline()
		while line: #iterate through all lines in INPUT
			g.write(line) # write same line to output
			if "real_variable" in line: #check if keyword in line
				for key in change_list: #iterate through keys
					if key in line: #check if key is also in line
						g.write(f.readline()) # read next line, should be "("
						f.readline()
						write_string = "  value            : %7.5f\n" % self.parameters[key]
						g.write(write_string) #write modified line to file

			line = f.readline()


		g.close()
		f.close()
		print("Done writing TOR")

	def exeGRASP(self):
		'''
		Executes GRASP processor
		'''
		print("EXECUTING GRASP")
		if (self.grasp_version == 10.6):
			command = [self.ticra_directory_name + "GRASP-10.6.0/bin/grasp-analysis" +self.grasp_analysis_extension, "batch.gxp", "out.out", "out.log"]
			print("Using Version 10.6.0")
		elif (self.grasp_version == 10.3):
			command = [self.ticra_directory_name + "GRASP-10.3.1/bin/grasp-analysis"+self.grasp_analysis_extension, "batch.gxp", "out.out", "out.log"]
			print("Using Version 10.3.1")
		else:
			command = ["grasp-analysis", "batch.gxp", "out.out", "out.log"]
			print("Using Version from PATH")
		print ("command: ", command)
		sys.stdout.flush()

		process = subprocess.Popen(command, stdout=subprocess.PIPE, cwd = self.GRASP_working_file)
		output, error = process.communicate()
		# print ("GRASP commuique: ", output, error)
		print("Done EXECUTING GRASP")
		sys.stdout.flush()


	def process_data_files(self, plot_feed = False):
		'''
		Collects all data from GRASP result files, plots everything
		'''

		## COLLECT ALL DATA FROM GRASP GENERATED FILES
		freq, s11 = process_grasp.process_par(self.GRASP_working_file + "S_parameters.par")
		dmax, cut = process_grasp.process_cut(self.GRASP_working_file + "Field_Data.cut", freq)
		if plot_feed:
			_ , feed_cut = process_grasp.process_cut(self.GRASP_working_file + "Feed_Data.cut", freq)

		## CALCULATE RELEVANT FIGURES OF MERIT
		mis = process_grasp.calc_mismatch(s11)
		aperture = process_grasp.calc_app_eff(freq, dmax)
		tsys = process_grasp.Tsys(freq)
		SEFD = process_grasp.SEFD(freq, dmax)

		## CALCULATE EFFICIENCIES AND LOSSES
		efficiency, mis_loss = self._loss(freq, aperture, mis)
		loss_val = efficiency + mis_loss
		self.EF.append(efficiency)
		self.LOSS.append(loss_val)

		## DEFINE CONVENIENT VARIABLE (REUSE CODE FROM OPT_FULL.PY)
		AntPos = self.parameters["z_dist"]

		## GENERATE FILE NAMES
		plots_directory = self.get_specific_results_path() + "/plots/"
		data_directory = self.get_specific_results_path() + "/data/"

		# print("NEW DIRECTORIES", plots_directory, data_directory)

		pattern_directory = plots_directory + "Patterns/AntPos=%4.2f_loss=%4.2f/"% (AntPos, efficiency)
		if plot_feed:
			feed_directory = plots_directory + "FeedPatterns/AntPos=%4.2f_loss=%4.2f/"% (AntPos, efficiency)

		## MAKE FILES
		if not os.path.exists(pattern_directory):
			os.mkdir(pattern_directory)

		if plot_feed:
			if not os.path.exists(feed_directory):
				os.mkdir(feed_directory)

		## PLOT DATA
		process_grasp.plot_pair_efficiencies(freq, s11, dmax, plots_directory + "Efficiency/Efficiencies Antenna Position = %4.2f Ef = %4.2f Loss = %4.3f.png" % (AntPos, efficiency, loss_val), AntPos)
		process_grasp.plot_SEFD(freq, dmax, plots_directory + "SEFD/SEFD Antenna Position = %4.2f Ef = %4.2f Loss = %4.3f.png" % (AntPos, efficiency, loss_val), AntPos)
		for frequency in freq:
			process_grasp.plot_cut(frequency, cut, AntPos, pattern_directory +"Radiation_Pattern_Position=%4.2f_Freq=%4.2f_Ef=%4.2f_Loss=%4.2f.png" %(AntPos, frequency, efficiency, loss_val))
			if plot_feed:
				process_grasp.plot_cut(frequency, feed_cut, AntPos, feed_directory +"Feed_Pattern_Position=%4.2f_Freq=%4.2f_Ef=%4.2f_Loss=%4.2f.png" %(AntPos, frequency, efficiency, loss_val), feed_pattern = True)
			# process_grasp.plot_cut(frequency, cut, AntPos, pattern_directory +"Radiation Pattern Position = %4.2f Freq = %4.2f Ef = %4.2f Loss = %4.2f.png" %(AntPos, frequency, efficiency, loss_val))

		## WRITE DATA TO FILES
		local_log = open(data_directory+"Results AntPos=%4.2f Ef = %4.2f Loss=%4.3f.csv" %(AntPos, efficiency, loss_val), 'w+')
		local_log.write(time.strftime("%c"))
		local_log.write("\n")
		local_log.write("Ef = %6.4f Loss = %6.4f.png\n" %(efficiency, loss_val))
		local_log.write("Frequency [MHz], S11 [dB], Dmax [dBi], Mismatch Efficiency, Aperture Efficiency, Tsys [1E3 K], SEFD [1E6 Jy]\n")
		for ii in range(len(freq)):
			local_log.write("%6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.4f\n" % (freq[ii], s11[ii], dmax[ii], mis[ii], aperture[ii], tsys[ii]/1E3, SEFD[ii]/1E6))
		local_log.close()

		## WRITE DATA TO LOG.CSV
		param = [self.parameters["z_dist"]]
		names = self.get_optimizable_parameter_names()

		for name in names:
			param.append(self.parameters[name])
		param.append(efficiency)
		param.append(loss_val)
		template = "%7.4f, "*len(param) + '\n'
		try:
			self.log = open(self.log_name, 'a')
			self.log.write(template % tuple(param))
			self.log.close()
		except:
			print("could not write to log.csv: ", template%tuple(param))


	def _loss(self, freq, effic, miss):
		#calculate the loss funciton
		ans = 0
		n = 0
		miss_loss = 0
		for i, f in enumerate(freq):
			if (f >= 60 and f <= 80):
				ans -= effic[i]
				n+= 1
				if (miss[i] < .4):
					miss_loss += (.4 - miss[i])
		return ans/n, miss_loss/n

	def wrap_up_simulation(self, ef, minloss):
		'''
		Does the following things:
		1. Changes name of specific results folder (i.e. folder for one specific configuration)
		2. writes configuration to global log file
		'''
		try:
			os.rename(self.specific_result_folder_path, self.specific_result_folder_path + "_EF = %4.3f_minLoss = %4.3f" % (ef, minloss))
		except:
			print("Error: Cannot rename file %s" %self.specific_result_folder_path)
	


	def simulate_single_configuration(self, parameters, parameter_names, plot_feed = False):
		'''
		must take in a vector of antenna parameters
		and return a loss function
		'''

		new_parameters = {}
		for ii, name in enumerate(parameter_names):
			# print (ii, name, new_parameters, parameters, parameter_names)
			new_parameters[name] = parameters[ii]
		self.set_parameters(new_parameters)
		self.gen_specific_result_folder()
		error = self.get_error_intersection()

		if error > 0.000000001:
			self.wrap_up_simulation(error, error)
			return error

		self.gen_sub_folders(plot_feed)
		self.edit_msh()

		self.EF = []
		self.LOSS = []
		for AntPos in np.linspace(self.bounds["z_dist"][0], self.bounds["z_dist"][1], self.number_of_z_dists):
			self.parameters["z_dist"] = AntPos
			self.edit_tor()
			self.exeGRASP()
			self.process_data_files(plot_feed)

		minLoss = np.min(self.LOSS)
		best_ef = self.EF[self.LOSS.index(minLoss)]

		self.wrap_up_simulation(best_ef, minLoss)

		# try:
		# 	os.rename(self.specific_result_folder_path, self.specific_result_folder_path + "_EF = %4.3f_minLoss = %4.3f" % (best_ef, minLoss))
		# except:
		# 	print("Error: Cannot rename file %s" % self.specific_result_folder_path)
		# template = "%7.4f, "*9 + '\n'
		# param = []
		# for name in self.parameter_names:
		# 	if name != "z_dist"
		# 		param.append(self.parameters[name])
		# param.append(best_ef)
		# param.append(minLoss)
		# self.Log.write(template % tuple(param))
		return minLoss
		

class LWA_like(antenna):
	def __init__(self, x = .77, y = .16, z = -.01,
				start_f = 60.0, end_f = 80.0, n_f = 5, alpha = 0,
				bnd_x = [0, 1.5], bnd_y = [0, .75], bnd_z = [-3.5, 1], 
				bnd_start_f = [0,100], bnd_end_f = [0,100], bnd_n_f =[1,100], bnd_alpha = [0, 360],
				grasp_version = 10.3):
		
		antenna.__init__(self, grasp_version = grasp_version)
		self.model_name = "40mLWA104"

		self.parameter_names += ["x", "y", "z", "start_f", "end_f", "n_f", "alpha"]
		self.parameters.update({"x":x, "y":y, "z":z, "start_f":start_f, "end_f":end_f, "n_f":n_f, "alpha":alpha })
		self.bounds.update({"x":bnd_x,"y":bnd_y, "z":bnd_z, "start_f":bnd_start_f, "end_f":bnd_end_f,"n_f":bnd_n_f, "alpha":bnd_alpha})
		# self.tor_line_numbers = {"z_dist":339, "start_f":346, "end_f":351, "n_f":356,"alpha":361} #checked
		
	def __str__(self):
		return "LWA like feed no Directors"



	def edit_msh(self):
		msh_template = "patch_leaf_inverted_V.msh"
		msh_out = self.GRASP_working_file + msh_template

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
			diff = np.array([self.parameters["x"], self.parameters["y"], self.parameters["z"]])
			if pt[0] == 0:
				diff[1] *= -1
			elif pt[0] == 1:
				diff[1] = 0
			diff *= pt[1]
			modified_line.append("%6.3f  %6.3f  %6.3f\n" % tuple(diff+origin[pt[0]]))

		f = open("bin/" + msh_template,'r')
		g = open(msh_out,'w+')

		for i, line in enumerate(f): #0-indexed line number
			if i in change_list:
				g.write(modified_line[change_list.index(i)])
			else:
				g.write(line)

		g.close()
		f.close()
		print("Done writing MSH")

class ELfeed(LWA_like):
	def __init__(self, sp = 1.2, start_f = 60.0, end_f = 80.0, n_f = 5, alpha = 0,
				bnd_sp = [0, 1.5], grasp_version = 10.3): #seperation is half the distance between dipoles
		
		LWA_like.__init__(self, start_f = start_f, end_f = end_f, n_f = n_f, alpha = 0, grasp_version = grasp_version)
		self.model_name = "40mQuadDipole"

		self.parameter_names += ["sp"]
		self.parameters.update({"sp":sp})
		self.bounds.update({"sp":bnd_sp})
		
		# self.tor_line_numbers = {"z_dist":489, "sp":333, "start_f":474, "end_f":479, "n_f":484,"alpha":510} #checked
		
	def __str__(self):
		return "Eleven Feed with no Directors"

	def get_error_ant_intersection(self):
		ans =  self.parameters["sp"] - self.parameters["x"] - self.parameters["y"] - 0.08 - 0.012
		if ans < 0:
			print("antenna antenna Instersection: ", ans)
			return ans**2
		else:
			return 0

	def get_error_intersection(self):
		return LWA_like.get_error_intersection(self) + self.get_error_ant_intersection()

	# def _gen_global_log_headers(self):
	# 	self.log = open(self.log_name, 'a')
	# 	self.log.write("antenna separation,antenna x,antenna y,antenna z,Efficiency,Loss\n")
	# 	self.log.close()

class ELfeedDir(ELfeed):
	def __init__(self, dl = 1.2, dw = .482, dsep = .25,
				start_f = 60.0, end_f = 80.0, n_f = 5, alpha = 0,
				bnd_dl = [0, 2.5], bnd_dw = [0, .75], bnd_dsep = [0, 1.5], #  positive dir_sep vals are directors
				grasp_version = 10.3): 						   #  Negative dir_sep vals are reflectors
		
		ELfeed.__init__(self,start_f = start_f, end_f = end_f, n_f = n_f, alpha = 0, grasp_version = grasp_version)
		self.model_name = "40mQuadDipoleWDir"

		self.parameter_names += ["dl", "dw", "dsep"]
		self.parameters.update({"dl":dl, "dw":dw, "dsep":dsep})
		self.bounds.update({"dl":bnd_dl, "dw":bnd_dw, "dsep":bnd_dsep})
		
		# self.tor_line_numbers = {"z_dist":489, "sp":333, "dl":502, "dw":507, "dsep":512, "start_f":457, "end_f":480, "n_f":485,"alpha":547} #checked
		#checked
	def __str__(self):
		return "Eleven Feed with one Director (dsep > 0)"

	def get_error_dir_intersection(self):
		ans =  2*self.parameters["sp"] - self.parameters["dw"] - self.parameters["dl"]
		if ans < 0:
			print("director Instersection: ", ans)
			return ans**2
		else:
			return 0

	def get_error_intersection(self):
		return ELfeed.get_error_intersection(self) + self.get_error_dir_intersection()

	def _gen_global_log_headers(self):
		self.log = open(self.log_name, 'a')
		self.log.write("antenna separation,antenna x,antenna y,antenna z,director length,director width,director sep,Efficiency,Loss\n")
		self.log.close()

class ELfeedRef(ELfeedDir):
	def __init__( dsep = -1.5, bnd_dsep = [-3.05, 0], #  positive dir_sep vals are directors
				start_f = 60.0, end_f = 80.0, n_f = 5, alpha = 0, grasp_version = 10.3): 	  #  Negative dir_sep vals are reflectors

		ELfeedDir.__init__(self, start_f = start_f, end_f = end_f, n_f = n_f, alpha = 0, grasp_version = grasp_version)
		self.model_name = "40mQuadDipoleWDir"

		self.parameters.update({"dsep":dsep})
		self.bounds.update({"dsep":bnd_dsep})
		

	def __str__(self):
		return "Eleven Feed with one Reflector (dsep < 0)"

	def get_error_dir_ant_intersection(self):
		ans =  self.parameters["z"] - self.parameters["dsep"]
		if ans < 0:
			return ans**2
		else:
			return 0

	def get_error_intersection(self):
		return ELfeedRef.get_error_intersection(self) +  self.get_error_dir_ant_intersection()

	def _gen_global_log_headers(self):
		self.log = open(self.log_name, 'a')
		self.log.write("antenna separation,antenna x,antenna y,antenna z,reflector length,reflector width,reflector sep,Efficiency,Loss\n")
		self.log.close()





















