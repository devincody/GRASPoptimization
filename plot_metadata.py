import pandas as pd 
import numpy as np 
import os
import matplotlib.pyplot as plt
from collections import defaultdict

def conv(x):
	return float(x.strip())

metadata_path = "../metadata"
if not os.path.exists(metadata_path):
	os.mkdir(metadata_path)
os.chdir(metadata_path)

units = defaultdict(lambda:'m', {'alpha':'deg', 'theta':'deg', 'efficiency':'%'})

cwd = os.getcwd()
items = os.listdir(cwd)
for i in items:
	if ".csv" in i[-5:]:
		print("Processing file ", i)
		a = pd.read_csv(i, skiprows = 2, header = 0, index_col=False)

		keys = a.keys()
		# print keys

		try:
			if "Unnamed" in keys[-1]:
				a = a.drop(keys[-1], axis = 1)
				print("dropping last column", keys[-1])
		except:
			pass

		try:
			loss = ' Loss'
			a[loss]
		except:
			loss = 'Loss'
			a[loss]

		a["radius"] = np.sqrt(a['x']**2 + a['z']**2)
		a["theta"] = np.arctan2(a['z'],a['x'])*180.0/np.pi

		a[loss] *= -1

		path = i[:-4]
		if not os.path.exists(path):
			os.mkdir(path)

		parameters = list(a.keys())
		parameters.remove(loss)
		for p in parameters:
			# print p
			if not "Unnamed:" in p:
				fig, ax = plt.subplots()
				ax.scatter(a[p],a[loss], alpha = .8, c = a['alpha']/45, cmap = 'seismic')
				ax.set_xlabel("%s [%s]" % (p, units[p.lower()]))
				ax.set_ylabel('Loss')
				ax.grid(linewidth = 1, linestyle = '--')
				ax.set_axisbelow(True)
				ax.set_title("Loss vs. %s"% p)
				plt.savefig(path + '/' + p + '.png')
				plt.close()

		cols = a.columns.tolist()
		print (cols)
		
		#reorder loss and eff with radius and theta, better formatting
		cols.remove(loss) 
		cols.remove("Efficiency")
		cols.append("Efficiency")
		cols.append(loss)

		a = a[cols]
		# a.sort_values([loss], inplace = True, ascending = False)
		a = a.nlargest(20,[loss], keep='first')
		a.to_csv("%s/%s_Top20.csv"%(path,path))


