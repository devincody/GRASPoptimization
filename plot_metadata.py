import pandas as pd 
import numpy as np 
import os
import matplotlib.pyplot as plt

metadata_path = "../metadata"
if not os.path.exists(metadata_path):
	os.mkdir(metadata_path)
os.chdir(metadata_path)

cwd = os.getcwd()
items = os.listdir(cwd)
for i in items:
	if ".csv" in i[-5:]:
		a = pd.read_csv(i, skiprows = 2)
		loss = ' Loss'
		a["radius"] = np.sqrt(a['x']**2 + a['z']**2)
		a["theta"] = np.arctan2(a['z'],a['x'])*180.0/np.pi
		a[loss] *= -1

		path = i[:-4]
		if not os.path.exists(path):
			os.mkdir(path)

			loss = ' Loss'
			parameters = list(a.keys())
			parameters.remove(loss)
			for p in parameters:
				print p
				if not "Unnamed:" in p:
					fig, ax = plt.subplots()
					ax.scatter(a[p],a[loss], alpha = .8, c = a['alpha']/45, cmap = 'seismic')
					ax.set_xlabel(p)
					ax.set_ylabel('Loss')
					plt.savefig(path + '/' + p + '.png')
					plt.close()


