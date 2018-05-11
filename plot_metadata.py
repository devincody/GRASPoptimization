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

cwd = os.getcwd()
items = os.listdir(cwd)
for i in items:
	if ".csv" in i[-5:]:
		print("item: ", i)
		a = pd.read_csv(i, skiprows = 2, header = 0, index_col=False)

		keys = a.keys()
		# print(keys)
		# print(a)

		if "Unnamed" in keys[-1] and " " in a[k[-1]][1]:
			a = a.drop(keys[-1])
			print("dropping last column")

		try:
			loss = ' Loss'
			a[loss]
		except:
			loss = 'Loss'
			a[loss]
		a["radius"] = np.sqrt(a['x']**2 + a['z']**2)
		a["theta"] = np.arctan2(a['z'],a['x'])*180.0/np.pi
		# print(type(a))
		# print(type(a[loss]))
		# print(type(a.loc[:,loss]))
		# print(type(a.loc[1,loss]))
		# a[loss] = a[loss].multiply(-1.0)
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
				ax.set_xlabel(p)
				ax.set_ylabel('Loss')
				plt.savefig(path + '/' + p + '.png')
				plt.close()


