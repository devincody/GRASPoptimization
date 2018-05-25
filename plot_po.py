import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt 
import os

start_f = 0.1818
end_f = 2
n_f = 20

po_path = "../po_opt"
if not os.path.exists(po_path):
	os.mkdir(po_path)

os.chdir(po_path)


cwd = os.getcwd()
items = os.listdir(cwd)
for it in items:
	if ".csv" in it[-5:]:
		print("Processing file ", it)

		path = it[:-4]
		if not os.path.exists(path):
			os.mkdir(path)

		a = pd.read_csv("log.csv", skiprows=2, index_col=False)
		a.dropna(inplace = True) #get rid of nans
		a["Efficiency"] *= -100
		z_dists = a["z_dist"].unique()
		number_of_z_dists = len(z_dists)
		# Append frequency information to dataFrame
		freq = np.linspace(start_f,end_f,n_f)
		for i in range(len(a)):
			a.at[i, "f"] = freq[int(i/number_of_z_dists)]
		labels = [] #keep track of good z_dists
		ax = plt.axes()
		good_dists = pd.DataFrame(columns = a.keys())
		for i in z_dists:
			# For each z_dist, create smaller dataFrame with only the one z_dist
			one_z_dist = a.where(a["z_dist"]==i)
			one_z_dist.dropna(inplace = True)
			if (np.mean(one_z_dist["Efficiency"]) > 40):
				# If sufficiently high sucess metric, plot the configuration
				one_z_dist.plot(x="f", y="Efficiency", ax = ax)
				labels.append("%6.3f"%i)
				good_dists = one_z_dist.append(good_dists)
		# ax.legend()
		ax.legend(labels)
		ax.set_title("Efficiency vs. Frequency at several Focus Lengths")
		ax.set_xlabel("Frequency [GHz]")
		ax.set_ylabel("Efficiency [%]")
		plt.savefig(path + "/" + path + ".png")

		keys = good_dists.keys()
		try:
			if "Unnamed" in keys[-2]:
				good_dists = good_dists.drop(keys[-2], axis = 1)
				print("dropping unnamed column", keys[-2])
		except:
			pass

		# good_dists.sort_values(by="f", inplace=True)
		good_dists.to_csv("%s/%s_best_focal_lengths.csv"%(path,path))