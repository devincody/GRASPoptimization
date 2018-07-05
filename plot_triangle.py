import numpy as np                 
import os                   
       
import matplotlib.pyplot as plt    
from collections import defaultdict
import pandas as pd
from scipy.interpolate import griddata
from scipy.signal import convolve2d

import scipy.stats as st

def gkern(kernlen=5, nsig=3):
    """Returns a 2D Gaussian kernel array."""
    interval = (2*nsig+1.)/(kernlen)
    x = np.linspace(-nsig-interval/2., nsig+interval/2., kernlen+1)
    kern1d = np.diff(st.norm.cdf(x))
    kernel_raw = np.sqrt(np.outer(kern1d, kern1d))
    kernel = kernel_raw/kernel_raw.sum()
    return kernel

sz_kernel = 7
kernel_sharpness = .3
mid = int((sz_kernel-1)/2)

# kernel = np.ones((sz_kernel,sz_kernel))*.8
# kernel[1:-1,1:-1] = np.ones((5,5))*.9
# kernel[2:-2,2:-2] = np.ones((3,3))*.95
# kernel[mid,mid] = 1
kernel = np.array(gkern(sz_kernel,kernel_sharpness))
kernel /= (kernel[mid,mid])
# kernel[mid,mid] = 1

print ("ker {} sum = {}".format(kernel, np.sum(kernel)))


metadata_path = "../metadata"
os.chdir(metadata_path)
cwd = os.getcwd()      
items = os.listdir(cwd)

partition = 40
for it in items:
	if len(it) > 4 and ".csv" in it[-5:]:
		print("Processing file ", it)
		path = it[:-4]
		if not os.path.exists(path):
			os.mkdir(path)

		a = pd.read_csv(it, skiprows = 2, header = 0, index_col=False)
		keys = a.keys()
		for i in range(len(keys)):
			try:
				if "Unnamed" in keys[i] or "alpha" in keys[i] or "100ohm" in keys[i]:
					a = a.drop(keys[i], axis = 1)
					print("dropping column", keys[i])
			except:
				pass
		keys = a.keys()

		try:
			loss = ' Loss'
			a[loss]
		except:
			loss = 'Loss'
			a[loss]

		fig = plt.figure(figsize=(15,15))

		plots = len(keys)-3

		for i in range (1,plots+1):
			for j in range(i,plots+1):



				maxi= []
				# a.sort_values(by=[keys[i], keys[j]], inplace = True)
				print("keys: {} and {}".format(keys[i], keys[j]))

				maxx = max(a[keys[i]])
				minx = min(a[keys[i]])
				miny = min(a[keys[j]])
				maxy = max(a[keys[j]])

				xbin = np.linspace(minx, maxx, partition)
				ybin = np.linspace(miny, maxy, partition)
				x_indicies = np.digitize(a[keys[i]], xbin, True)
				y_indicies = np.digitize(a[keys[j]], ybin, True)

				array = np.zeros((partition,partition))

				for c in range(len(x_indicies)):
					# print("c {}".format(c))
					array[x_indicies[c], y_indicies[c]] = min(a.iloc[c].values[-1], array[x_indicies[c], y_indicies[c]])

				# np.savetxt("a.npy", array)
				# print("saved")



				# array = np.minimum(convolve2d(array, kernel, mode = "same"), array)
				big = np.zeros((partition+sz_kernel-1, partition+sz_kernel-1))
				big[mid:-mid,mid:-mid] = array

				out_array = np.zeros((partition, partition))

				for xx in range(partition):
					for yy in range(partition):
						out_array[xx,yy] = np.min(big[xx:xx+sz_kernel, yy:yy+sz_kernel]*kernel)

				array = out_array
				# print(array)

				array[array == 0] = 'nan'

				plt.subplot(plots,plots,i+plots*(j-1))
				if (i == j):
					# Plot scatter plot if on diagonal
					# print("lens {} and {}".format(len(a[keys[i]]), len(a[loss])))
					plt.scatter(a[keys[i]], -a[loss], marker='.', alpha=.7)

				else:
					# Else plot interpolated phase space
					# grid_x, grid_y = np.mgrid[minx:maxx:partition*1j, miny:maxy:partition*1j]
					# try:
					# 	method = 'linear'
					# 	grid_z0 = griddata(locations, values, (grid_x, grid_y), method=method)
					# except:
					# 	method = 'nearest'
					# 	grid_z0 = griddata(locations, values, (grid_x, grid_y), method=method)

					plt.imshow(-array.T, extent=(minx,maxx, miny,maxy),
								 origin='lower', cmap = 'plasma', 
								 interpolation = 'quadric', aspect='auto')

				# plot label if on the edge
				if (i == 1):
					plt.ylabel(keys[j])
				
				if (j == plots):
					plt.xlabel(keys[i])
				


		plt.subplots_adjust(left=0.1, bottom=0.05, right=.95, top=.95, wspace=0.35, hspace=0.35)
		plt.savefig(path + '/' + "triangle({},{})".format(sz_kernel, kernel_sharpness) + '.png')
		plt.close()
		plt.show()




