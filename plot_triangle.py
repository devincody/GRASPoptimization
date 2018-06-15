import numpy as np                 
import os                          
import matplotlib.pyplot as plt    
from collections import defaultdict
import pandas as pd
from scipy.interpolate import griddata


metadata_path = "../metadata"
os.chdir(metadata_path)
cwd = os.getcwd()      
items = os.listdir(cwd)

partition = 20
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
					print("dropping last column", keys[i])
			except:
				pass
		keys = a.keys()
		fig = plt.figure(figsize=(10,10))

		plots = len(keys)-3

		for i in range (1,plots+1):
			for j in range(i,plots+1):

				maxi= []
				maxi.append(a.iloc[0].values)
				for c in range(len(a[keys[0]])):
					diff = np.abs(a.iloc[c].values - maxi[-1])
					if (diff[i]+diff[j] < 1E-2):
						if (a.iloc[c].values[-1] < maxi[-1][-1]):
							maxi[-1] = a.iloc[c].values
					else:
						maxi.append(a.iloc[c].values)


				maxi = np.array(maxi)
				locations = maxi[:,(i,j)]
				values = maxi[:,-1]


				maxx = max(locations[:,0])
				minx = min(locations[:,0])
				miny = min(locations[:,1])
				maxy = max(locations[:,1])

				plt.subplot(plots,plots,i+plots*(j-1))
				if (i == j):

					plt.scatter(locations[:,0], -values, marker='.', alpha=.7)

				else:
					
					grid_x, grid_y = np.mgrid[minx:maxx:partition*1j, miny:maxy:partition*1j]
					try:
						method = 'linear'
						grid_z0 = griddata(locations, values, (grid_x, grid_y), method=method)
					except:
						method = 'nearest'
						grid_z0 = griddata(locations, values, (grid_x, grid_y), method=method)

					plt.imshow(-grid_z0.T, extent=(minx,maxx, miny,maxy),
								 origin='lower', cmap = 'plasma', 
								 interpolation = 'quadric', aspect='auto')

				if (i == 1):
					plt.ylabel(keys[j])
				
				if (j == plots):
					plt.xlabel(keys[i])
				


		plt.subplots_adjust(left=0.1, bottom=0.05, right=.95, top=.95, wspace=0.3, hspace=0.3)
		plt.savefig(path + '/' + "triangle" + '.png')
		plt.close()
		# plt.show()



# plt.subplot(1,2,2)
# plt.scatter(locations[:,0], locations[:,1], c= -values, cmap = 'plasma')
# plt.colorbar()
# plt.xlabel(keys[i])
# plt.ylabel(keys[j])
# plt.title("scatter")
# 











# maxx = max(a['x'])            
# minx = min(a['x'])            
# miny = min(a['y'])            
# maxy = max(a['y'])            
# X = np.linspace(minx, maxx, partition+1)
# Y = np.linspace(miny, maxy, partition+1)
# A = np.zeros((partition,partition))          





# plt.imshow(A, cmap='hot', interpolation='bicubic')  
# plt.colorbar()
# plt.xlabel("X")
# plt.ylabel("Y")
# plt.xticks(np.linspace(0,partition,np.min((partition/2, 12))),["%4.2f"%x for x in np.linspace(minx,maxx,np.min((partition/2, 12)))] )
# plt.yticks(np.linspace(0,partition,np.min((partition/2, 12))),["%4.2f"%x for x in np.linspace(miny,maxy,np.min((partition/2, 12)))] )
# plt.show()