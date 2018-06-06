############ plots a heatmap of the data that is inserted as a csv file #########


################################
#		input:	-csv file with the counts of flies per bin in an m x m grid at diferent frames of a video
#				the i, j entry of this matrix corresponds to the number of flies jth bin in the ith frame of the video
#				-value of m
#
#
#		output:  2 graphs
#						-heatmap of the data
#						-3d histogram of the data
#
#################################

##### measure the runtime of the program
import time
start = time.time()

#### three arguments, the name of the file of the the counts of flies per bin in an m x m grid at diferent frames of a video, the value of m and the lenght of a side of a bin in mm
import sys
if len(sys.argv) != 3:
	print "this program needs (2) arguments to run"
	sys.exit()


#### importing libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D

############## Loading data into a numpy array (ij entry corresponds to the count of flies ith bin at the jth frame of the video)##########################
import csv
datafile=sys.argv[1]
data=np.array(list(csv.reader(open(datafile))), dtype=np.float).T
Nbins=int(sys.argv[2])

print np.shape(data)

data2=np.reshape(data,(Nbins,Nbins),order='C')
fig, ax = plt.subplots()
for y in range(data2.shape[0]):
	for x in range(data2.shape[1]):
		plt.text(x + 0.5, y + 0.5, '%.4f' % data2[y, x],
			 horizontalalignment='center',
			 verticalalignment='center',
			 )
heatmap = ax.pcolor(data2,cmap=plt.cm.Blues)
plt.show()


####plotting the 3d histogram of the spatial correlations

thedata=data2
x = thedata[:,0]    # data from the first column
y = thedata[:,1]    # data from the second column
hist, xedges, yedges = np.histogram2d(x, y, bins=np.arange(0, Nbins+1 , 1))

elements = (len(xedges) - 1) * (len(yedges) - 1)    # number of boxes
xpos, ypos = np.meshgrid(xedges[:-1]+0.25, yedges[:-1]+0.25)
xpos = xpos.flatten()           # x-coordinates of the bars
ypos = ypos.flatten()           # y-coordinates of the bars
zpos = np.zeros(elements)       # zero-array
dx = 0.5 * np.ones_like(zpos)   # length of the bars along the x-axis
dy = dx.copy()                  # length of the bars along the y-axis
dz = thedata.flatten()             # height of the bars



fig = plt.figure()
ax = Axes3D(fig)
#ax.set_title('$Time$ $Average$ $of$ $the$ $number$ $of$ $flies$ $per$ $box$ $(data)$')

ax.bar3d(xpos, ypos, zpos,      # lower corner coordinates
		 dx, dy, dz,            # width, depth and height
		 color='r',
		 alpha=0.7          # transparency of the bars
		 )
plt.show()



end = time.time()

print "time elapsed in this run " +str(end - start)
