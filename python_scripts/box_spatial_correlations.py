

############ Calculates the spatial correlations between boxes averaged over time #########

################################
#		input:	-csv file with the counts of flies per bin in an m x m grid at diferent frames of a video
#				the i, j entry of this matrix corresponds to the number of flies jth bin in the ith frame of the video
#				-value of m
#
#
#		output:  3 graphs
#						-heatmap of the spatial correlation matrix <(n(xi)-<n(xi)>)(n(xj)-<n(xj)>)>_{t} (where i and j are diferent bins in the m x m grid)
#						-heatmap of the spatial correlation as a function of distance, averaging over all correlation functions between bins taking into account orientation and dividing by the relevat number of correlation functions that can be averaged (e.g there is just one correlation function which has a neighbour m units to the right and m units down, namely the correlation function of box 0 with the rest)
#						-histogram of the spatial correlation as a function of distance calculated the same way as described before
#################################

##### measure the runtime of the program
import time
start = time.time()


###### first argument is the data with the counts of flies per bin in an m x m grid , the second argument is m
import sys
if len(sys.argv) != 3:
	print "this program needs (2) arguments to run"
	sys.exit()

### importing libraries
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


#### Loading data into a python list (ij entry corresponds to the count of flies jth bin in the ith frame of the video)
import csv
datafile=sys.argv[1]
prob = list(csv.reader(open(datafile)))
prob2=[[0 for i in range(np.shape(prob)[0])] for j in range(np.shape(prob)[1])]


#transposing list
for i in range(np.shape(prob)[0]):
	for j in range(np.shape(prob)[1]):
		prob2[j][i]=float(prob[i][j])

#substracting the mean to each bin
for i in range(np.shape(prob2)[0]):      ##do not substract mean if only one time data set, eg vex_dif
	prob2[i]-=np.mean(prob2[i])

#calculating spatial correlation as the average of the diference <(n(x)-<n(x)>)(n(x')-<n(x')>)>
SPCorr=np.dot(np.array(prob2,dtype=np.float),np.array(prob2,dtype=np.float).T)/float(np.size(prob2[0]))

#ploting the matrix in wich the entry ij corresponds to the time average of the product (n(xi)-<n(xi)>)(n(xj)-<n(xj)>)
data = SPCorr
fig, ax = plt.subplots()
heatmap = ax.pcolor(data,cmap=plt.cm.Blues)
plt.show()

#reshaping one of the lines in the matrix taking into account the shape of the grid
m=int(sys.argv[2])
onebox=np.reshape(data.T[0],(m,m),order='C')


#defining the matrix that will display the correlations as a function of distance, wich will be the weighted average of the correlation cuntion of each bin with the others
len_newbox=2*np.shape(onebox)[0]-1
Spat=np.zeros((len_newbox,len_newbox),dtype=np.float)
counter=np.zeros((len_newbox,len_newbox),dtype=np.float)

#averaging the correlation function of each box (6x6) with the others by shifting them to get center of the matrix Spat as the average of each bin with itself
for i in range(np.shape(data)[0]):
	wall1 = np.zeros((len_newbox,len_newbox),dtype=np.float)
	wall2 = np.zeros((len_newbox,len_newbox),dtype=np.float)
	block = np.reshape(data.T[i],(np.shape(onebox)[0],np.shape(onebox)[0]),order='C')
	x = m-1-(i-(i%(np.shape(block)[0])))/(np.shape(block)[0])
	y =-(i%(np.shape(block)[0]))+np.shape(block)[1]-1
	wall1[x:x+block.shape[0], y:y+block.shape[1]] = block
	wall2[x:x+block.shape[0], y:y+block.shape[1]] = np.ones((np.shape(onebox)[0],np.shape(onebox)[0]),dtype=np.float)
	Spat+=wall1
	counter+=wall2

### normalizing by the amount of 6x6 for the data that was analized for the presentation
Spat/=counter
fig, ax = plt.subplots()
heatmap = ax.pcolor(Spat,cmap=plt.cm.Blues)
plt.show()


####plotting the 3d histogram of the spatial correlations

thedata=Spat
x = thedata[:,0]    # data from the first column
y = thedata[:,1]    # data from the second column
hist, xedges, yedges = np.histogram2d(x, y, bins=np.arange(0, len_newbox+1 , 1))

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
