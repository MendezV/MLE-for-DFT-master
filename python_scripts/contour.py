import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D





import time
start = time.time()

#### three arguments, the name of the file of the the counts of flies per bin in an m x m grid at diferent frames of a video, the value of m and the lenght of a side of a bin in mm
import sys
if len(sys.argv) != 2:
	print "this program needs (1) arguments to run"
	sys.exit()

import csv
datafile=sys.argv[1]
data=np.array(list(csv.reader(open(datafile))), dtype=np.float).T
data=data.flatten()


#-- Generate Data -----------------------------------------

Nbins=data.size

r1=np.zeros(Nbins)
theta1=np.zeros(Nbins)
ZS=[]

for i in range(Nbins):
	ZS.append(2*i+1)
	if(Nbins==np.sum(np.array(ZS))):
		break
#ZS=np.array(ZS)
ZS=2*np.array(ZS)    #for Yunus bins


def determine_ring(i,ZS):
	x=0
	counter=0.0
	r=0.0
	for j in ZS:
		x=i-np.sum(ZS[:(list(ZS).index(j)+1)])
		if(x<=0):
			counter=list(ZS).index(j)
			break
	return counter
def determine_angle(r,i,ZS):
	Latbin=np.sum(ZS[:(r+1)])-i
	return (360*Latbin/np.float(ZS[r]))%360

for i in range(Nbins):
	r1[i]=determine_ring(i+1,ZS)
	theta1[i]=determine_angle(r1[i],i+1,ZS)



X=3.025*r1*np.cos(np.radians(theta1))/6.0
Y=3.025*r1*np.sin(np.radians(theta1))/6.0



fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
c = np.abs(data)
ax.scatter(X, Y, data,c=c ,cmap=cm.YlOrRd)

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

m = cm.ScalarMappable(cmap=cm.YlOrRd)
m.set_array(data)
plt.colorbar(m)

plt.show()

print data





