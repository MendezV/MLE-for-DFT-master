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
ZS=np.array(ZS)
Nmat=np.zeros(Nbins)
#print ZS


rings=ZS.size

Nmat[0]=data[0]

for i in range(rings-1):
	a=np.sum(data[ZS[i]:ZS[i]+ZS[i+1]])/(ZS[i+1])
	#print np.sum(ZS[:i+1]),np.sum(ZS[:i+2])
	for j in range(np.sum(ZS[:i+1]),np.sum(ZS[:i+2])):
		Nmat[j]=a



file_Av = open('sim'+datafile, 'w')
for i in range(np.size(Nmat)-1):
	s=str(Nmat[i])
	file_Av.write(s+',')
s=str(Nmat[np.size(Nmat)-1])
file_Av.write(s)
file_Av.close()
