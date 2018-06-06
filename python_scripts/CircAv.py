import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D





import time
start = time.time()


import sys
if len(sys.argv) != 3:
	print "this program needs (2) arguments to run"
	sys.exit()

import csv
datafile=sys.argv[1]
data=np.array(list(csv.reader(open(datafile))), dtype=np.float).T
data=data.flatten()

datafile2=sys.argv[2]
Err=np.array(list(csv.reader(open(datafile2))), dtype=np.float).T
Err=Err.flatten()


#-- Generate Data -----------------------------------------

Nbins=data.size

r1=np.zeros(Nbins)
theta1=np.zeros(Nbins)
ZS=[]
ZS.append(0)

for i in range(Nbins):
	ZS.append(2*i+1)
	if(Nbins==2*np.sum(np.array(ZS))):
		break
ZS=2*np.array(ZS)
Nmat=np.zeros(Nbins)
print ZS, Nbins


rings=ZS.size

Nmat[0]=data[0]
bars=np.zeros(Nbins)
bars[0]=Err[0]

for i in range(rings-1):
	a=np.sum(data[ZS[i]:ZS[i]+ZS[i+1]])/(ZS[i+1])
	b=np.sqrt(np.sum(np.square(Err[ZS[i]:ZS[i]+ZS[i+1]])))/(ZS[i+1])
	print np.sum(ZS[:i+1]),np.sum(ZS[:i+2])
	for j in range(np.sum(ZS[:i+1]),np.sum(ZS[:i+2])):
		Nmat[j]=a
		bars[j]=b



file_Av = open('sim'+datafile, 'w')
for i in range(np.size(Nmat)-1):
	s=str(Nmat[i])
	file_Av.write(s+',')
s=str(Nmat[np.size(Nmat)-1])
file_Av.write(s)
file_Av.close()

file_E = open('avtot_Errb'+datafile, 'w')
for i in range(np.size(bars)-1):
	s=str(bars[i])
	file_E.write(s+',')
s=str(bars[np.size(bars)-1])
file_E.write(s)
file_E.close()

