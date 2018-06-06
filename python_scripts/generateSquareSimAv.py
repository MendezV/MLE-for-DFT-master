##### measure the runtime of the program
import time
start = time.time()



#### outputs error if the name of the file with the data for vexations is given as argument
import sys
if len(sys.argv) != 2:
	print "this program needs (1) arguments to run"
	sys.exit()



#importing libraries
import numpy as np
from scipy.optimize import brentq
from scipy import special
from scipy.interpolate import interp1d


#Loading Vexations into the program
import csv
datafile=sys.argv[1]
vex= list(csv.reader(open(datafile)))[0]
Vexation=np.zeros(len(vex))
for i in range(len(vex)):
	Vexation[i]=float(vex[i])

import csv
datafile=sys.argv[1]
Nmat1=np.array(list(csv.reader(open(datafile), delimiter=',')), dtype=np.float)

#print Nmat1
Nmat=list(Nmat1)
Tframes=np.shape(Nmat)[0]
Nbins=np.shape(Nmat)[1]
side=np.sqrt(Nbins)
test=np.arange(Nbins)
test2=np.zeros(Nbins)

Nmatsim=list(np.zeros((Tframes,Nbins)))
for i in range(Nbins):
	#indices that are going to be averaged in for this particular bin
	i0=i%side
	j0=(i-i0)/side
	i1=side-i0-1
	j1=side-i0-1
	i2=side-i0-1
	j2=j0
	i3=i0
	j3=side-j0-1
	#print i0,j0,i1,j1,i2,j2,i3,j3
	test2[i]=test[i]+test[i1+j1*side]+test[i2+j2*side]+test[i3+j3*side]
	test2[i]/=4.0


for j in range(Tframes):
	for i in range(Nbins):
	#indices that are going to be averaged in for this particular bin
		i0=i%side
		j0=(i-i0)/side
		i1=side-i0-1
		j1=side-j0-1
		i2=side-i0-1
		j2=j0
		i3=i0
		j3=side-j0-1
		Nmatsim[j][i]=Nmat[j][i]+Nmat[j][i1+j1*side]+Nmat[j][i2+j2*side]+Nmat[j][i3+j3*side]
		Nmatsim[j][i]/=4.0


#print Nmatsim

### printin everything
f = open('sim'+datafile, 'w')
for j in range(Tframes):
	
	for i in range(Nbins-1):
		s = str(float(Nmatsim[j][i]))+","
		f.write(s)
	s = str(float(Nmatsim[j][Nbins-1]))
	f.write(s)
	f.write("\n")
f.close()



end = time.time()
print "time elapsed in this run " +str(end - start)