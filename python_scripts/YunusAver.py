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


Nmat=list(Nmat1)
Tframes=np.shape(Nmat)[0]
Nbins=9
side=np.sqrt(Nbins)
test0=np.arange(Nbins)
test=np.arange(Nbins)
test2=np.zeros(Nbins)
test3=np.zeros(Nbins)
test4=np.zeros(Nbins)
test5=np.zeros(Nbins)
test6=np.zeros(Nbins)
test7=np.zeros(Nbins)
test8=np.zeros(Nbins)


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
	print i0,j0,i1,j1,i2,j2,i3,j3
	test[i]=test0[i0+j0*side]
	test2[i]=test0[i1+j1*side]
	test3[i]=test0[i2+j2*side]
	test4[i]=test0[i3+j3*side]
	test5[i]=test0[j0+i0*side]
	test6[i]=test0[j1+i1*side]
	test7[i]=test0[j2+i2*side]
	test8[i]=test0[j3+i3*side]
	


print test.reshape([3,3]),"identity"
print test2.reshape([3,3]), "rot 180"
print test3.reshape([3,3]), "fliplr"
print test4.reshape([3,3]), "flipud"
print test5.reshape([3,3]), "transpose"
print test6.reshape([3,3]), "transpose non-prinpcipal diagonal"
print test7.reshape([3,3]), "rot 90"
print test8.reshape([3,3]), "rot 270"


Nmat=list(Nmat1)
Tframes=np.shape(Nmat)[0]
Nbins=np.shape(Nmat)[1]
side=np.sqrt(Nbins)
Nmatsim=list(np.zeros((Tframes,Nbins)))

Val=0


################# average transpose
for j in range(Tframes):
	for i in range(Nbins):
	#indices that are going to be averaged in for this particular bin and assign to both indices
		i0=i%side
		j0=(i-i0)/side
		Val=(Nmat[j][i0+j0*side]+Nmat[j][j0+i0*side])/2.0
		Nmatsim[j][i0+j0*side]=Val
		Nmatsim[j][j0+i0*side]=Val

###################


################# flip lr
auxArr=[]
for j in range(Tframes):
	for i in range(Nbins):
		auxArr.append(Nmatsim[j][i])

for j in range(Tframes):
	for i in range(Nbins):
		#indices that are going to be averaged in for this particular bin
		i0=i%side
		j0=(i-i0)/side
		i2=side-i0-1
		j2=j0
		Nmatsim[j][i]=auxArr[int(i2+j2*side)]

###################

################# average flip lr with rot 90 and assign to both indices
for j in range(Tframes):
	for i in range(Nbins):
		#indices that are going to be averaged in for this particular bin
		i0=i%side
		j0=(i-i0)/side
		Val=(Nmatsim[j][i0+j0*side]+Nmatsim[j][j0+i0*side])/2.0
		Nmatsim[j][i0+j0*side]=Val
		Nmatsim[j][j0+i0*side]=Val

###################

################# flip ud
auxArr=[]
for j in range(Tframes):
	for i in range(Nbins):
		auxArr.append(Nmatsim[j][i])

for j in range(Tframes):
	for i in range(Nbins):
		#indices that are going to be averaged in for this particular bin
		i0=i%side
		j0=(i-i0)/side
		i3=i0
		j3=side-j0-1
		Nmatsim[j][i]=auxArr[int(i3+j3*side)]

###################


################# average flip ud with array before flip ud and assign just to
for j in range(Tframes):
	for i in range(Nbins):
		#indices that are going to be averaged in for this particular bin
		i0=i%side
		j0=(i-i0)/side
		i3=i0
		j3=side-j0-1
		Val=(Nmatsim[j][i0+j0*side]+Nmatsim[j][i3+j3*side])/2.0
		Nmatsim[j][i3+j3*side]=Val

###################
# overwritting the aray while averaging
# averaging just three times , too little symmetries
# averaging with pairs of bins, there may be even 4 or 8 equivalent ones
# some symmetries are not being averaged, the array is just being modified

# things to try, plot heatmap with values on top to double check
# do not handpick the bins that should be symmetric









### printin everything
f = open('simYun'+datafile, 'w')
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
