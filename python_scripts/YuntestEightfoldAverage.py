##### measure the runtime of the program
import time
start = time.time()



#### outputs error if the name of the file with the data for vexations is given as argument
import sys
if len(sys.argv) != 3:
	print "this program needs (2) arguments to run"
	sys.exit()



#importing libraries
import numpy as np
from scipy.optimize import brentq
from scipy import special
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


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
Nbins=np.shape(Nmat)[1]
side=np.sqrt(Nbins)
Nmatsim=list(np.zeros((Tframes,Nbins)))


########## error load ##########
datafileErr=sys.argv[2]
ERR= list(csv.reader(open(datafileErr)))[0]
Ebar=np.zeros(len(ERR))
for i in range(len(ERR)):
	Ebar[i]=float(ERR[i])*float(ERR[i])

NERROR=np.zeros(Nbins)
#########################
Variable=0.0
Variable2=0.0
testside=5
Nbins=testside*testside
for r in range(Nbins):
	test=np.zeros(Nbins)
	test[r]=1
	simtest=np.zeros(Nbins)
	print np.reshape(test,(testside,testside))
	
	
	#########

################# average transpose

	for i in range(Nbins):
		#indices that are going to be averaged in for this particular bin and assign to both indices
		i0=i%testside
		j0=(i-i0)/testside
		Val=(test[i0+j0*testside]+test[j0+i0*testside])/2.0
		simtest[i0+j0*testside]=Val
		simtest[j0+i0*testside]=Val

###################


################# flip lr
	auxArr=[]

	for i in range(Nbins):
		auxArr.append(simtest[i])

	for i in range(Nbins):
		#indices that are going to be averaged in for this particular bin
		i0=i%testside
		j0=(i-i0)/testside
		i2=testside-i0-1
		j2=j0
		simtest[i]=auxArr[int(i2+j2*testside)]

###################

################# average flip lr with rot 90 and assign to both indices
	for i in range(Nbins):
		#indices that are going to be averaged in for this particular bin
		i0=i%testside
		j0=(i-i0)/testside
		Val=(simtest[i0+j0*testside]+simtest[j0+i0*testside])/2.0
		simtest[i0+j0*testside]=Val
		simtest[j0+i0*testside]=Val

###################

################# flip ud
	auxArr=[]

	for i in range(Nbins):
		auxArr.append(simtest[i])

	for i in range(Nbins):
		#indices that are going to be averaged in for this particular bin
		i0=i%testside
		j0=(i-i0)/testside
		i3=i0
		j3=testside-j0-1
		simtest[i]=auxArr[int(i3+j3*testside)]

###################


################# average flip ud with array before flip ud and assign just to

	for i in range(Nbins):
		#indices that are going to be averaged in for this particular bin
		i0=i%testside
		j0=(i-i0)/testside
		i3=i0
		j3=testside-j0-1
		Val=(simtest[i0+j0*testside]+simtest[i3+j3*testside])/2.0
		simtest[i3+j3*testside]=Val
	##########
	
	print np.reshape(simtest,(testside,testside))





end = time.time()
print "time elapsed in this run " +str(end - start)
