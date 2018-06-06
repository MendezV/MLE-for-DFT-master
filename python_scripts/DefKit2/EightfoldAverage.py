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
for j in range(Tframes):
	ind=[]
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
		ind=[i0+j0*side,i1+j1*side,i2+j2*side,i3+j3*side,j0+i0*side,j1+i1*side,j2+i2*side,j3+i3*side]
		setind=set(ind)
		ind2=list(setind)
		print len(ind2)
		Variable=0.0
		Variable2=0.0
		for k in ind2:
			Variable+=Nmat[j][k]
			Variable2+=Ebar[k]
		for k in ind2:      #### doing this to avoid a machine precision bug
			Nmatsim[j][k]=Variable/float(len(ind2))
			NERROR[k]=np.sqrt(Variable2)/float(len(ind2))



print NERROR,np.array(Nmatsim[0]), len(list(set(np.array(Nmatsim[0])))), len(list(set(NERROR)))
plt.xlim(-1,Nbins +1)
plt.errorbar(np.arange(Nbins),np.array(Nmatsim[0]), yerr=NERROR, fmt='o')
plt.show()



### printin' everything
f = open('sim'+datafile, 'w')
for i in range(Nbins-1):
	s = str(float(Nmatsim[j][i]))+","
	f.write(s)
s = str(float(Nmatsim[j][Nbins-1]))
f.write(s)
f.write("\n")
f.close()


file_vex = open('avtotErrorb_'+datafile, 'w')
for i in range(np.size(NERROR)-1):
	s=str(NERROR[i])
	file_vex.write(s+',')
s=str(NERROR[np.size(NERROR)-1])
file_vex.write(s)
file_vex.close()


end = time.time()
print "time elapsed in this run " +str(end - start)
