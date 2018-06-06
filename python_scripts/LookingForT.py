##### measure the runtime of the program
import time
start = time.time()

#### outputs error if the name of the file with the data, the chemical potential are not given as arguments and the area of a fly
import sys
if len(sys.argv) != 6:
	print "this program needs (5) arguments to run"
	sys.exit()


#### importing libraries
import numpy as np
from scipy.optimize import brentq
from scipy import special
import matplotlib.pyplot as plt
from scipy import stats


#### Loading data into a numpy array (ij entry corresponds to the count of flies ith bin at the jth frame of the video)
import csv
datafile1=sys.argv[1]
Nmat1=np.array(list(csv.reader(open(datafile1))), dtype=np.int).T



#Loading Vexations 1 into the program
datafileV1=sys.argv[2]
vex1= list(csv.reader(open(datafileV1)))[0]
Vexation1=np.zeros(len(vex1))
for i in range(len(vex1)):
	Vexation1[i]=float(vex1[i])


datafile2=sys.argv[3]
Nmat2=np.array(list(csv.reader(open(datafile2))), dtype=np.int).T



#Loading Vexations 2 into the program

datafileV2=sys.argv[4]
vex2= list(csv.reader(open(datafileV2)))[0]
Vexation2=np.zeros(len(vex2))
for i in range(len(vex2)):
	Vexation2[i]=float(vex2[i])



N_fluc1=np.power(np.std(Nmat1, axis=1),2)
Average1=np.mean(Nmat1,axis=1)



N_fluc2=np.power(np.std(Nmat2, axis=1),2)
Average2=np.mean(Nmat2,axis=1)

alpha=float(sys.argv[5])/(4*4)
def mu_der(n):
	return 2*alpha/(1-alpha*n)**2+(alpha**2)*n/((1-alpha*n)**2)

inversebeta1=np.zeros(np.size(Average1))
for i in range(np.size(Average1)):
	inversebeta1[i]=N_fluc1[i]*mu_der(Average1[i])/(1-N_fluc1[i]*special.polygamma(1, 1+Average1[i]))
	print alpha*Average1[i]/(1-alpha*Average1[i])-np.log(1-alpha*Average1[i])+special.digamma(1+Average1[i])+Vexation1[i]

inversebeta2=np.zeros(np.size(Average2))
for i in range(np.size(Average2)):
	inversebeta2[i]=N_fluc2[i]*mu_der(Average2[i])/(1-N_fluc2[i]*special.polygamma(1, 1+Average2[i]))
	print alpha*Average2[i]/(1-alpha*Average2[i])-np.log(1-alpha*Average2[i])+special.digamma(1+Average2[i])+Vexation2[i]

print inversebeta1, np.mean(inversebeta1)
print inversebeta2,np.mean(inversebeta2)

slope, intercept, r_value, p_value, std_err = stats.linregress(inversebeta1*Vexation1, inversebeta2*Vexation2)

plt.scatter(inversebeta1*Vexation1,inversebeta2*Vexation2)
plt.plot(inversebeta1*Vexation1,slope*inversebeta1*Vexation1+intercept)

plt.show()

print " for this two vexations, the relative factor is ", slope
print "\n the r value of the fit is ", r_value
print "\n the standard error of the fit is ", std_err
print "\n the p value of the fit is ", p_value






end = time.time()

print "time elapsed in this run " +str(end - start)

