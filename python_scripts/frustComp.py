############  compares two vexations to see if they differ by scaling and a constant shift#########


################################
#		input:	-Csv file with Vexation calculated from least squares fit in least_squares_fit.py
#				-Csv file with another Vexation calculated from low density data
#
#
#		output:	-scatter plot of the two vexations with linear regression
#				-parameters of the linear regression
#				-linear plot fitting the scatter plot
#
#
#################################

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
from scipy import stats
import matplotlib.pyplot as plt
from scipy import interpolate


#importing frustrations
datos1 = np.loadtxt(sys.argv[1])
F=datos1[:,1]
print np.shape(datos1[:,1])

#importing frustrations
datos2 = np.loadtxt(sys.argv[2])
F2=datos2[:,1]
print np.shape(datos2[:,1])

Nmax=np.min([np.max(datos1[:,0]),np.max(datos2[:,0])])
print Nmax,datos1[:Nmax,0],datos1[:Nmax,1]



fig, ax = plt.subplots()
SCF1=ax.scatter(datos1[:,0], F[:])
SCF2=ax.scatter(datos2[:,0], F2[:],c='g')
Diff=ax.scatter(datos1[:Nmax,0], F[:Nmax]-F2[:Nmax],c='r', marker='x')
ax.legend((SCF1, SCF2,Diff), ('Circle', 'Circle Yun','C-CY'))
plt.show()

#linear plot to the plot of vexation1 vs vexation2
slope, intercept, r_value, p_value, std_err = stats.linregress(datos1[:4,0], F[:4]-F2[:4])

fig, ax = plt.subplots()
SCF1=ax.scatter(datos1[:,0], F[:])
SCF2=ax.scatter(datos2[:,0], F2[:],c='g')
Diff=ax.scatter(datos1[:Nmax,0], F[:Nmax]-F2[:Nmax],c='r', marker='x')
ax.plot(datos1[:Nmax,0], slope*datos1[:Nmax,0]+intercept)
ax.legend((SCF1, SCF2,Diff), ('Circle', 'Square','C-S'))
plt.show()

plt.scatter(datos1[:Nmax,0], F[:Nmax]-F2[:Nmax])
plt.plot(datos1[:Nmax,0], slope*datos1[:Nmax,0]+intercept)
plt.show()

print " for this two vexations, the relative factor is ", slope
print "\n the r value of the fit is ", r_value
print "\n the standard error of the fit is ", std_err
print "\n the p value of the fit is ", p_value


end = time.time()

print "time elapsed in this run " +str(end - start)
