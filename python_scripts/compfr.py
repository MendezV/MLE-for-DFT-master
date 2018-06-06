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

import seaborn as sns
sns.set(font_scale=2)
sns.set_style("whitegrid")


#### outputs error if the name of the file with the data for vexations is given as argument
import sys
if len(sys.argv) != 7:
	print "this program needs (6) arguments to run"
	sys.exit()



#importing libraries
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from scipy import interpolate


#Loading Vexations into the program
datos = np.loadtxt(sys.argv[1])
#F = interp1d(datos[:,0], datos[:,1], kind='cubic')
delta=datos[1,0]-datos[0,0]
F=datos[:,1]
N=np.arange(np.size(F))
slope, intercept, r_value, p_value, std_err = stats.linregress(N[0:2], F[0:2])
plt.plot(N, F-N*slope-intercept)


datos = np.loadtxt(sys.argv[2])
#F = interp1d(datos[:,0], datos[:,1], kind='cubic')
delta=datos[1,0]-datos[0,0]
F=datos[:,1]
N=np.arange(np.size(F))
slope, intercept, r_value, p_value, std_err = stats.linregress(N[0:2], F[0:2])
plt.plot(N, F-N*slope-intercept)


datos = np.loadtxt(sys.argv[3])
#F = interp1d(datos[:,0], datos[:,1], kind='cubic')
delta=datos[1,0]-datos[0,0]
F=datos[:,1]
N=np.arange(np.size(F))
slope, intercept, r_value, p_value, std_err = stats.linregress(N[0:2], F[0:2])
plt.plot(N, F-N*slope-intercept)


datos = np.loadtxt(sys.argv[4])
#F = interp1d(datos[:,0], datos[:,1], kind='cubic')
delta=datos[1,0]-datos[0,0]
F=datos[:,1]
N=np.arange(np.size(F))
slope, intercept, r_value, p_value, std_err = stats.linregress(N[0:2], F[0:2])
plt.plot(N, F-N*slope-intercept)

datos = np.loadtxt(sys.argv[5])
#F = interp1d(datos[:,0], datos[:,1], kind='cubic')
delta=datos[1,0]-datos[0,0]
F=datos[:,1]
N=np.arange(np.size(F))
slope, intercept, r_value, p_value, std_err = stats.linregress(N[0:2], F[0:2])
plt.plot(N, F-N*slope-intercept)

datos = np.loadtxt(sys.argv[6])
#F = interp1d(datos[:,0], datos[:,1], kind='cubic')
delta=datos[1,0]-datos[0,0]
F=datos[:,1]
N=np.arange(np.size(F))
slope, intercept, r_value, p_value, std_err = stats.linregress(N[0:2], F[0:2])
plt.plot(N, F-N*slope-intercept)

plt.show()


end = time.time()
print "time elapsed in this run " +str(end - start)
