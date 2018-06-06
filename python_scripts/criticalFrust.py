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
from scipy import stats
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


#importing frustrations

datos = np.loadtxt(sys.argv[1])

F = interp1d(datos[:,0], datos[:,1], kind='cubic')
slope, intercept, r_value, p_value, std_err = stats.linregress(np.log(datos[950:,0]), np.log(datos[950:,1]))

plt.scatter(np.log(datos[950:,0]), np.log(datos[950:,1]))
plt.plot(np.log(datos[950:,0]), slope*np.log(datos[950:,0])+intercept)
plt.show()

print slope , intercept
end = time.time()

print "time elapsed in this run " +str(end - start)


