

############ predicts vexation as a function of position from data consisting of counts per bin for all bins at different times #########

################################
#		input:	-csv file with the counts of flies per bin in an m x m grid at diferent frames of a video
#				the i, j entry of this matrix corresponds to the number of flies jth bin in the ith frame of the video
#			    -Area of the bins
#
#				-average area of a fly from the given dataset calculated in least_squares_fit.py (2.812 for 2.csv)
#
#		output:	-csv file with predicted vexation in each bin
#			
#################################

##### measure the runtime of the program
import time
start = time.time()

#### outputs error if the name of the file with the data, the chemical potential are not given as arguments and the area of a fly
import sys
if len(sys.argv) != 3:
	print "this program needs (2) arguments to run"
	sys.exit()


#### importing libraries
import numpy as np
from scipy.optimize import brentq
from scipy import special



#### Loading data into a numpy array (ij entry corresponds to the count of flies ith bin at the jth frame of the video)
import csv
datafile=sys.argv[1]
Nmat=np.array(list(csv.reader(open(datafile))), dtype=np.int).T



##### introducing parameters of the system

alpha=float(sys.argv[3])/float(sys.argv[2])**2



#### Calculating the average number of flies per bin
Average=np.mean(Nmat,axis=1)



print  -alpha*Average/(1-alpha*Average)+np.log(1-alpha*Average)-special.digamma(1+Average)

Vex=-alpha*Average/(1-alpha*Average)+np.log(1-alpha*Average)-special.digamma(1+Average)


#print a file with the vexation of the given dataset
file_vex = open('vex_Dft_'+datafile, 'w')
for i in range(np.size(Vex)-1):
	s=str(Vex[i])
	file_vex.write(s+',')
s=str(Vex[np.size(Vex)-1])
file_vex.write(s)
file_vex.close()




end = time.time()

print "time elapsed in this run " +str(end - start)

