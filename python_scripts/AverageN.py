

############ predicts vexation as a function of position from data consisting of counts per bin for all bins at different times #########

################################
#		input:	-csv file with the counts of flies per bin in an m x m grid at diferent frames of a video
#				the i, j entry of this matrix corresponds to the number of flies jth bin in the ith frame of the video
#
#
#		output:	-csv file with average number of flies n each bin
#
#################################

##### measure the runtime of the program
import time
start = time.time()

#### outputs error if the name of the file with the data, the chemical potential are not given as arguments and the area of a fly
import sys
if len(sys.argv) != 2:
	print "this program needs (1) arguments to run"
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



#### Calculating the average number of flies per bin
Average=np.mean(Nmat,axis=1)



print Average


#print a file with the vexation of the given dataset
file_Av = open('Average_N_'+datafile, 'w')
for i in range(np.size(Average)-1):
	s=str(Average[i])
	file_Av.write(s+',')
s=str(Average[np.size(Average)-1])
file_Av.write(s)
file_Av.close()




end = time.time()

print "time elapsed in this run " +str(end - start)

