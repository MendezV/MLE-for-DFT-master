

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
from scipy import special
import matplotlib.pyplot as plt



#### Loading data into a numpy array (ij entry corresponds to the count of flies ith bin at the jth frame of the video)
import csv
datafile=sys.argv[1]
Nmat=np.array(list(csv.reader(open(datafile))), dtype=np.int).T

Av=[]
MaxFlies=np.max(np.array(Nmat).flatten())
ns=np.arange(0,MaxFlies+1,1)
Nbins=np.shape(Nmat)[0]
Tframes=np.shape(Nmat)[1]
print MaxFlies, ns


#### Calculating the average number of flies per bin
for i in range(Nbins):
	hist, bin_edges = np.histogram(Nmat[i] ,bins=np.arange(0,MaxFlies+2,1), density=True)
	Av.append(np.sum(ns*hist))

Average=np.array(Av)



print Average


#print a file with the vexation of the given dataset
file_Av = open('Average_N'+datafile, 'w')
for i in range(np.size(Average)-1):
	s=str(Average[i])
	file_Av.write(s+',')
s=str(Average[np.size(Average)-1])
file_Av.write(s)
file_Av.close()

################ Error bars ##################
NERROR=np.zeros(Average.size)

def sigmaDiri(x,Tframes):
	return np.sqrt(Tframes*x*(Tframes-Tframes*x)/(Tframes*Tframes*(Tframes+1)))

for i in range(Nbins):
	hist, bin_edges = np.histogram(Nmat[i] ,bins=np.arange(0,MaxFlies+2,1), density=True)
	NERROR[i]=np.sqrt(np.sum(np.square(ns*sigmaDiri(hist,Tframes))))



print NERROR
plt.xlim(-1,Average.size +1)
plt.errorbar(np.arange(Average.size),Average,yerr=NERROR, fmt='o')

plt.show()





file_vex = open('avErrorb_'+datafile, 'w')
for i in range(np.size(NERROR)-1):
	s=str(NERROR[i])
	file_vex.write(s+',')
s=str(NERROR[np.size(NERROR)-1])
file_vex.write(s)
file_vex.close()


end = time.time()

print "time elapsed in this run " +str(end - start)

