############ obtains vexation from low density data by looking at the ratio of independent probabilities for 0 and 1 fly #########


################################
#		input:	-csv file with the counts of flies per bin in an m x m grid at diferent frames of a video
#				the i, j entry of this matrix corresponds to the number of flies jth bin in the ith frame of the video
#				-value of m
#
#
#
#		output:  1 graphs
#
#						-heatmap of interpolated vexaton from the fit
#
#				 1 files
#						-csv with vexation
#
#
#				  prints
#						-Vexation
#
#
#################################


##### measure the runtime of the program
import time
start = time.time()

#### three arguments, the name of the file of the the counts of flies per bin in an m x m grid at diferent frames of a video, the value of m and the lenght of a side of a bin in mm
import sys
if len(sys.argv) != 3:
	print "this program needs (2) arguments to run"
	sys.exit()


#### importing libraries
import numpy as np
from scipy import special
import matplotlib.pyplot as plt
from scipy import interpolate



############## Loading data into a numpy array (ij entry corresponds to the count of flies ith bin at the jth frame of the video)##########################
import csv
datafile=sys.argv[1]
Nmat=list(np.array(list(csv.reader(open(datafile))), dtype=np.int).T)


MaxFlies=np.max(np.array(Nmat).flatten())   #maximum number of flies
Nbins=int(sys.argv[2])*int(sys.argv[2])


#calculating the matrix of histograms (i j entry corresponds to the relative frequency of the occurrence of j flies in the ith bin in the video)
Sis=[]
Vex=[]
for i in range(Nbins):
	hist, bin_edges = np.histogram(Nmat[i] ,bins=np.arange(0, MaxFlies+1 , 1), density=True)
	hist=hist+1e-17    #adding machine precision value to prevemnt infs later
	Sis.append(hist)
	Vex.append(np.log(hist[1]/hist[2])) ##taking the vexation as the log of the ratio of por()0 and prob(1)
####################################################################################################################
Vex=list(reversed(Vex))

#print a file with the vexation of the given dataset
file_vex = open('vex_'+datafile, 'w')
for i in range(np.size(Vex)-1):
	s=str(Vex[i])
	file_vex.write(s+',')
s=str(Vex[np.size(Vex)-1])
file_vex.write(s)
file_vex.close()

#ploting the low density vexation


#plot vexation interpolated as a heat map
x = np.arange(0, np.sqrt(Nbins), 1)
y = np.arange(0, np.sqrt(Nbins), 1)
f = interpolate.interp2d(x, y, Vex, kind='cubic')
xnew = np.arange(0, np.sqrt(Nbins)-1, 1e-2)
ynew = np.arange(0, np.sqrt(Nbins)-1, 1e-2)
znew = f(xnew, ynew)
fig = plt.figure(1, figsize=(10.5,8.5))
plt.imshow(znew)
plt.show()


print Vex
end = time.time()

print "time elapsed in this run " +str(end - start)
