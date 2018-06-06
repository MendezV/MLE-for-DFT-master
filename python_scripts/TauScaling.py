############ makes a least squares fit to the input data for a Boltzmann distribution consistent with a non linear correction given by volume exclusion #########


################################
#		input:	-csv file with the counts of flies per bin in an m x m grid at diferent frames of a video
#				the i, j entry of this matrix corresponds to the number of flies jth bin in the ith frame of the video
#				-value of m
#				-length of the side of a square bin in mm
#
#
#		output:  1 graph
#						-frustration as a function of the number of flies in a bin
#						-frustration as a function of the density of flies in mm^2
#						-heatmap of interpolated vexaton from the fit
#						-time correlations of the system
#				 1 files
#						-csv with vexation from fit
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
#import matplotlib.pyplot as plt




############## Loading data into a numpy array (ij entry corresponds to the count of flies ith bin at the jth frame of the video)##########################
import csv
datafile=sys.argv[1]
Nmat=list(np.array(list(csv.reader(open(datafile))), dtype=np.int).T)


MaxFlies=np.max(np.array(Nmat).flatten())   #maximum number of flies
Nbins=int(sys.argv[2])*int(sys.argv[2])

####################################################################################################################



##############calculating the time correlations and the time constant to get true error in the model#################
for kT in range (100,3600,100):

#substracting the mean number of flies to each bin in order to calculate the correlation function
	delta_Nmat=[[0 for i in range(kT)] for j in range(np.shape(Nmat)[0])]
	for i in range(np.shape(Nmat)[0]):
		delta_Nmat[i]=Nmat[i][:kT]-np.mean(Nmat[i][:kT])



	Time_frames=np.shape(delta_Nmat)[1]#defining a variable which corresonds to the total number of frames in the data
	TcorrConv=np.zeros(Time_frames)        # declaring array that will contain the time correlations as a function of the difference of times between frames averaged over all bins, details in the pdf that is attached to these programs

#calculating the time correlation matrix averaged over all bins as the matrix product of the counts per bin minus their average in every frame transposed which is a Time_frames x Nbins matrix and delta_Nmat which is a Nbins x Time_frames matrix. the ij entry of Tcorr is the correlations of the system between the ith and jth time frame averaged over all bins
	TCorr=np.dot(np.array(delta_Nmat).T,np.array(delta_Nmat))/float(Nbins)

	TCorr+=TCorr.T # symmetrizing the TCorr matrix to simplify further calculations
	counter=np.zeros(Time_frames)  #counter will weight the number of points that will be averaged (we will be averaging each point on a given off diagonal of the TCorr matrix, this is because the difference between the indices of a given off diagonal is constant)

	for i in range(Time_frames):
		for j in range(i,Time_frames):
			TcorrConv[j-i]+=TCorr[i,j]    #averaging all elements in the TCorr matrix that are equally spaced in time
			counter[j-i]+=2

	TcorrConv/=counter						#diving by the relevsant number of points that were averaged
	MaxT= int(2.0*np.where(TcorrConv<0.0)[0][0]/3.0)			#time frame in which the first negative entry appears, were noise takes over the exponential decay in the correlation function


	tau=np.sum(TcorrConv[:MaxT]/TcorrConv[0])  #the time constant of the system calculated as the integral of the normalized correlation function
	print float(kT), " ", tau
#########################################################################################################################



end = time.time()

print "time elapsed in this run " +str(end - start)
