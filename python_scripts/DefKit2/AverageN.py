

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







##############calculating the time correlations and the time constant to get true error in the model#################


#substracting the mean number of flies to each bin in order to calculate the correlation function
delta_Nmat=[[0 for i in range(np.shape(Nmat)[1])] for j in range(np.shape(Nmat)[0])]
for i in range(np.shape(Nmat)[0]):
	delta_Nmat[i]=Nmat[i]-np.mean(Nmat[i])



Time_frames=np.shape(delta_Nmat)[1]    #defining a variable which corresonds to the total number of frames in the data
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
tau=sum(TcorrConv[:MaxT]/TcorrConv[0])  #the time constant of the system calculated as the integral of the normalized correlation function
#########################################################################################################################









print Average, np.sum(Average)


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
	#NERROR[i]=np.sqrt(np.sum(np.square(ns*sigmaDiri(hist,Tframes))))  #good ones
	NERROR[i]=np.std(Nmat[i])/np.sqrt(Tframes/tau)



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

