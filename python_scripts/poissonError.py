########### makes a least squares fit to the input data for a Boltzmann distribution consistent with a non linear correction given by volume exclusion #########


################################
#		input:	-csv file with the counts of flies per bin in an m x m grid at diferent frames of a video
#				the i, j entry of this matrix corresponds to the number of flies jth bin in the ith frame of the video
#				-value of m
#				-length of the side of a square bin in mm
#
#
#		output:  3 graphs
#						-frustration as a function of the number of flies in a bin
#						-frustration as a function of the density of flies in mm^2
#						-heatmap of interpolated vexaton from the fit
#						-time correlations of the system
#				 2 files
#						-csv with vexation from fit
#						-csv with frustration
#
#				  prints
#						-Time constant of the correlation function
#						-area of a fly according to the van der waals fit
#						-chi squared before and after introducing frustration to th model
#
#################################

##### measure the runtime of the program
import time
start = time.time()

#### three arguments, the name of the file of the the counts of flies per bin in an m x m grid at diferent frames of a video, the value of m and the lenght of a side of a bin in mm
import sys
if len(sys.argv) != 4:
	print "this program needs (3) arguments to run"
	sys.exit()


#### importing libraries
import numpy as np
from scipy import special
import matplotlib.pyplot as plt
from numpy import linalg as la
from scipy import stats
from scipy.optimize import curve_fit
from scipy import interpolate
import math
from scipy.interpolate import interp1d
import random




############## Loading data into a numpy array (ij entry corresponds to the count of flies ith bin at the jth frame of the video)##########################
import csv
datafile=sys.argv[1]
Nmat=list(np.array(list(csv.reader(open(datafile))), dtype=np.int).T)

Solsillas=[]

MaxFlies=np.max(np.array(Nmat).flatten())   #maximum number of flies
Nbins=int(sys.argv[2])*int(sys.argv[2])

#print "the maximum number of flies per bin is ", MaxFlies, " \n"





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



def randomMat(Sig):
	n_side_N=np.shape(Sig)[0]
	n_side_M=np.shape(Sig)[1]
	RandMat=np.zeros(np.shape(Sig))
	for i in range(n_side_N):
		for j in range(n_side_M):
			RandMat[i,j]=np.sqrt(Sig)*np.random.randn()


###########iterating the process of doing the fit to get error bars#######################

repeats=10000

for k in range (repeats):

	#############calculating the matrix of histograms (i j entry corresponds to the relative frequency of the occurrence of j flies in the ith bin in the video)##############
	Sis=[]
	for i in range(Nbins):
		hist, bin_edges = np.histogram(Nmat[i] ,bins=np.arange(0, MaxFlies+1 , 1), density=True)
		for j in range(np.size(hist)):
			a=hist[j]
			if( a >= 0.0):
				hist[j]=a+1e-17
		Sis.append(hist)
#print np.sum(hist)

	####################################################################################################################







	#################doing the least squares fit to get frustrations and vexations, to check how this system of linear equations arise look at the attached pdf##################


	Sig=np.array(Sis)*Time_frames/tau        # error matrix must be in terms of counts not normalized histogram so we multiply by Time_frames
	ns=[[i for i in range(np.shape(Sig)[1])] for j in range(np.shape(Sig)[0])]
	nss=[[i**2 for i in range(np.shape(Sig)[1])] for j in range(np.shape(Sig)[0])]
	nfact=[[special.gamma(i+1) for i in range(np.shape(Sig)[1])] for j in range(np.shape(Sig)[0])]
	Sig[Sig<10E-5]=10E15
	rho=-np.log(nfact*np.array(Sis)/np.size(Nmat[1])+1e-17) +np.random.normal(0,np.sqrt(1/Sig),(np.shape(Sig)[0],np.shape(Sig)[1]))   #adding machine precision value to prevemnt infs later
	Sig[Sig>10E10]=0
	#building the least squeares matrix
	ff=np.diag(np.sum(Sig, axis=0))
	fA=Sig
	fV=Sig*ns
	AA=np.diag(np.sum(Sig, axis=1))
	AV=np.diag(np.sum(ns*Sig, axis=1))
	VV=np.diag(np.sum(nss*Sig, axis=1))
	Af=fA.T
	Vf=fV.T
	VA=AV.T

	A=np.bmat([[ff,Af,Vf], [fA, AA, VA],[fV,AV, VV]])

	#building the solution vector
	bf=np.diag(np.dot(rho.T,Sig))
	bA=np.diag(np.dot(rho,Sig.T))
	bV=np.diag(np.dot(rho,(ns*Sig).T))

	B=np.bmat([[bf],[bA],[bV]]).T

	Sol=la.lstsq(A,B)

	Solsillas.append(np.ravel(Sol[0][0:]))



for k in range(repeats):
	plt.scatter(np.arange(np.size(Solsillas[k][(np.shape(ff)[1]+np.shape(AA)[0]):])),Solsillas[k][(np.shape(ff)[1]+np.shape(AA)[0]):],s=0.1)
plt.title('Vexation with gaussian noise')
plt.ylabel('Vexation')
plt.xlabel('bin number')
plt.show()


FR=[]
for k in range(repeats):
	F=np.ravel(Solsillas[k][0:np.shape(ff)[1]])
	N=np.arange(np.size(F))
	Ncont=np.linspace(0,np.size(F)-1,1000)
	slope, intercept, r_value, p_value, std_err = stats.linregress(N[0:2], F[0:2])
	plt.scatter(N,F,s=0.1)
	#plt.scatter(N,F-N*slope-intercept,s=0.1)
	FR.append(F)
	#FR.append(F-N*slope-intercept)
plt.title('Frustration with gaussian noise')
plt.ylabel('Frustration')
plt.xlabel('Number of flies in a bin')
plt.show()

Means=np.mean(Solsillas,axis=0)
STDS=np.std(Solsillas,axis=0)
slope, intercept, r_value, p_value, std_err = stats.linregress(N[0:2], Means[0:2])
print Means ,STDS

F=Means[0:np.shape(ff)[1]]
plt.errorbar(N,F-N*slope-intercept,yerr=2*STDS[0:np.size(F)], fmt='o')
plt.title('Frustration with gaussian noise')
plt.ylabel('Frustration')
plt.xlabel('Number of flies in a bin')
plt.xlim( (-0.1, np.shape(ff)[1]+0.1) )
plt.show()

Fmea=np.mean(FR,axis=0)
FRer=np.std(FR,axis=0)
plt.errorbar(N,Fmea,yerr=2*FRer, fmt='o')
plt.title('Frustration with gaussian noise')
plt.ylabel('Frustration')
plt.xlabel('Number of flies in a bin')
plt.xlim( (-0.1, np.shape(ff)[1]+0.1) )
plt.show()

Vex=Means[(np.shape(ff)[1]+np.shape(AA)[0]):]
plt.errorbar(np.arange(len(Vex)),Vex,yerr=2*STDS[(np.shape(ff)[1]+np.shape(AA)[0]):], fmt='o')
plt.title('Vexation with gaussian noise')
plt.ylabel('Vexation')
plt.xlabel('bin number')
plt.show()





end = time.time()

print "time elapsed in this run " +str(end - start)