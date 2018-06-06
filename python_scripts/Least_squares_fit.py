############ makes a least squares fit to the input data for a Boltzmann distribution consistent with a non linear correction given by volume exclusion #########


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
from scipy import fftpack



############## Loading data into a numpy array (ij entry corresponds to the count of flies ith bin at the jth frame of the video)##########################
import csv
datafile=sys.argv[1]
Nmat=list(np.array(list(csv.reader(open(datafile))), dtype=np.int).T)


MaxFlies=np.max(np.array(Nmat).flatten())   #maximum number of flies
Nbins=int(sys.argv[2])*int(sys.argv[2])
				
print "the maximum number of flies per bin is ", MaxFlies, " \n"
#calculating the matrix of histograms (i j entry corresponds to the relative frequency of the occurrence of j flies in the ith bin in the video)
Sis=[]
for i in range(Nbins):
	hist, bin_edges = np.histogram(Nmat[i] ,bins=np.arange(0, MaxFlies+1 , 1), density=True)
	hist=hist+1e-17    #adding machine precision value to prevemnt infs later
	Sis.append(hist)

####################################################################################################################








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
plt.plot(np.arange(MaxT),TcorrConv[:MaxT]/TcorrConv[0])
plt.ylim([-0.1,1.1])
plt.show()

tau=sum(TcorrConv[:MaxT]/TcorrConv[0])  #the time constant of the system calculated as the integral of the normalized correlation function

#time correlations with fft
Tcorr=np.zeros(np.size(Nmat[1]))
for i in range(np.shape(Nmat)[0]):
	Tcorr += np.abs(fftpack.ifft(fftpack.fft(Nmat[i])*np.conjugate(fftpack.fft(Nmat[i]))))/(np.shape(Nmat)[0]*np.size(Nmat[1]))-np.mean(Nmat[i])**2/(np.shape(Nmat)[0])

def expfit(x,amp,tau,c):
	return (amp*np.exp(-x/tau))+c
x=np.arange(np.size(Nmat[1])/2)

popt , pcov= curve_fit(expfit , x , Tcorr[:np.size(Nmat[1])/2] )


#tau=popt[1] to switch between time constants, just uncomment this tau

plt.show()
#########################################################################################################################






#################doing the least squares fit to get frustrations and vexations, to check how this system of linear equations arise look at the attached pdf##################


Sig=np.array(Sis)*Time_frames/tau       # error matrix must be in terms of counts not normalized histogram so we multiply by Time_frames
ns=[[i for i in range(np.shape(Sig)[1])] for j in range(np.shape(Sig)[0])]
nss=[[i**2 for i in range(np.shape(Sig)[1])] for j in range(np.shape(Sig)[0])]
nfact=[[special.gamma(i+1) for i in range(np.shape(Sig)[1])] for j in range(np.shape(Sig)[0])]
rho=-np.log(nfact*np.array(Sis)/np.size(Nmat[1])+1e-17)    #adding machine precision value to prevemnt infs later

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

L=np.bmat([[ff,Af,Vf], [fA, AA, VA],[fV,AV, VV]])

#building the solution vector
bf=np.diag(np.dot(rho.T,Sig))
bA=np.diag(np.dot(rho,Sig.T))
bV=np.diag(np.dot(rho,(ns*Sig).T))

B=np.bmat([[bf],[bA],[bV]]).T

Sol=la.lstsq(L,B)
F=np.ravel(Sol[0][0:np.shape(ff)[1]])    # F corresponds to the first parameters that come from the fit

#print "\n Frustration is given by ",F

#printing Frustrations to file
file_frust = open('frust'+datafile, 'w')
for i in range(np.size(F)-1):
	s=str(F[i])
	file_frust.write(s+',')
s=str(F[np.size(F)-1])
file_frust.write(s)
file_frust.close()

#plotting frustrations with errorbars and van der waals fit

Sig[Sig<10e-10]=10e18  #taking error off the poits that do not contribute


N=np.arange(np.size(F))
Ncont=np.linspace(0,np.size(F)-1,1000)
slope, intercept, r_value, p_value, std_err = stats.linregress(N[0:2], F[0:2])


#fitting to the Van Der Waals functional
def fit(xx,A,B,C):
	f= A + B*xx - xx*np.log (1 - C*xx/1000.0+1e-17)
	return f

popt, pcov = curve_fit(fit ,N, F-N*slope-intercept)
plt.xlim([-0.1,np.size(F)+0.1])
plt.plot(Ncont,fit(Ncont,popt[0],popt[1],popt[2]),'r-', lw=2)
plt.errorbar(N,F-N*slope-intercept,yerr=np.sqrt(np.sum(1/Sig, axis=0)/float(Nbins)**2), fmt='o')
plt.show()

Abin=float(sys.argv[3])*float(sys.argv[3])
print "time constant of the system ",tau
print "\n fit parameter for the area of a fly is " , Abin*(popt[2]/(1000.0))

Sig[Sig>10e10]=10e-18    #returning to the original Sig matrix

#plotting and printing vexations
Vex=np.reshape(Sol[0][(np.shape(ff)[1]+np.shape(AA)[0]):],(np.sqrt(Nbins),np.sqrt(Nbins)),order='C') #Assuming square chamber

#print a file with the vexation of the given dataset
file_vex = open('vex'+datafile, 'w')
for i in range(np.size(Vex)-1):
	s=str(np.ravel(Sol[0][np.shape(ff)[1]+np.shape(AA)[0]+i])[0])
	file_vex.write(s+',')
s=str(np.ravel(Sol[0][np.shape(ff)[1]+np.shape(AA)[0]+np.size(Vex)-1])[0])
file_vex.write(s)
file_vex.close()
#print Vex




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
########################################################




################Calculating Chi-squared#################


#checking Chi-Squeared
Xisq=0.0
for i in range(np.shape(ff)[1]):
	for j in range(np.shape(AA)[0]):
		Xisq+=tau*Sig[j,i]*(Sol[0][j+np.shape(ff)[1]]+Sol[0][j+np.shape(ff)[0]+np.shape(AA)[0]]*i+Sol[0][i]-rho[j,i])**2/Time_frames  #divided by the number of relevant observations Time_frames/tau
print "Chi squeared " +str(np.ravel(np.ravel(Xisq))[0])


#checking Chi-Squeared without frustration
L2=np.bmat([ [AA, VA],[AV, VV]])
B2=np.bmat([[bA],[bV]]).T

Sol2=la.lstsq(L2,B2)
Xisq=0.0
for i in range(np.shape(ff)[1]):
	for j in range(np.shape(AA)[0]):
		Xisq+=tau*Sig[j,i]*(Sol2[0][j]+Sol2[0][j+np.shape(AA)[0]]*i-rho[j,i])**2/Time_frames   #divided by the number of relevant observations Time_frames/tau
print "Chi squeared without frustration " +str(np.ravel(np.ravel(Xisq))[0])

#########################################################


end = time.time()

print "time elapsed in this run " +str(end - start)