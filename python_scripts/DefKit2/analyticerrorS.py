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
if len(sys.argv) != 3:
	print "this program needs (2) arguments to run"
	sys.exit()


#### importing libraries
import numpy as np
np.set_printoptions(threshold=np.inf)
from scipy import special
import matplotlib.pyplot as plt
from numpy import linalg as la
from scipy import stats
from scipy.optimize import curve_fit
from scipy import interpolate
import math
from scipy.interpolate import interp1d




############## Loading data into a numpy array (ij entry corresponds to the count of flies ith bin at the jth frame of the video)##########################
import csv
datafile=sys.argv[1]
Nmat=np.array(list(csv.reader(open(datafile), delimiter=',')), dtype=np.int).T


print np.shape(Nmat)


MaxFlies=np.max(np.array(Nmat).flatten())   #maximum number of flies
Nbins=int(np.shape(Nmat)[0])

print "the maximum number of flies per bin is ", MaxFlies, " \n"
#calculating the matrix of histograms (i j entry corresponds to the relative frequency of the occurrence of j flies in the ith bin in the video)
Sis=[]
for i in range(Nbins):
	hist, bin_edges = np.histogram(Nmat[i] ,bins=np.arange(0, MaxFlies+2 , 1), density=True)
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
tau=sum(TcorrConv[:MaxT]/TcorrConv[0])  #the time constant of the system calculated as the integral of the normalized correlation function
#########################################################################################################################






#################doing the least squares fit to get frustrations and vexations, to check how this system of linear equations arise look at the attached pdf##################


Sig=1/( special.polygamma(1,np.array(Sis)*Time_frames/tau )  - special.polygamma(1,Time_frames/tau ) )  # error matrix must be in terms of counts not normalized histogram so we multiply by Time_frames
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

A=np.bmat([[ff,Af,Vf], [fA, AA, VA],[fV,AV, VV]])

#building the solution vector
bf=np.diag(np.dot(rho.T,Sig))
bA=np.diag(np.dot(rho,Sig.T))
bV=np.diag(np.dot(rho,(ns*Sig).T))

B=np.bmat([[bf],[bA],[bV]]).T

Sol=la.lstsq(A,B)

Apinv=la.pinv(A)

############building the linear operator that acts on rho on the r.h.s to produce the b vector#####
Nrowsrho=np.shape(rho)[0]
Ncolsrho=np.shape(rho)[1]
testrho=np.mat(rho.flatten()).T   #rho goes from the vector space of matrices Nrowsrho xNcolsrho  to R^n (n=dimrho)
dimrho=np.size(testrho)


#matrix rep of the linear operation that turns rho into bf evaluating in the canonical basis of matrices the same dimension as rho
RepbF=[]
for i in range(dimrho):
	row=np.zeros(dimrho)
	row[i]=1
	RepbF.append(np.diag(np.dot(np.reshape(row,(Nrowsrho, Ncolsrho)).T, Sig)).flatten())

matbF=np.array(RepbF).T
#has=np.dot(matbF,testrho).T-np.mat(bf)
#print np.shape(matbF),np.shape(has),np.shape(np.mat(bf)),np.shape(np.dot(matbF,testrho).T),has

#matrix rep of the linear operation that turns rho into bA evaluating in the canonical basis of matrices the same dimension as rho
RepbA=[]
for i in range(dimrho):
	row=np.zeros(dimrho)
	row[i]=1
	RepbA.append(np.diag(np.dot(np.reshape(row,(Nrowsrho, Ncolsrho)),Sig.T)).flatten())

matbA=np.array(RepbA).T
#has=np.dot(matbA,testrho).T-np.mat(bA)
#print np.shape(matbA),np.shape(has),np.shape(np.mat(bA)),np.shape(np.dot(matbA,testrho).T),has

#matrix rep of the linear operation that turns rho into bV evaluating in the canonical basis of matrices the same dimension as rho
RepbV=[]
for i in range(dimrho):
	row=np.zeros(dimrho)
	row[i]=1
	RepbV.append(np.diag(np.dot(np.reshape(row,(Nrowsrho, Ncolsrho)),(ns*Sig).T)).flatten())

matbV=np.array(RepbV).T
#has=np.dot(matbV,testrho).T-np.mat(bV)
#print np.shape(matbV),np.shape(has),np.shape(np.mat(bV)),np.shape(np.dot(matbV,testrho).T),has


#full linear operator that acts on rho to produce B by gluing the above representations of the linear transformation as a direct sum and taking the product with three identity matrices stacked upon each other to fit the dimensions of the vector space from which the linear transformation departs

L1=np.bmat([[matbF,np.zeros((np.shape(matbF)[0],dimrho)),np.zeros((np.shape(matbF)[0],dimrho))],[np.zeros((np.shape(matbA)[0],dimrho)),matbA,np.zeros((np.shape(matbA)[0],dimrho))],[np.zeros((np.shape(matbV)[0],dimrho)),np.zeros((np.shape(matbV)[0],dimrho)),matbV]])

tripId=np.bmat([[np.identity(dimrho)],[np.identity(dimrho)],[np.identity(dimrho)]])
L=np.dot(L1,tripId)

#error thar comes from the propagation in a linear transformation is just the square of the coefficients

M=np.dot(Apinv, L)  #multiplying L by the moore penrose pseudoinverse
Sigmas=1/np.mat(Sig.flatten()).T  #Sig is the inverse of the squares of the errors in rho

Error=np.ravel(np.sqrt(np.dot(np.square(M),Sigmas)))  #full propagation

print np.shape(M), np.shape(Error), Error







#############################continuing to the plots of the least squares fit#######################



#plotting frustrations with errorbars and van der waals fit

F=np.ravel(Sol[0][0:np.shape(ff)[1]])
N=np.arange(np.size(F))
Ncont=np.linspace(0,np.size(F)-1,1000)
slope, intercept, r_value, p_value, std_err = stats.linregress(N[0:2], F[0:2])

print np.shape(F)

chich=0.01
Abin=float(sys.argv[2])*float(sys.argv[2])
factor=(2/np.pi)*special.ellipk(chich**2)/np.sqrt(1-chich*chich)


#fitting to the VanDer Waals functional
#def frust2(xx):
#	f= (1.5*(xx)/((1-xx/float(MaxFlies+0.01))**2) - 1.5*xx- xx*np.log (1 - xx/float(MaxFlies+0.01)))*factor
#	return f
#xnew = np.linspace(N[0]+0.1, N[-1]-0.1 , 1000)
#plt.scatter(N, F-N*slope-intercept)
#plt.plot(xnew,frust2(xnew))
#plt.show()

#def fit(xx,R):
#	f= R*(xx)**2
#	return f

def frust1(xx):
	f= - xx*np.log (1 - xx/float(MaxFlies))
	return f

xnew = np.linspace(N[0]+0.1, N[-1]-0.1 , 1000)
F2=F-N*slope-intercept
plt.scatter(N, F2)
plt.plot(xnew,frust1(xnew))
plt.show()


#def fit2(xx,R,B,C,D):
#	f= - math.fabs(R)*xx*np.log(1 - xx/(Abin*1000.0))+B*xx*xx/(1 - xx/(Abin*1000.0))+C*xx*xx*xx/(1 - xx/(Abin*1000.0))**2
#	return f
#def fit(xx,R,E,C):
#	f= - xx*np.log(1 - R*xx)+E*(xx**2)/(1 - R*xx)+C*(xx**3)/(1 - R*xx)**2
#	return f
#popt, pcov = curve_fit(frust1 ,N, F-N*slope-intercept)
#popt2, pcov2 = curve_fit(fit2 ,N, F-N*slope-intercept)
#f2 = interp1d(N, F-N*slope-intercept, kind='cubic')


plt.xlim([-0.1,np.size(F)+0.1])
#Fspl=f2(xnew)
Fspl=F
Fspl/=Abin
#en=xnew/(Abin)
en=N/(Abin)
plt.errorbar(N,F,yerr=2*Error[0:np.size(F)], fmt='o')
plt.title('Frustration with analytic error bars')
plt.ylabel('Frustration')
plt.xlabel('Number of flies in a bin')
plt.show()

plt.plot(en,Fspl,'r-', lw=2)
plt.show()

print "time constant of the system ",tau
#print "\n fit parameter for the area of a fly is ", popt2, popt



file_frust = open('frustSPL'+datafile[:-4]+'.dat', 'w')

for i in range(np.size(Fspl)):
	file_frust.write("%f %f\n"%(en[i], Fspl[i]))
file_frust.close()


#returning to the original Sig matrix

#plotting and printing vexations
Vex=np.reshape(Sol[0][(np.shape(ff)[1]+np.shape(AA)[0]):],(np.sqrt(Nbins),np.sqrt(Nbins)),order='C') #Assuming square chamber


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



#ploting vexation with errorbars
Vex=np.ravel(Vex.flatten())
print np.shape(Vex)
plt.errorbar(np.arange(len(Vex)),Vex,yerr=Error[(np.shape(ff)[1]+np.shape(AA)[0]):], fmt='o')
plt.title('Vexation with analytic error bars')
plt.ylabel('Vexation')
plt.xlabel('bin number')
plt.show()

#print a file with the vexation of the given dataset
file_vex = open('vex'+datafile, 'w')
for i in range(np.size(Vex)-1):
	s=str(np.ravel(Sol[0][np.shape(ff)[1]+np.shape(AA)[0]+i])[0])
	file_vex.write(s+',')
s=str(np.ravel(Sol[0][np.shape(ff)[1]+np.shape(AA)[0]+np.size(Vex)-1])[0])
file_vex.write(s)
file_vex.close()


#print a file with the errorbars
file_vex = open('Errorb_'+datafile, 'w')
for i in range(np.size(Error)-1):
	s=str(Error[i])
	file_vex.write(s+',')
s=str(Error[np.size(Error)-1])
file_vex.write(s)
file_vex.close()

#print a file with the errorbars for vexation alone
vexerror=Error[(np.shape(ff)[1]+np.shape(AA)[0]):]

file_vex = open('Errfitvex_'+datafile, 'w')
for i in range(np.size(vexerror)-1):
	s=str(vexerror[i])
	file_vex.write(s+',')
s=str(vexerror[np.size(vexerror)-1])
file_vex.write(s)
file_vex.close()


########################################################



end = time.time()

print "time elapsed in this run " +str(end - start)
