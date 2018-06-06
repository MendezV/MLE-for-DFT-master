############ makes a least squares fit to the input data for a Boltzmann distribution consistent with a non linear correction given by volume exclusion #########


################################
#		input:	-csv file with the counts of flies per bin in an m x m grid at diferent frames of a video
#				the i, j entry of this matrix corresponds to the number of flies jth bin in the ith frame of the video
#				-value of m
#
#
#		output:  3 graphs
#						-frustration as a function of the number of flies in a bin
#						-frustration as a function of the density of flies in mm^2
#						-heatmap of interpolated vexaton from the fit
#						-time correlations of the system
#				 2 files
#						-csv with four columns: kx,ky,kx^2+ky^2,tau
#
#
#
#################################

##### measure the runtime of the program
import time
start = time.time()

#### three arguments, the name of the file of the the counts of flies per bin in an m x m grid at diferent frames of a video, the value of m and the lenght of a side of a bin in mm
import sys
if len(sys.argv) != 2:
	print "this program needs (1) arguments to run"
	sys.exit()


#### importing libraries
import numpy as np
from scipy import special
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import fftpack
from scipy import stats



############## Loading data into a numpy array (ij entry corresponds to the count of flies ith bin at the jth frame of the video)##########################
import csv
datafile=sys.argv[1]
Nmat=np.array(list(csv.reader(open(datafile))), dtype=np.int)
T_frames=np.shape(Nmat)[0]
Side=int(np.sqrt(np.shape(Nmat)[1]))
Nbins=Side*Side
print Nbins, T_frames


####2d fast fourier transform of each time frame and index extraction####
Kmat=[]
for i in range(T_frames):
	Nmat_prime=np.reshape(Nmat[i],(Side,Side),order='C')
	Nmat_ft = fftpack.fft2(Nmat_prime)
	Kmat.append(Nmat_ft.flatten())
FreqCompRows = np.fft.fftfreq(Nmat_ft.shape[0],d=2) #frequency index in the transformed index (ky)
FreqCompCols = np.fft.fftfreq(Nmat_ft.shape[1],d=2) #frequency index in the transformed index (kx)

###########################################






####### calculating the correlation functions with FFt and convolution, relating tau with the norm of the k vector #########

Kmat_prime=list(np.array(Kmat).T)
delta_Kmat=[[0 for i in range(T_frames)] for j in range(Nbins)]
for i in range(Nbins):
	delta_Kmat[i]=Kmat_prime[i]-np.mean(Kmat_prime[i])


nus1=np.zeros(Nbins) ### nus calsulated using convolution
k_squares=np.zeros(Nbins)

def autocorr(x):
	result = np.correlate(x,x, mode='full')
	return result[result.size/2:]

for k in range(Nbins):
	TcorrConv=np.zeros(T_frames)
	counter=np.zeros(T_frames)  #counter will weight the number of points that will be averaged
	
	#######with convolution#####
	
	TcorrConv=np.abs(autocorr(delta_Kmat[k]))

	#diving by the relevsant number of points that were averaged
	MaxT=30
	#MaxT=T_frames
	tau1=sum(TcorrConv[:MaxT]/(TcorrConv[0]))  #the time constant of the system calculated as the integral of the normalized correlation function
	k_squares[k]=(FreqCompCols[k%Side])**2+(FreqCompRows[(k-k%Side)/Side])**2
	nus1[k]=1/tau1
	plt.title("logcorrelation vs time(frames)")
	plt.scatter(np.arange(np.size(TcorrConv[:MaxT])),np.log(TcorrConv[:MaxT]/TcorrConv[0]),c=np.random.rand(3,1))
	print nus1[k],k_squares[k]
#   print TcorrConv[:MaxT]
#	print TcorrConv[:MaxT]/TcorrConv[0]

	############################


plt.show()


### plots##
plt.scatter(nus1,k_squares)
plt.title("$k^{2}$ vs $\nu$")
plt.show()
print nus1,k_squares
#x=np.linspace(0,2,1000)
#slope, intercept, r_value, p_value, std_err = stats.linregress(nus2[1:], k_squares[1:])
#plt.scatter(nus2[1:],k_squares[1:])
#plt.plot(x,x*slope+intercept)
#plt.show()

#print r_value,slope

end = time.time()

print "time elapsed in this run " +str(end - start)



