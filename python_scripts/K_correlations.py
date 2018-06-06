############ makes a least squares fit to the input data for a Boltzmann distribution consistent with a non linear correction given by volume exclusion #########


################################
#		input:	-csv file with the counts of flies per bin in an m x m grid at diferent frames of a video
#				the i, j entry of this matrix corresponds to the number of flies jth bin in the ith frame of the video
#				-value of m
#
#
#		output:  2 graphs
#						-
#				 1 files
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
if len(sys.argv) != 3:
	print "this program needs (2) arguments to run"
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
Side=int(sys.argv[2])
Nbins=Side*Side



####2d fast fourier transform of each time frame and index extraction####
Kmat=[]
for i in range(T_frames):
	Nmat_prime=np.reshape(Nmat[i],(Side,Side),order='C')
	Nmat_ft = fftpack.fft2(Nmat_prime)
	Kmat.append(Nmat_ft.flatten())
FreqCompRows = np.fft.fftfreq(Nmat_ft.shape[0],d=2) #frequency index in the transformed index (ky)
FreqCompCols = np.fft.fftfreq(Nmat_ft.shape[1],d=2) #frequency index in the transformed index (kx)

###########################################



#testing the fft2d
#psd2D= np.abs(Nmat_ft)**2
#fig = plt.figure(1, figsize=(10,10))
#plt.imshow(np.log(psd2D[:,:]))
#plt.show()




########## Exponential fit for the FFT method of calculating tau#####
def expfit(x,amp,tau,c):
	return (amp*np.exp(-x/tau))+c
##########



####### calculating the correlation functions with FFt and convolution, relating tau with the norm of the k vector #########

Kmat_prime=list(np.array(Kmat).T)
delta_Kmat=[[0 for i in range(T_frames)] for j in range(Nbins)]
for i in range(Nbins):
	delta_Kmat[i]=Kmat_prime[i]-np.mean(Kmat_prime[i])


nus1=np.zeros(Nbins) ### nus calsulated using convolution
nus2=np.zeros(Nbins) ### nus calculated using fft
k_squares=np.zeros(Nbins)


for k in range(Nbins):
	TcorrConv=np.zeros(T_frames)
	TcorrFFT=np.zeros(T_frames)
	counter=np.zeros(T_frames)  #counter will weight the number of points that will be averaged
	
	
	
	#######with convolution#####
	for i in range(T_frames):
		for j in range(T_frames-i):
			TcorrConv[j]+=np.conjugate(delta_Kmat[k][i])*delta_Kmat[k][i+j]
			counter[j]+=1
	TcorrConv=np.real(TcorrConv)
	TcorrConv/=counter						#diving by the relevsant number of points that were averaged
	MaxT=6
	tau1=sum(TcorrConv[:MaxT]/TcorrConv[0])  #the time constant of the system calculated as the integral of the normalized correlation function
	k_squares[k]=(FreqCompCols[k%Side])**2+(FreqCompRows[(k-k%Side)/Side])**2
	nus1[k]=1/tau1
	print nus1[k],k_squares[k]
	############################


	#######test plots and print statements###########
	#print TcorrConv/TcorrConv[0]
	#Plotting the correlation function
	#plt.plot(np.arange(530),(TcorrConv[:530])/TcorrConv[0])
	#plt.ylim([-0.1,1.1])
	#plt.show()
	#MaxT= int(np.where(TcorrConv<0.0)[0][0])			#time frame in which the first negative entry appears, were noise takes over the exponential decay in the correlation function
	#print MaxT

	#######################



	########withfft##########
	MaxTFFT=6  ### this is ignoring to much data for me personally
	TcorrFFT += np.abs(fftpack.ifft(fftpack.fft(delta_Kmat[k])*np.conjugate(fftpack.fft(delta_Kmat[k]))))/T_frames
	x=np.arange(MaxTFFT)
	popt , pcov= curve_fit(expfit , x , TcorrFFT[:MaxTFFT] )
	tau2=popt[1] #to switch between time constants, just uncomment this tau
	nus2[k]=1/tau2
	print nus2[k],k_squares[k]
	#####################


### plots##
plt.scatter(nus1,k_squares)
plt.show()

x=np.linspace(0,2,1000)
slope, intercept, r_value, p_value, std_err = stats.linregress(nus2[1:], k_squares[1:])
plt.scatter(nus2[1:],k_squares[1:])
plt.plot(x,x*slope+intercept)
plt.show()

print r_value,slope

end = time.time()

print "time elapsed in this run " +str(end - start)



