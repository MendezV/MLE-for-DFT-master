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
import matplotlib.pyplot as plt




############## Loading data into a numpy array (ij entry corresponds to the count of flies ith bin at the jth frame of the video)##########################
import csv
datafile=sys.argv[1]
Nmat=list(np.array(list(csv.reader(open(datafile))), dtype=np.int).T)


MaxFlies=np.max(np.array(Nmat).flatten())   #maximum number of flies
m=int(sys.argv[2])
Nbins=m*m

#substracting the mean number of flies to each bin in order to calculate the correlation function
delta_Nmat=[[0 for i in range(np.shape(Nmat)[1])] for j in range(np.shape(Nmat)[0])]
for i in range(np.shape(Nmat)[0]):
	delta_Nmat[i]=Nmat[i]-np.mean(Nmat[i])



Neigh=[]
for i in range(0,m):
	for j in range(i,m):
		Neigh.append(i**2+j**2)
Neigh=list(set(Neigh))
Neigh.sort()
Num_neigh=np.size(Neigh)


T_frames=np.shape(delta_Nmat)[1]    #defining a variable which corresonds to the total number of frames in the data
counter=np.zeros((T_frames,Num_neigh))
corrConv=np.zeros((T_frames,Num_neigh))   # declaring array that will contain the correlations as a function of the difference of times between frames averaged over all bins and the distance between them, details in the pdf that is attached to these programs


def distance(k,r,m,Neigh):
	i1=(k-k%m)/m
	i2=(r-r%m)/m
	j1=k%m
	j2=r%m
	return Neigh.index((i1-i2)**2+(j1-j2)**2)

Distances=[]
for k in range(Nbins):
	for i in range(T_frames):
		for j in range(T_frames-i):
			for r in range(Nbins):
				w=distance(k,r,m,Neigh)
				corrConv[j,w]+=(delta_Nmat[k][i])*delta_Nmat[r][i+j]
				counter[j,w]+=1
				Distances.append(w)


corrConv/=counter

hist, bin_edges = np.histogram(Distances ,bins=np.arange(0, Num_neigh+1 , 1), density=True)
hist=hist+1e-17
plt.scatter(np.arange(np.size(hist)),hist)
print np.sum(hist)

data = counter
fig, ax = plt.subplots()
heatmap = ax.pcolor(data,cmap=plt.cm.Blues)
plt.show()

data = corrConv
fig, ax = plt.subplots()
heatmap = ax.pcolor(data,cmap=plt.cm.Blues)
plt.show()

#print a file with the correlations of the given dataset
file_corr = open('Rcorr_'+datafile, 'w')
for i in range(T_frames):
	for k in range(Num_neigh-1):
		s=str(corrConv[i,k])
		file_corr.write(s+',')
	s=str(corrConv[i,Num_neigh-1])
	file_corr.write(s)
file_corr.close()

end = time.time()
print "time elapsed in this run " +str(end - start)







