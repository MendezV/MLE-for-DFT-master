############ predicts the average number of flies per bin in a m x m grid by explicitly calculating the average for the fruit fly system#########


################################
#		input:	-Csv file with Vexation
#				-File 2 columns, frustration
#				-number of flies in the experiment (65 for 2.csv)
#				-maximum number of flies per bin in the experiment, for faster calculation of the chemical potential
#
#		output:	-csv file with predicted average number of flies n each bin
#				-prints the predicted number of flies in each bin,
#				-the chemical potential inferred from the data ,
#				-the total number of flies in the experiment and the number of steps that it took for the binary search to find the chemical potential
#################################

##### measure the runtime of the program
import time
start = time.time()



#### outputs error if the name of the file with the data for vexations is given as argument
import sys
if len(sys.argv) != 5:
	print "this program needs (4) arguments to run"
	sys.exit()



#importing libraries
import numpy as np
from scipy.optimize import brentq
from scipy import special
from scipy.interpolate import interp1d


#Loading Vexations into the program
import csv
datafile=sys.argv[1]
vex= list(csv.reader(open(datafile)))[0]
Vexation=np.zeros(len(vex))
for i in range(len(vex)):
	Vexation[i]=float(vex[i])

#importing frustrations

datos = np.loadtxt(sys.argv[2])

F = interp1d(datos[:,0], datos[:,1], kind='cubic')
delta=datos[1,0]-datos[0,0]


#remu=real mu of the system that gives the correct number of flies Nflies for which the parameter alpha defines the ratio between the area of a fly and the area of a bin
remu=0.0



Nflies=float(sys.argv[3])
Nbins=np.size(Vexation) #total number of bins


#Defining the equation that needs to be solved in order to minimize the Energy functional for all boxes
def g(n,mu,V):
	return -mu + V + (F(n/10.0+delta)-F(n/10.0))/delta +special.digamma(1+n)

def g2(n,mu,V):
	return (-mu + V)*n + F(n/10.0) +special.digamma(1+n)
#array that will contain the predicted mean number of flies per bin
nb=np.zeros(Nbins)


########## binary search for the chemical potential and brent method for the average number of flies per bin #########

NmaxSteps=130    #number of steps before giving up on the search
tolerance= 1e-7    #tolerance in the difference between the actual number of flies and the one that comes up with the guessed mu
counter=0
#a=5 #guess for the highest lower bound for the chemical potential
#b=10   #guess for the lowest higher bound for the chemical potential
mu=np.max(Vexation)
NfliesGuess=0.0
a=mu+0.5
b=mu-0.5
MaxFlies=float(sys.argv[4])

ns=np.arange(0,MaxFlies,1)+0.1
#performing the binary search
while(counter<NmaxSteps):
	
	mu=(a+b)/2.0
	
	
	#solving the equation for each bin with the brent method fixing mu asuming that the average lies somewhere between 0 flies and the maximum number of flies per bin (6 in the case of 2.csv)
	for i in range(Nbins):
		print np.sign(g(0.5,mu,Vexation[i])),np.sign(g(MaxFlies,mu,Vexation[i])),g(0.5,mu,Vexation[i]),g(MaxFlies,mu,Vexation[i]),mu
		
		#nb[i] =brentq(g, 0.5, MaxFlies ,args=(mu,Vexation[i]))
		nb[i]=np.sum(ns*np.exp(-g2(ns,mu,Vexation[i])))/np.sum(np.exp(-g2(ns,mu,Vexation[i])))
	
	
	NfliesGuess=np.sum(nb)  #checking if the inferred average sums to the total number of flies, if not continue searching with a different mu
	if (NfliesGuess==Nflies)|(np.abs(NfliesGuess-Nflies)<tolerance) :
		remu=mu
		break
	counter+=1
	print NfliesGuess
	#adjusting the chemical potential acordingly, if NfliesGuess is higher than Nflies we need a lower mu, if NfliesGuess is lower than Nflies we need a higher mu
	if(NfliesGuess<Nflies):
		b=mu
		print a,b
	else:
		a=mu
		print a,b


print nb, "\n the chemical potential is ", remu, "\n and the total number of flies is ", np.sum(nb), "\n number of steps was in binary search was ", counter

#print a file with the average of the given dataset
file_nb = open('n_Dft_'+datafile, 'w')
for i in range(np.size(nb)-1):
	s=str(nb[i])
	file_nb.write(s+',')
s=str(nb[np.size(nb)-1])
file_nb.write(s)
file_nb.close()

end = time.time()

print "time elapsed in this run " +str(end - start)




