import time
start = time.time()

#### three arguments, the name of the file of the the counts of flies per bin in an m x m grid at diferent frames of a video, the value of m and the lenght of a side of a bin in mm



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
from scipy.optimize import brentq
nside=30
V=np.zeros((nside,nside))

for i in range(nside):
	for j in range(nside):
		V[i,j]=-(i**2+j**2)/100.0

MaxFlies=5.0
alpha=1/MaxFlies

def f(n,alpha):
	return -n*np.log(1-alpha*n)

ns=np.arange(MaxFlies)

def g(n,alpha,V):
	return  V +alpha*n/(1-alpha*n)-np.log(1-alpha*n)+special.digamma(1+n)
def Prob(N,V):
	return np.exp(-V*N-f(N,alpha))/(special.gamma(1+N)*np.sum(np.exp(-V*ns-f(ns,alpha))/special.gamma(1+ns)))


Average1=[]
Average2=[]

print np.sum(ns*Prob(ns,-3)),brentq(g, 0, MaxFlies-0.1 ,args=(alpha,-3))
for Vi in V.flatten():
	Average1.append(np.sum(ns*Prob(ns,Vi)))
	Average2.append(brentq(g, 0, MaxFlies-0.1 ,args=(alpha,Vi)))

plt.scatter(Average1,Average2)
plt.title('Maximum of the distribution Vs Average for several V')
plt.ylabel('Maximum')
plt.xlabel('Average')
plt.plot(np.linspace(0,MaxFlies-0.1,1000),np.linspace(0,MaxFlies-0.1,1000))
plt.show()

plt.plot(np.linspace(0,MaxFlies-0.1,1000),Prob(np.linspace(0,MaxFlies-0.1,1000),-16))
plt.title('Distribution for N_b with V_b=-4')
plt.ylabel('P_b(N)')
plt.xlabel('N')
plt.show()



end = time.time()

print "time elapsed in this run " +str(end - start)