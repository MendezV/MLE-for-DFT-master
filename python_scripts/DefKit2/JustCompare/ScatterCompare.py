############  compares two vexations to see if they differ by scaling and a constant shift#########


################################
#		input:	-Csv file with Vexation calculated from least squares fit in least_squares_fit.py
#				-Csv file with another Vexation calculated from low density data
#
#
#		output:	-scatter plot of the two vexations with linear regression
#				-parameters of the linear regression
#				-linear plot fitting the scatter plot
#
#
#################################

##### measure the runtime of the program
import time
start = time.time()



#### outputs error if the name of the file with the data for vexations is given as argument
import sys
if len(sys.argv) != 6:
	print "this program needs (5) arguments to run"
	sys.exit()



#importing libraries
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from scipy import interpolate


#Loading Vexations into the program
import csv
datafile=sys.argv[1]
vex= list(csv.reader(open(datafile)))[0]
Vexation=np.zeros(len(vex))
for i in range(len(vex)):
	Vexation[i]=float(vex[i])

datafile2=sys.argv[2]
vex2= list(csv.reader(open(datafile2)))[0]
Vexation2=np.zeros(len(vex2))
for i in range(len(vex2)):
	Vexation2[i]=float(vex2[i])


########## error load ##########
datafileErr=sys.argv[3]
ERR= list(csv.reader(open(datafileErr)))[0]
Ebar1=np.zeros(len(ERR))
for i in range(len(ERR)):
	Ebar1[i]=float(ERR[i])


datafileErr=sys.argv[4]
ERR= list(csv.reader(open(datafileErr)))[0]
Ebar2=np.zeros(len(ERR))
for i in range(len(ERR)):
	Ebar2[i]=float(ERR[i])
Nbins=int(sys.argv[5])*int(sys.argv[5])
print np.shape(Vexation), np.shape(Vexation2)

#linear plot to the plot of vexation1 vs vexation2
slope, intercept, r_value, p_value, std_err = stats.linregress(Vexation, Vexation2)

#plt.scatter(Vexation, Vexation2)
plt.errorbar(Vexation,Vexation2, xerr=Ebar1,yerr=Ebar2, fmt='o')
plt.plot(Vexation,Vexation)
plt.savefig('vex_dif'+datafile[:-4]+'_'+datafile2[:-4], bbox_inches='tight')
plt.show()



#############plot both
plt.errorbar(np.arange(Nbins),Vexation,yerr=Ebar1, fmt='o')
plt.errorbar(np.arange(Nbins),Vexation2,yerr=Ebar2, fmt='o')
plt.show()


print " for this two vexations, the relative factor is ", slope
print "\n the r value of the fit is ", r_value
print "\n the standard error of the fit is ", std_err
print "\n the p value of the fit is ", p_value
print np.sum(Vexation),np.sum(Vexation2)

Vex_diff=slope*Vexation+intercept-Vexation2
Nbins=int(sys.argv[5])*int(sys.argv[5])


#plot difference of vexations interpolated as a heat map
x = np.arange(0, np.sqrt(Nbins), 1)
y = np.arange(0, np.sqrt(Nbins), 1)
f = interpolate.interp2d(x, y, Vex_diff, kind='cubic')
xnew = np.arange(0, np.sqrt(Nbins)-1, 1e-2)
ynew = np.arange(0, np.sqrt(Nbins)-1, 1e-2)
znew = f(xnew, ynew)
fig = plt.figure(1, figsize=(10.5,8.5))
plt.imshow(znew)
plt.savefig('vex_dif_heatmap'+datafile[:-4]+'_'+datafile2[:-4], bbox_inches='tight')
plt.show()


#print a file with the vexation of the given dataset
file_vex = open('vex_dif'+datafile[:-4]+'_'+datafile2[:-4]+'.csv', 'w')
for i in range(np.size(Vexation)-1):
	s=str(Vex_diff[i])
	file_vex.write(s+',')
s=str(Vex_diff[np.size(Vexation)-1])
file_vex.write(s)
file_vex.close()

end = time.time()

print "time elapsed in this run " +str(end - start)
