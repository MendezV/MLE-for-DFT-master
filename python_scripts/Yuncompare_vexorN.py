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
if len(sys.argv) != 3:
	print "this program needs (2) arguments to run"
	sys.exit()



#importing libraries
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from scipy import interpolate

#[1 2 3 4 9 10 11 17 18 25]
Yunindex=[0,1,2,3,8,9,10,16,17,24]
#Loading Vexations into the program
import csv
datafile=sys.argv[1]
vex= list(csv.reader(open(datafile)))[0]
Vexation=np.zeros(len(Yunindex))
for i in Yunindex:
	Vexation[Yunindex.index(i)]=float(vex[i])

datafile2=sys.argv[2]
vex2= list(csv.reader(open(datafile2)))[0]
Vexation2=np.zeros(len(Yunindex))
for i in Yunindex:
	Vexation2[Yunindex.index(i)]=float(vex2[i])

print np.shape(Vexation), np.shape(Vexation2)

#linear plot to the plot of vexation1 vs vexation2
slope, intercept, r_value, p_value, std_err = stats.linregress(Vexation, Vexation2)

plt.scatter(Vexation, Vexation2)
plt.plot(Vexation,Vexation)
plt.savefig('vex_dif'+datafile[:-4]+'_'+datafile2[:-4], bbox_inches='tight')
plt.show()

print " for this two vexations, the relative factor is ", slope
print "\n the r value of the fit is ", r_value
print "\n the standard error of the fit is ", std_err
print "\n the p value of the fit is ", p_value
print np.sum(Vexation),np.sum(Vexation2)

Vex_diff=slope*Vexation+intercept-Vexation2
Nbins=np.shape(Vexation)[0]


#plot difference of vexations interpolated as a heat map
x = np.arange(0, Nbins, 1)
plt.scatter(x,Vexation)
plt.scatter(x,Vexation2,color='r')
plt.ylim(0, 10)
plt.show()


#print a file with the vexation of the given dataset
file_vex = open('vex_dif'+datafile[:-4]+'_'+datafile2[:-4]+'.csv', 'w')
for i in range(np.size(Vexation)-1):
	s=str(Vex_diff[i])
	file_vex.write(s+',')
s=str(Vex_diff[np.size(Vexation)-1])
file_vex.write(s)
file_vex.close()

############plot both########
print (list(set(Vexation))),(list(set(Vexation2)))

def f7(seq):
	seen = set()
	seen_add = seen.add
	return [x for x in seq if not (x in seen or seen_add(x))]

N = len(list(set(Vexation)))
menMeans = f7(Vexation)
menStd = f7(Vexation*0.01)

ind = np.arange(N)  # the x locations for the groups
width = 0.35       # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(ind, menMeans, width, color='r', yerr=menStd)

womenMeans = f7(Vexation2)
womenStd =  f7(Vexation2*0.01)
rects2 = ax.bar(ind + width, womenMeans, width, color='y', yerr=womenStd)

# add some text for labels, title and axes ticks
ax.set_ylabel('Mean Number of flies')
ax.set_title('DFT Prediction for 1 to 127 flies')
ax.set_xticks(ind + width)
ax.set_xticklabels(('1', '2', '3', '4', '5', '6', '7', '8', '9', '10'))

ax.legend((rects1[0], rects2[0]), ('Measured', 'Predicted'))


#def autolabel(rects):
# attach some text labels
#	for rect in rects:
#	height = rect.get_height()
#	ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
#		'%f' % float(height),
#		ha='center', va='bottom')

#autolabel(rects1)
#autolabel(rects2)

plt.show()



#############################


end = time.time()

print "time elapsed in this run " +str(end - start)
