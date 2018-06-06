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
if len(sys.argv) != 8:
	print "this program needs (7) arguments to run"
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

datafile3=sys.argv[3]
vex3= list(csv.reader(open(datafile3)))[0]
Vexation3=np.zeros(len(vex3))
for i in range(len(vex3)):
	Vexation3[i]=float(vex3[i])


########## error load ##########
datafileErr=sys.argv[4]
ERR= list(csv.reader(open(datafileErr)))[0]
Ebar1=np.zeros(len(ERR))
for i in range(len(ERR)):
	Ebar1[i]=float(ERR[i])


datafileErr=sys.argv[5]
ERR= list(csv.reader(open(datafileErr)))[0]
Ebar2=np.zeros(len(ERR))
for i in range(len(ERR)):
	Ebar2[i]=float(ERR[i])

datafileErr=sys.argv[6]
ERR= list(csv.reader(open(datafileErr)))[0]
Ebar3=np.zeros(len(ERR))
for i in range(len(ERR)):
	Ebar3[i]=float(ERR[i])

Nbins=int(sys.argv[7])*int(sys.argv[7])
print np.shape(Vexation), np.shape(Vexation2)

#linear plot to the plot of vexation1 vs vexation2
slope, intercept, r_value, p_value, std_err = stats.linregress(Vexation, Vexation2)

#plt.scatter(Vexation, Vexation2)
plt.errorbar(Vexation,Vexation2, xerr=Ebar1,yerr=Ebar2, fmt='o')
plt.plot(Vexation,Vexation)
plt.savefig('vex_dif'+datafile[:-4]+'_'+datafile2[:-4], bbox_inches='tight')
plt.show()


#############plot both########
print (list(set(Vexation))),(list(set(Vexation2))),len(list(set(Ebar1))),len(list(set(Ebar2)))

def f7(seq):
	seen = set()
	seen_add = seen.add
	return [x for x in seq if not (x in seen or seen_add(x))]

N = len(list(set(Vexation)))
menMeans = f7(Vexation)
menStd = f7(Ebar1)

ind = np.arange(N)  # the x locations for the groups
width = 0.35       # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(ind, menMeans, width, color='r', yerr=menStd)

womenMeans = f7(Vexation2)
womenStd =  f7(Ebar2)
rects2 = ax.bar(ind + width, womenMeans, width, color='y', yerr=womenStd)

transMeans = f7(127*Vexation3)
transStd =  f7(127*Ebar3)
rects3 = ax.bar(ind + 2*width, transMeans, width, color='g', yerr=transStd)


# add some text for labels, title and axes ticks
ax.set_ylabel('Mean Number of flies')
ax.set_title('DFT Prediction for 1 to 127 flies')
ax.set_xticks(ind + width)
ax.set_xticklabels(('1', '2', '3', '4', '5', '6', '7', '8', '9', '10'))

ax.legend((rects1[0], rects2[0],rects3[0]), ('Measured', 'Predicted','Naive'))


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
