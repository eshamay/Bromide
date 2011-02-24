import sys, glob, numpy, PlotUtility, Smoothing
from ColumnDataFile import ColumnDataFile as CDF
from scipy import interpolate
import matplotlib.pyplot as plt
from pylab import *




file = sys.argv[1]
cdf = CDF(file)

xi = numpy.array(cdf[0])
yi = numpy.array(cdf[1])
yi = yi / yi.sum()
#histo,edges = numpy.histogram(xi, weights=yi, normed=True, bins=400)

fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
axs = fig.add_subplot(1,1,1)
axs.plot(xi,yi,linestyle=':')
	
# cubic fitting
'''
tck = interpolate.splrep(xi,yi,s=0)
xnew = numpy.arange(min(xi),max(xi),0.01)
ynew = interpolate.splev(xnew,tck,der=0)
'''
ynew = smooth(yi,window_len=25,window='hamming')
axs.plot(xi,ynew,linewidth=4,color='k')

xticks(fontsize=48)
plt.xlim(-1.0,1.0)

yticks([])

name = cdf.filename
if "theta" in name or "oh" in name:
	plt.xlabel (r'$\cos(\theta)$', fontsize=64)
elif "phi" in name:
	plt.xlabel (r'$\cos(\phi)$', fontsize=64)


'''
# also plot out a 1d version that reduces all the data across the rows or columns
fig_1d = plt.figure(num=2, facecolor='w', edgecolor='w', frameon=True)
axs = fig_1d.add_subplot(1,1,1)
data_1d = zi.sum(0)
# normalize it all to a max value of 1
#data_1d = data_1d / data_1d.max()
if axis == 1:
	axs.plot(yi,data_1d)
else:
	axs.plot(xi,data_1d)

# also save the data so it can be read in with grace
numpy.savetxt(sys.argv[1]+".1d.txt", (xi,ynew), fmt="%12.6e")
'''

import csv
writer = csv.writer(open(sys.argv[1]+".1d.txt", "wb"),delimiter=' ')
writer.writerows(zip(xi,ynew))


plt.show()
