import sys
import glob
from ColumnDataFile import ColumnDataFile as CDF
import matplotlib.pyplot as plt
from pylab import *
import numpy
import operator
import PlotUtility

#new_bins = []
#for i in range(4):
#  for j in range(4):
#new_bins.append(i+10*j)
#new_bins.sort()
#N = len(new_bins)
ind = range(7)
#ind = numpy.arange(N)[1:]
width = 0.35

#files_cold = glob.glob(sys.argv[1])
files_cold = glob.glob('[1-5]/'+sys.argv[1]+'*')
files_hot = glob.glob('[6-9]/'+sys.argv[1]+'*')
files_hot = files_hot + glob.glob('10/'+sys.argv[1]+'*')

fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)


cdfs = [CDF(f) for f in files_cold]
data = [c[0] for c in cdfs]
data = reduce(operator.add, data)
data = [(int(i) % 10) + (int(i)/10) for i in data]

cold_histo, cold_bin_edges = numpy.histogram (data, bins=ind)
plt.bar(ind[:-1], cold_histo, width, color='b')

cdfs = [CDF(f) for f in files_hot]
data = [c[0] for c in cdfs]
data = reduce(operator.add, data)
data = [(int(i) % 10) + (int(i)/10) for i in data]

hot_histo, hot_bin_edges = numpy.histogram (data, bins=ind)
plt.bar([i+width for i in ind[:-1]], hot_histo, width, color='r')


xticks(arange(len(ind)), ind)

plt.show()
