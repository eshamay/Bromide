import sys
from glob import glob
from ColumnDataFile import ColumnDataFile as CDF
import matplotlib.pyplot as plt
from pylab import *
from numpy import array, interp, linspace, zeros
from operator import itemgetter
import PlotUtility

# load all of the files for the analysis and parse out the columns of data
files = glob(sys.argv[1:])
#files = glob(sys.argv[1])
cdf = [CDF(f) for f in files]

fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
axs = fig.add_subplot(1,1,1)

# create a common x axis for interpolation across all the data sets
num_columns = 3
num_rows = len(cdf[0][0])
x_min = 0.0
x_max = 40.0
xnew = linspace(x_min,x_max,num_rows*2)

# create an empty container to hold all the incoming data
data = [zeros(len(xnew)) for i in range(num_columns-1)]

# in every file, take the first column to be the x-axis.. and sort all the data by the x-axis
for c in cdf:
	c_data = [c[i] for i in c][1:]	# only grab the columns of interest
	c_data = zip(*c_data)
	c_data = sorted(c_data, key=itemgetter(0))	# sort the columns by the values in one of them
	c_data = zip(*c_data)
	c_data = [interp (xnew,c_data[0],c_data[i]) for i in range(1,num_columns)]

	for d in range(num_columns-1):
		data[d] = data[d] + c_data[d]
	
data = [data[i]/len(files) for i in range(num_columns-1)]

data.insert(0,xnew)
data = zip(*data)
for d in data:
	for i in d:
		print i,
	print


