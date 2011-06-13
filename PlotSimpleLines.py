import sys
import glob
from ColumnDataFile import ColumnDataFile as CDF
import matplotlib.pyplot as plt
import numpy
import operator
import PlotUtility

# load all of the files for the analysis and parse out the columns of data
files_cold = glob.glob('[1-5]/'+sys.argv[1]+'*')
files_hot = glob.glob('[6-9]/'+sys.argv[1]+'*')
files_hot = files_hot + glob.glob('10/'+sys.argv[1]+'*')

cdfs_cold = [CDF(f) for f in files_cold]
cdfs_hot = [CDF(f) for f in files_hot]
# now we assume a single x-axis and equally spaced in all the data files
x = cdfs_cold[0][0]

data_cold = [numpy.array(c[1]) for c in cdfs_cold]
data_hot = [numpy.array(c[1]) for c in cdfs_hot]

avg_data_cold = numpy.array(reduce(operator.add, data_cold))/len(data_cold)
avg_data_hot = numpy.array(reduce(operator.add, data_hot))/len(data_hot)


fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
axs = fig.add_subplot(1,1,1)
axs.plot(x,avg_data_cold)
axs.plot(x,avg_data_hot)
plt.show()
