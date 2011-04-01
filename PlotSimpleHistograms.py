import sys
from glob import glob
from ColumnDataFile import ColumnDataFile as CDF
import matplotlib.pyplot as plt
import matplotlib as mpl
from pylab import *
from numpy import array, interp, linspace, zeros, arange
from operator import itemgetter, add
import PlotUtility

def FileData(path):

	# load all of the files for the analysis and parse out the columns of data
	files = glob(path+'*')
	cdfs = [CDF(f) for f in files]

	# grab the common x axis
	x_axis = cdfs[0][0]

	# sum all the values of the histograms
	cdf_data = [array(c[1]) for c in cdfs]
	y_axis = reduce(add, cdf_data)
	
  	return (x_axis, y_axis)
	
parms = mpl.figure.SubplotParams(bottom=0.15)
fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True, subplotpars=parms)
axs = fig.add_subplot(1,1,1)

O2S = FileData (sys.argv[1])
OSO = FileData (sys.argv[2])
CYCLE = FileData (sys.argv[3])

axs.plot(O2S[0], O2S[1], label=r'O$_2$S - OH$_2$', linewidth=2.5, linecolor='r')
axs.plot(OSO[0], OSO[1], label=r'OSO - HOH', linewidth=2.5, linecolor='b')
axs.plot(OSO[0], OSO[1], label=r'3-Water Cycle', linewidth=2.5, linecolor='k')

x_max = 10.0
x_min = -10.0
axs.set_xlabel(r'Distance Above Surface / $\AA$', fontsize='64')
xticks(arange(x_min,x_max,2.5),fontsize=48)
xlim(x_min, x_max)

yticks([])

PlotUtility.ShowLegend(axs)
plt.show()
