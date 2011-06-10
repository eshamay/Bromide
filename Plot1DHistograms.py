import sys
import glob
from ColumnDataFile import ColumnDataFile as CDF
from scipy import interpolate
from scipy import signal
import matplotlib.pyplot as plt
from pylab import *
import numpy
import PlotUtility

def PlotFiles(files):
	cdfs = [CDF(f) for f in files]
	num = len(cdfs[0][0])
	xi = numpy.array(cdfs[0][0])
	yi = numpy.zeros(num)
	
	for c in cdfs:
		yi = yi + c[1]
	
	#yi = yi / yi.max()
	#histo,edges = numpy.histogram(xi, weights=yi, normed=True, bins=400)
	
	axs.plot(xi,yi)

#files_cold = glob.glob(sys.argv[1])
files_cold = glob.glob('[1-5]/'+sys.argv[1]+'*')
#files_hot = glob.glob('[6-9]/'+sys.argv[1]+'*')
#files_hot = files_hot + glob.glob('10/'+sys.argv[1]+'*')

fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
axs = fig.add_subplot(1,1,1)

PlotFiles (files_cold)
#PlotFiles (files_hot)

xticks(fontsize=48)
yticks([])

plt.show()
