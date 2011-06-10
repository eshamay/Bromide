import sys
import glob
from ColumnDataFile import ColumnDataFile as CDF
from scipy import interpolate
from scipy import signal
import matplotlib.pyplot as plt
from pylab import *
import numpy
import PlotUtility
import Smoothing

def PlotFiles(files):
	cdfs = [CDF(f) for f in files]
	num = len(cdfs[0][0])
	xi = numpy.array(cdfs[0][0])
	yi = numpy.zeros(num)
	
	for c in cdfs:
		yi = yi + c[1]
	
	#yi = yi / yi.max()
	#histo,edges = numpy.histogram(xi, weights=yi, normed=True, bins=400)
	
	yi = Smoothing.window_smooth(yi,window_len=10)
	axs.plot(xi,yi,linewidth=3.0)

#files_cold = glob.glob(sys.argv[1])
files_cold = glob.glob('[1-5]/'+sys.argv[1]+'*')
files_hot = glob.glob('[6-9]/'+sys.argv[1]+'*')
files_hot = files_hot + glob.glob('10/'+sys.argv[1]+'*')

fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
axs = fig.add_subplot(1,1,1)

PlotFiles (files_cold)
PlotFiles (files_hot)

axs.set_xlabel(r'SO$_2$ Orientation: $\phi$ / degrees', fontsize='64')
#axs.set_ylabel('Bondlength', fontsize='64')

xticks(fontsize=48)
yticks([])
axs.set_ylim(0,1350)

plt.show()
