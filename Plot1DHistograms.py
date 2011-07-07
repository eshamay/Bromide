import sys
import glob
from ColumnDataFile import ColumnDataFile as CDF
from ColumnDataFile import CDFGroup as CGroup
from scipy import interpolate
from scipy import signal
import matplotlib.pyplot as plt
from pylab import *
import numpy
import PlotUtility
import Smoothing

def PlotFiles(files,lbl=""):

	group = CGroup(files)
	xi,yi = group.Reduce1D(1)
	'''
	cdfs = [CDF(f) for f in files]
	num = len(cdfs[0][0])
	xi = numpy.array(cdfs[0][0])
	yi = numpy.zeros(num)
	
	for c in cdfs:
		yi = yi + c[1]
	
	'''
	yi = yi / group.num
	#yi = yi / yi.sum()
	#yi = yi / yi.max()
	#histo,edges = numpy.histogram(xi, weights=yi, normed=True, bins=400)
	
	#yi = Smoothing.window_smooth(yi,window_len=8)
	axs.plot(xi,yi,linewidth=3.0,label=lbl)

#files_cold = glob.glob(sys.argv[1])
files_cold = glob.glob('[1-5]/'+sys.argv[1]+'*')
files_hot = glob.glob('[6-9]/'+sys.argv[1]+'*')
files_hot = files_hot + glob.glob('10/'+sys.argv[1]+'*')

fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
axs = fig.add_subplot(1,1,1)

for arg in sys.argv[1:]:
	files_cold = glob.glob('[1-5]/'+arg+'*')
	files = glob.glob('[6-9]/'+arg+'*')
	files = files + glob.glob('10/'+arg+'*')
	PlotFiles (files_cold, "Cold")
	PlotFiles (files, "Hot")
#PlotFiles (files_hot)

axs.set_xlabel(r'$\phi$', fontsize='64')
axs.set_ylabel('', fontsize='64')

xticks(fontsize=48)
yticks(fontsize=48)
axs.set_xlim(0.0,7.0)

PlotUtility.ShowLegend(axs)
plt.show()
