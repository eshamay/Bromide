import sys
import glob
from ColumnDataFile import ColumnDataFile as CDF
import matplotlib.pyplot as plt
from pylab import *
import numpy
import PlotUtility

def PlotFiles (files,fig):

	cdfs = [CDF(f) for f in files]

	xmin = min(cdfs[0][0])
	xmax = max(cdfs[0][0])
	ymin = min(cdfs[0][1])
	ymax = max(cdfs[0][1])

	xnum=180
	ynum=60
	xi = linspace (xmin, xmax, num=xnum)
	yi = linspace (ymin, ymax, num=ynum)

	zi = numpy.zeros((ynum,xnum))

	for c in cdfs:
		new_data = numpy.array(griddata(c[0], c[1], c[2], xi, yi))
		zi = zi + new_data

	# normalize to max = 1, min = 0
	#zi = zi / zi.max(0)

	im = plt.imshow(zi, extent=(xmin,xmax,ymax,ymin), interpolation='bilinear', figure=fig, aspect='auto')

	return zi, im


#files = glob.glob(sys.argv[1])
files_cold = glob.glob('[1-5]/'+sys.argv[1]+'*')
files_hot = glob.glob('[6-9]/'+sys.argv[1]+'*')
files_hot = files_hot + glob.glob('10/'+sys.argv[1]+'*')

fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
zi_cold,im_cold = PlotFiles (files_cold,fig)
xticks(fontsize=48)
yticks(fontsize=48)
axs = plt.gca()
axs.set_xlabel(r'$\theta$ / degrees', fontsize='64')
axs.set_ylabel(r'O$_{SO_2}$-H$_{H_2O}$ Distance / $\AA$', fontsize='64')

fig = plt.figure(num=2, facecolor='w', edgecolor='w', frameon=True)
zi_hot,im_hot = PlotFiles (files_hot,fig)

#fig = plt.figure(num=3, facecolor='w', edgecolor='w', frameon=True)
#im = plt.imshow(zi_hot - zi_cold, extent=(0.0,180.0,3.0,0.0), interpolation='bilinear', figure=fig, aspect='auto')

xticks(fontsize=48)
yticks(fontsize=48)
axs = plt.gca()
axs.set_xlabel(r'$\theta$ / degrees', fontsize='64')
axs.set_ylabel(r'O$_{SO_2}$-H$_{H_2O}$ Distance / $\AA$', fontsize='64')


#xticks(fontsize=48)
#plt.xlim(-10.0,10.0)
#plt.xlabel (r'Distance to water slab surface / $\AA$', fontsize=64)


#name = cdfs[0].filename
#if "theta" in name:
	#plt.ylabel (r'$\cos(\theta)$', fontsize=64)
	#plt.ylim(-1.0,1.0)
	#yticks([-1.0,-0.5,0.0,0.5,1.0], fontsize=48)
#elif "phi" in name:
	#plt.ylabel (r'$\cos(\phi)$', fontsize=64)
	#plt.ylim(0.0,1.0)
	#yticks([0.0,0.25,0.5,0.75,1.0], fontsize=48)


plt.show()
