import sys
import glob
from ColumnDataFile import ColumnDataFile as CDF
import matplotlib.pyplot as plt
from pylab import *
import numpy
import PlotUtility

#files = glob.glob(sys.argv[1])
#files = glob.glob('[1-5]/'+sys.argv[1]+'*')
files = glob.glob('[6-9]/'+sys.argv[1]+'*')
files = files + glob.glob('10/'+sys.argv[1]+'*')
cdfs = [CDF(f) for f in files]

xmin = min(cdfs[0][0])
xmax = max(cdfs[0][0])
ymin = min(cdfs[0][1])
ymax = max(cdfs[0][1])

xnum=100
ynum=100
xi = linspace (xmin, xmax, num=xnum)
yi = linspace (ymin, ymax, num=ynum)

zi = numpy.zeros((ynum,xnum))

for c in cdfs:
	new_data = numpy.array(griddata(c[0], c[1], c[2], xi, yi))
	zi = zi + new_data

# normalize to max = 1, min = 0
#zi = zi / zi.max(0)

fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
im = plt.imshow(zi, extent=(xmin,xmax,ymax,ymin), interpolation='bilinear', figure=fig, aspect='auto')
xticks(fontsize=48)
yticks(fontsize=48)
axs = plt.gca()
axs.set_xlabel(r'$\theta$ / degrees', fontsize='64')
axs.set_ylabel(r'$\phi$ / degrees', fontsize='64')


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
