import sys
import glob
from ColumnDataFile import ColumnDataFile as CDF
import matplotlib.pyplot as plt
from pylab import *
import numpy
import PlotUtility
#from mpl_toolkits.axes_grid import make_axes_locatable
#import  matplotlib.axes as maxes


#files = glob.glob(sys.argv[1])
files = glob.glob('dat/'+sys.argv[1]+'*.dat')
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
	#new_data = new_data / new_data.max(0)	# old normalization scheme - doesn't really work or make sense
	zi = zi + new_data

# normalize to max = 1, min = 0
#zi = zi / zi.max(0)

fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
im = plt.imshow(zi, extent=(xmin,xmax,ymax,ymin), interpolation='bilinear', figure=fig, aspect='auto')


#divider = make_axes_locatable(ax)
#cax = divider.new_horizontal("5%", pad=0.05, axes_class=maxes.Axes)
#fig.add_axes(cax) 
#cbar = fig.colorbar(im, ticks=[0, .5, 1.0])
#cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='vertical',ticks=[0.,0.5,1.0],shrink=0.5)


xticks(fontsize=48)
plt.xlim(-10.0,10.0)
plt.xlabel (r'Distance to water slab surface / $\AA$', fontsize=64)


name = cdfs[0].filename
if "theta" in name:
	plt.ylabel (r'$\cos(\theta)$', fontsize=64)
	plt.ylim(-1.0,1.0)
	yticks([-1.0,-0.5,0.0,0.5,1.0], fontsize=48)
elif "phi" in name:
	plt.ylabel (r'$\cos(\phi)$', fontsize=64)
	plt.ylim(0.0,1.0)
	yticks([0.0,0.25,0.5,0.75,1.0], fontsize=48)

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
numpy.savetxt(sys.argv[1]+".1d.txt", data_1d, fmt="%12.6e")
'''

plt.show()
