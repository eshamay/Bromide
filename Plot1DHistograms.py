import sys
import glob
from ColumnDataFile import ColumnDataFile as CDF
from scipy import interpolate
from scipy import signal
import matplotlib.pyplot as plt
from pylab import *
import numpy
import PlotUtility


def smooth(x,window_len=11,window='hanning'):
	"""smooth the data using a window with requested size.
    
	This method is based on the convolution of a scaled window with the signal.
	The signal is prepared by introducing reflected copies of the signal 
	(with the window size) in both ends so that transient parts are minimized
	in the begining and end part of the output signal.
    
	input:
		x: the input signal 
		window_len: the dimension of the smoothing window; should be an odd integer
		window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
		flat window will produce a moving average smoothing.
		
	output:
		the smoothed signal
        
	example:

		t=linspace(-2,2,0.1)
		x=sin(t)+randn(len(t))*0.1
		y=smooth(x)
    
	see also: 
    
		numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
		scipy.signal.lfilter
 
		TODO: the window parameter could be the window itself if an array instead of a string   
	"""

	if x.ndim != 1:
		raise ValueError, "smooth only accepts 1 dimension arrays."

	if x.size < window_len:
		raise ValueError, "Input vector needs to be bigger than window size."


	if window_len<3:
		return x


	if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
		raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


	s=numpy.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
	#print(len(s))
	if window == 'flat': #moving average
		w=numpy.ones(window_len,'d')
  if window == 'gaussian':
    w=scipy.signal.gaussian(window_len,5.0)
	else:
		w=eval('numpy.'+window+'(window_len)')

	y=numpy.convolve(w/w.sum(),s,mode='same')
	return y[window_len:-window_len+1]




files = glob.glob('dat/'+sys.argv[1]+'.*')
cdfs = [CDF(f) for f in files]

num = len(cdfs[0][0])
xi = numpy.array(cdfs[0][0])
yi = numpy.zeros(num)

for c in cdfs:
	yi = yi + c[1]

#yi = yi / yi.sum()
#histo,edges = numpy.histogram(xi, weights=yi, normed=True, bins=400)

fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
axs = fig.add_subplot(1,1,1)
axs.plot(xi,yi)
	
# cubic fitting
'''
tck = interpolate.splrep(xi,yi,s=0)
xnew = numpy.arange(min(xi),max(xi),0.01)
ynew = interpolate.splev(xnew,tck,der=0)
'''
#ynew = smooth(yi,window_len=25,window='hamming')
#axs.plot(xi,ynew,linewidth=4,color='k')

xticks(fontsize=48)
#plt.xlim(-1.0,1.0)

yticks([])

'''
name = cdfs[0].filename
if "theta" in name or "oh" in name:
	plt.xlabel (r'$\cos(\theta)$', fontsize=64)
elif "phi" in name:
	plt.xlabel (r'$\cos(\phi)$', fontsize=64)
'''


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
numpy.savetxt(sys.argv[1]+".1d.txt", (xi,ynew), fmt="%12.6e")
'''

#import csv
#writer = csv.writer(open(sys.argv[1]+".1d.txt", "wb"),delimiter=' ')
#writer.writerows(zip(xi,ynew))

plt.show()
