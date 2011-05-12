import numpy, sys, operator, PlotUtility, scipy.signal, Smoothing
from ColumnDataFile import ColumnDataFile as CDF
import matplotlib.pyplot as plt
import Smoothing
import AutoCorrelation
import operator

### Some constants for the calculations ###
# speed of light in cm/s
c = 29979245800.0 
# planck's
k = 1.3806504e-23
hbar = 1.05457148e-34
# timestep of the simulation
dt = 0.75e-15
correlation_tau = 8000
T = 300.0	# in kelvin
Beta = 1/(k*T)

#########********** Do the time-domain work here - calculate the ACF/TCFs **********############
#########                                                                           ############

def WindowFunction(data):
	return numpy.hanning(len(data))

def SmoothFunction(data):
	return WindowFunction(data) * data

def FreqAxis(N,dt):
  return numpy.array(numpy.fft.fftfreq(n=N, d=dt))/c

def FFT(tcf):
  tcf = numpy.array(tcf)
  N = len(tcf)
  window = numpy.hanning(N)
  tcf = tcf*window
  fft = numpy.fft.fft(tcf)
  freqs = FreqAxis(len(tcf),dt)
  return numpy.array(freqs*freqs*abs(fft)*abs(fft))

''' compute either the cross-correlation or the autocorrelation of time functions '''
''' equivalent to the methods used in the NIST engineering handbook:
		http://itl.nist.gov/div898/handbook/eda/section3/eda35c.htm
'''
def Correlate(x,y=None):
	ret = []
	x = x - numpy.average(x)
	if y == None:	# the autocorrelation case
		ret = numpy.array(numpy.correlate(x,x,"full"))[::-1]	# reverse the list with [::-1]
	else:	# cross-correlate x and y
		y = y - numpy.average(y)
		ret = numpy.array(numpy.correlate(x,y,"full"))[::-1]

	ret = ret[len(ret)/2:]	# only take the 1st half of the result (2nd half is just reversed)
	#ret = [n / (len(ret)-tau) for n,tau in zip(ret,range(len(ret)))]	# do the ensemble average over the time lags
	return ret/ret[0]


''' 	Brute-force direct calculation from the nist reference above
def autocovariance (x,lag,mean_x,y=None,mean_y=None):
	x = numpy.array(x) - mean_x

	sum = 0.0
  	if y == None:
		for i in range(len(x)-lag):
	  		sum = sum + x[i]*x[i+lag]
	else:
		y = numpy.array(y) - mean_y
		for i in range(len(x)-lag):
		  	sum = sum + x[i]*y[i+lag]
	
	return sum/(len(x)-lag)
'''


<<<<<<< HEAD
'''
def covariance(x,y=None):
  	ret = 0.0
	if y == None:
	  	ret = numpy.std(x) * numpy.std(x)
	else:
	  	ret = numpy.std(x) * numpy.std(y)
'''

'''
	ret = 0.0
	x = numpy.array(x) - mean_x
	sum_x = 0.0
  	for i in range(len(x)):
		sum_x = sum_x + x[i]*x[i]
	ret = sum_x

	sum_y = 0.0
	if y != None:
	  	y = numpy.array(y) - mean_y
		for i in range(len(y)):
		  	sum_y = sum_y + y[i]*y[i]
		ret = numpy.sqrt(sum_x*sum_y)
	return ret
'''
''' The original brute-force method -- pretty slow '''
'''
def ManualCorrelate(func,tau,x,y=None):
=======
	''' The original brute-force method -- pretty slow '''
def ManualCorrelate(func,tau,mean,x):
>>>>>>> e8e27416c9c4723ecf9faf9aab914b5f4043279e
	Nstep = float(len(x))
	Nav = int(Nstep-tau)

	d = 0.0
	for n in range(Nav):
<<<<<<< HEAD
		if y == None:
			d = d + func(x[n], x[Ncorr+n])	# Autocorrelation
		else:
			d = d + func(x[n], y[Ncorr+n])	# Cross correlation
	return numpy.array(d)/float(Nav)		# ensemble averaging
'''
=======
		d = d + func(x[n]-mean, x[tau+n]-mean)	# Autocorrelation
	return d/float(Nav)		# ensemble averaging
>>>>>>> e8e27416c9c4723ecf9faf9aab914b5f4043279e


def TCFAxis(fig_num=1):
	# plot the tcf
	fig = plt.figure(num=fig_num, facecolor='w', edgecolor='w', frameon=True)
	axs = fig.add_subplot(1,1,1)
	axs.set_ylabel(r'ACF', size='xx-large')
	axs.set_xlabel(r'Timestep', size='xx-large')
	return axs

def AverageTCFs(tcfs):
	avg_tcf = numpy.array(reduce(operator.add,tcfs))/len(tcfs)
	return avg_tcf

def PowerSpectrumAxis(n):
	# Set up the figure for the spectral plots
	fig = plt.figure(num=n, facecolor='w', edgecolor='w', frameon=True)
	#axs = fig.add_subplot(1,1,1, autoscaleon=False)
	axs = plt.subplot(1,1,1)
	axs.set_ylabel(r'Spectrum', size='xx-large')
	axs.set_xlabel(r'Frequency / cm$^{-1}$', size='xx-large')
	return axs

def PowerSpectrum(data):
	# power spectrum
	smooth_data = SmoothFunction(data)
	freqs = numpy.array(numpy.fft.fftfreq(n=len(data), d=dt))/c
	fft = numpy.array(numpy.fft.fft(smooth_data))
	fft = fft * freqs
	spectrum = abs(fft) * abs(fft)
	smooth_spectrum = Smoothing.window_smooth(spectrum, window_len=5)
  	return (freqs,spectrum,smooth_spectrum)
