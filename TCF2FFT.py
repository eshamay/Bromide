import numpy, sys, operator, PlotUtility, scipy.signal, Smoothing
from ColumnDataFile import ColumnDataFile as CDF
import matplotlib.pyplot as plt
import Smoothing

### Some constants for the calculations ###
# speed of light in cm/s
c = 29979245800.0 
# planck's
k = 1.3806504e-23
hbar = 1.05457148e-34
# timestep of the simulation
#dt = 0.75e-15
dt = 0.75e-15
T = 300.0	# in kelvin
Beta = 1/(k*T)

# Given a list of discrete data points, returns the ACF using numpy's built-in
def autocorr(x):
	result = numpy.correlate(x, x, mode='full')
	return result[result.size/2:]

def fft_correlate(A,B):
	#return scipy.signal.fftconvolve(A,B[::-1,::-1,...],*args,**kwargs)
	result = scipy.signal.fftconvolve(numpy.array(A),numpy.array(B),mode='full')
	return result[result.size/2:]

def fft_autocorr(x):
	result = fft_correlate(x, x)
	return result

# Using the Wiener-Kintchine theorem we can obtain the autocorrelation via FFT
def WK_autocorr(x):
	s = numpy.fft.fft(x)
	result = numpy.array(numpy.real(numpy.fft.ifft(s*numpy.conjugate(s)))/numpy.var(x))
	return result

# Given a TCF/ACF list or array, this average each point (tau) by the number of elements going in to calculate it. i.e. 1/T * Sum(...)
def ensembleAverage(x):
	return map(lambda i, tau: float(i)/float(len(x)-tau), x, range(len(x)))

def simpleTCF(x):
	return map (lambda f: f*x[0], x)

def timeAxis(x):
	return [i*dt/1000.0 for i in range(len(x))]

# load the column-data file
cdf = CDF (sys.argv[1])
tcf = cdf[0]

# do the fft here

#########********** Do the time-domain work here - calculate the ACF/TCFs **********############
#########                                                                           ############
# best convolve method = WK_acorr w/o ensemble avg and taking the 'full' without returning half the tcf
#acorr = [WK_autocorr(cdf[i]) for i in range(len(cdf))]	# use numpy's autocorrelation function via fft in the WK theorem
#tcf = acorr
#tcf = [simpleTCF(cdf[i]) for i in range(len(cdf))]	# just multiply evevrything by the first value
#tcf = [ensembleAverage(i) for i in acorr]	# finish the autocorr with ensemble averaging


#fig = plt.figure(num=None, facecolor='w', edgecolor='w', frameon=True)
#axs = fig.add_subplot(1,1,1)
#axs.set_ylabel(r'TCF', size='xx-large')
#axs.set_xlabel(r'Time / ps', size='xx-large')
#map (lambda x: axs.plot(timeAxis(tcf[x]), tcf[x]), range(len(tcf)))	


# simplest FT
def fft(x):
	return numpy.fft.fft(x)

# simple FT with a windowing function to clean up the data
def windowFFT(x,window_fn):
	# window to clean up the fft
	window = map (operator.mul, x, window_fn((len(x))))
	return fft(window)

# norm-squared of the fft
def squareFFT(x):
	return [abs(i)*abs(i) for i in x]

# multiply by the prefactor
def prefactorFFT(x,freq):
	#return [numpy.tanh(Beta*f*hbar/2.0)*4.0*numpy.pi*f/3.0/hbar/c*i for (i,f) in zip(x,freq)]
	return [numpy.tanh(f)*f*i for (i,f) in zip(x,freq)]

def freqAxis():
	#return [float(i)/len(x)/dt/c/2.0 for i in range(len(x))]		# old way
	axis = numpy.array(numpy.fft.fftfreq(n=len(tcf), d=dt))/c
	return axis

#########********** Do the frequency-domain work here - calculate the FFT of the TCFs **********############
#########                                                                                       ############

# Set up the figure for the spectral plots
fig = plt.figure(num=None, facecolor='w', edgecolor='w', frameon=True)
axs = fig.add_subplot(1,1,1)
axs.set_ylabel(r'Spectrum', size='xx-large')
axs.set_xlabel(r'Frequency / cm$^{-1}$', size='xx-large')

# use the same axis for all the spectra
freq_axis = freqAxis()

#window_fft = [windowFFT(tcf[x], numpy.hamming) for x in range(len(tcf))]
fft = squareFFT(fft(tcf))
#window_fft = [prefactorFFT(x,freq_axis) for x in window_fft]
#fft = [prefactorFFT(x,freq_axis) for x in fft]

# plots out the spectrum (or power spectrum if using abs(x)**2)
#map (lambda x: axs.plot (freq_axis, numpy.abs(x)), fft)
#test_tcf = map (lambda x,y,z: numpy.sqrt(x**2+y**2+z**2), tcf[0], tcf[1], tcf[2])
#avg_fft = [numpy.abs(prefactorFFT(fft(x), freq_axis))**2 for x in tcf]
#avg_fft = reduce (operator.add, avg_fft)
#print len(avg_fft)

#axs.plot(freq_axis, avg_fft, color='k', linewidth=0.3, label="Power Spectrum")
#hamming = numpy.abs(prefactorFFT(windowFFT(test_tcf,numpy.hamming),freq_axis))
#axs.plot(freq_axis, hamming, linewidth=1.0, color='r', label="Hamming")

smooth_gaussian = Smoothing.window_smooth(numpy.array(fft),window_len=20,window='gaussian')
axs.plot(freq_axis, fft, linewidth=2.0, color='r', label="Raw")
axs.plot(freq_axis, smooth_gaussian, linewidth=2.0, color='k', label="Smooth")
axs.set_xlim(1000,4000)
#axs.plot(freq_axis, numpy.abs(windowFFT(test_tcf))
#plt.xlim(1000,4000)


PlotUtility.ShowLegend(axs)
plt.show()

