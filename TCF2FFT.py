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
T = 250.0	# in kelvin
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

def vectorAutocorr(x,tau):
  Nstep = float(len(x))
  Ncorr = tau
  Nshift = 1
  Nav = int((Nstep-Ncorr)/Nshift)

  d = 0.0
  for n in range(Nav):
    d = d + numpy.dot(x[n], x[Ncorr+n])
  return d/float(Nav)


def d_fft_n_range(n):
  return range(-n,n,1)[1:]

def d_wk(k,N,dtau):
  return k*2*numpy.pi/(2*(N-1)*dtau)

def alpha_k(k,N,C,dtau):
  w_k = d_wk(k,N,dtau)
  w_k = w_k * w_k

  sum = 0.0
  for n in d_fft_n_range(N):
    ex = numpy.exp(-2.0*numpy.pi*(1.0j)*k*n/(2*(N-1))) * C[n] * dtau
    sum = sum + ex

  return sum * w_k

def alpha_fft(C,dtau):
  return [alpha_k(k,len(C),C,dtau) for k in range(len(C))]

# simplest FT
def fft(x):
	return numpy.fft.rfft(x)

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

#########********** Do the time-domain work here - calculate the ACF/TCFs **********############
#########                                                                           ############

# load the column-data file
cdf = CDF (sys.argv[1])
dipoles = [numpy.array([cdf[0][i], cdf[1][i], cdf[2][i]]) for i in range(len(cdf[0]))]
tcf = [vectorAutocorr(dipoles,tau) for tau in range(3500)]

# plot the tcf
fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
axs = fig.add_subplot(1,1,1)
axs.set_ylabel(r'ACF', size='xx-large')
axs.set_xlabel(r'Timestep', size='xx-large')
axs.plot (range(len(tcf)), tcf)


# apply a windowing function
tcf = [t * hwin for t,hwin in zip(tcf,numpy.hanning(len(tcf)))]


#########********** Do the frequency-domain work here - calculate the FFT of the TCFs **********############
#########                                                                                       ############

# Set up the figure for the spectral plots
fig = plt.figure(num=None, facecolor='w', edgecolor='w', frameon=True)
axs = fig.add_subplot(1,1,1)
axs.set_ylabel(r'Spectrum', size='xx-large')
axs.set_xlabel(r'Frequency / cm$^{-1}$', size='xx-large')

# use the same axis for all the spectra
freq_axis = numpy.array(numpy.fft.fftfreq(n=len(tcf), d=dt))/c
#freq_axis = freq_axis[:len(freq_axis)/2+1]
print len(freq_axis)

#fft_tcf = fft(tcf)
#fft_tcf = alpha_fft(tcf,dt)
fft_tcf = numpy.fft.fft(tcf)
print len(fft_tcf)
#fft_tcf = fft_tcf[len(fft_tcf)/2:]
#fft_tcf = [abs(f)*abs(f) for f in fft_tcf]
fft_tcf = [w*w*abs(f)*abs(f) for w,f in zip(freq_axis,fft_tcf)]

axs.plot(freq_axis, fft_tcf, linewidth=2.0, color='r')
plt.xlim(1000,4000)

#PlotUtility.ShowLegend(axs)
plt.show()

