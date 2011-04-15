import numpy, sys, operator, PlotUtility, scipy.signal, Smoothing
from ColumnDataFile import ColumnDataFile as CDF
import matplotlib.pyplot as plt
import Smoothing
import Autocorr

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

def vectorAutocorr(x,tau):
  Nstep = float(len(x))
  Ncorr = tau
  Nshift = 1
  Nav = int((Nstep-Ncorr)/Nshift)

  d = 0.0
  for n in range(Nav):
    d = d + numpy.dot(x[n], x[Ncorr+n])
  return d/float(Nav)

def autocorr(x,tau):
  Nstep = float(len(x))
  Ncorr = tau
  Nshift = 1
  Nav = int((Nstep-Ncorr)/Nshift)

  d = 0.0
  for n in range(Nav):
    d = d + x[n]*x[Ncorr+n]
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

# simple FT with a windowing function to clean up the data
def windowFFT(x,window_fn):
	# window to clean up the fft
	window = map (operator.mul, x, window_fn((len(x))))
	return fft(window)


#########********** Do the time-domain work here - calculate the ACF/TCFs **********############
#########                                                                           ############

def numpyTCF(data):
  result = numpy.array(numpy.correlate(data, data, mode='same'))
  result = result[result.size/2:]
  return result

def webTCF(data):
  return Autocorr.fftCyclicAutocorrelation1D(data)

def TCF(vals):
  tcf = [autocorr(vals,tau) for tau in range(3000)]
  return tcf

# plot the tcf
def PlotTCF(tcf,axs):
  axs.plot (range(len(tcf)), tcf)




#########********** Do the frequency-domain work here - calculate the FFT of the TCFs **********############
#########                                                                                       ############
def WindowTCF(tcf):
  tcf = [t * hwin for t,hwin in zip(tcf,numpy.hanning(len(tcf)))]
  return tcf

def PlotFFT(data,axs):
  data = WindowTCF(data)
  freq_axis = numpy.array(numpy.fft.fftfreq(n=len(data), d=dt))/c
  fft = numpy.fft.fft(data)
  fft = numpy.array([abs(f)*abs(f) for w,f in zip(freq_axis,fft)])
  N = len(freq_axis)
  axs.plot(freq_axis[10:N/2], fft[10:N/2], linewidth=2.0)

def PlotTCFFFT(tcf, axs):
  tcf = WindowTCF(tcf)
  # use the same axis for all the spectra
  freq_axis = numpy.array(numpy.fft.fftfreq(n=len(tcf), d=dt))/c
  #freq_axis = freq_axis[:len(freq_axis)/2+1]

  fft_tcf = numpy.fft.fft(tcf)
  #fft_tcf = fft_tcf[len(fft_tcf)/2:]
  #fft_tcf = [abs(f)*abs(f) for f in fft_tcf]
  fft_tcf = numpy.array([abs(f)*abs(f) for f in fft_tcf])
  #fft_tcf = numpy.array([w*w*abs(f)*abs(f) for w,f in zip(freq_axis,fft_tcf)])
  N = len(freq_axis)
  axs.plot(freq_axis[10:N/2+1], fft_tcf[10:N/2+1], linewidth=2.0)

  #fft_smooth = Smoothing.window_smooth(fft_tcf,window_len=10,window='gaussian')
  #axs.plot(freq_axis[10:N/2], fft_smooth[10:N/2], linewidth=2.0)



data = [numpy.array(CDF(file)[1]) for file in sys.argv[1:]]
tcf = [webTCF(d) for d in data]

# set up the tcf graph
fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
axs = fig.add_subplot(1,1,1)
axs.set_ylabel(r'ACF', size='xx-large')
axs.set_xlabel(r'Timestep', size='xx-large')
for t in tcf:
  PlotTCF(t,axs)

# Set up the figure for the spectral plots
fig = plt.figure(num=None, facecolor='w', edgecolor='w', frameon=True)
axs = fig.add_subplot(1,1,1)
axs.set_ylabel(r'Spectrum', size='xx-large')
axs.set_xlabel(r'Frequency / cm$^{-1}$', size='xx-large')

for t in tcf:
#for d in data:
  #PlotFFT(d,axs)
  PlotTCFFFT(t,axs)


plt.xlim(0,6000)
#PlotUtility.ShowLegend(axs)
plt.show()

