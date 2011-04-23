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
#k = 1.3806504e-23
#hbar = 1.05457148e-34
# timestep of the simulation
#dt = 0.75e-15
dt = 0.75e-15
correlation_tau = 8000
#T = 340.0	# in kelvin
#Beta = 1/(k*T)

#########********** Do the time-domain work here - calculate the ACF/TCFs **********############
#########                                                                           ############

def AutoCorr1d(x):
  return AutoCorrelation.fftCyclicAutocorrelation1D(x)

def AverageTCF(tcf):
  N = len(tcf)
  tcf = [tcf[i] / (N-i) for i in range(len(tcf))]
  return numpy.array(tcf)

def FreqAxis(N,dt):
  return numpy.array(numpy.fft.fftfreq(n=N, d=dt))/c

def FFT(tcf):
  tcf = numpy.array(tcf)
  N = len(tcf)
  window = numpy.hanning(N)
  tcf = tcf*window
  fft = numpy.fft.fft(tcf)
  freqs = FreqAxis(len(tcf),dt)
  return freqs*freqs*abs(fft)*abs(fft)

def vectorAutocorr(x,tau):
  Nstep = float(len(x))
  Ncorr = tau
  Nshift = 1
  Nav = int((Nstep-Ncorr)/Nshift)

  d = 0.0
  for n in range(Nav):
    d = d + numpy.dot(x[n], x[Ncorr+n])
  return numpy.array(d)/float(Nav)


# load the column-data files
cdfs = [CDF(i) for i in sys.argv[1:]]

# grab the different components of the dipole moment
data_x = [i[0] for i in cdfs]
data_y = [i[1] for i in cdfs]
data_z = [i[2] for i in cdfs]

# get the dipoles teased out of the data files
dipoles = [[numpy.array([x,y,z]) for x,y,z in zip(d_x,d_y,d_z)] for d_x,d_y,d_z in zip(data_x,data_y,data_z)]

# perform the autocorrelation on each data point to get the time correlation functions
tcfs_x = [numpy.array(AutoCorr1d(d)) for d in data_x]
tcfs_y = [numpy.array(AutoCorr1d(d)) for d in data_y]
tcfs_z = [numpy.array(AutoCorr1d(d)) for d in data_z]

# perform the autocorrelation using a vector inner-product
vector_tcfs = [[vectorAutocorr(d,tau) for tau in range(correlation_tau)] for d in dipoles]

# average the time correlation functions
tcfs = [(x+y+z)/3.0 for x,y,z in zip(tcfs_x,tcfs_y,tcfs_z)]
#tcfs = tcfs_z
avg_tcf = numpy.array(reduce(operator.add,tcfs))/len(tcfs)
avg_tcf = avg_tcf[:len(avg_tcf)/2]

avg_vector_tcf = numpy.array(reduce(operator.add, vector_tcfs))/len(vector_tcfs)
#avg_vector_tcf = avg_vector_tcf[:len(avg_tcf)/2]

# plot the tcf
fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
axs = fig.add_subplot(1,1,1)
axs.set_ylabel(r'ACF', size='xx-large')
axs.set_xlabel(r'Timestep', size='xx-large')

#for t in tcfs:
	#axs.plot(range(len(t)), t)
#axs.plot(range(len(avg_tcf)), avg_tcf, linewidth=3.0, color='k')
axs.plot(range(len(avg_tcf)), avg_tcf, linewidth=3.0, color='k')
axs.plot(range(len(avg_vector_tcf)), avg_vector_tcf, linewidth=3.0, color='k')

#########********** Do the frequency-domain work here - calculate the FFT of the TCFs **********############
#########                                                                                       ############

# Set up the figure for the spectral plots
fig = plt.figure(num=None, facecolor='w', edgecolor='w', frameon=True)
axs = fig.add_subplot(1,1,1)
axs.set_ylabel(r'Spectrum', size='xx-large')
axs.set_xlabel(r'Frequency / cm$^{-1}$', size='xx-large')

# use the same axis for all the spectra
freq_axis = FreqAxis(len(tcfs[0])/2,dt)

fft_tcfs_x = [FFT(t) for t in tcfs_x]
fft_tcfs_y = [FFT(t) for t in tcfs_y]
fft_tcfs_z = [FFT(t) for t in tcfs_z]

fft_tcfs = [(x+y+z)/3.0 for x,y,z in zip(fft_tcfs_x,fft_tcfs_y,fft_tcfs_z)]
#avg_fft = reduce(operator.add, fft_tcfs_z)/len(fft_tcfs_z)
avg_fft = FFT(avg_tcf)
smooth = Smoothing.window_smooth(avg_fft)


#for f in fft_tcfs:
	#axs.plot(freq_axis, f, linewidth=2.0)
axs.plot(freq_axis, smooth, linewidth=4.0, linestyle='-', color='k')

vector_tcf_ffts = [FFT(t) for t in vector_tcfs]
avg_vector_fft = numpy.array(reduce(operator.add,vector_tcf_ffts))/len(vector_tcf_ffts)

smooth_vector_fft = Smoothing.window_smooth(avg_vector_fft)
vector_freq_axis = FreqAxis(correlation_tau,dt)
print len(vector_freq_axis)
print len(smooth_vector_fft)
axs.plot(vector_freq_axis, smooth_vector_fft, linewidth=4.0, linestyle='-', color='r')

plt.xlim(0,6000)
#PlotUtility.ShowLegend(axs)
plt.show()
