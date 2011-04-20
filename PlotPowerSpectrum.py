import numpy, sys, operator, PlotUtility, scipy.signal, Smoothing
from ColumnDataFile import ColumnDataFile as CDF
import matplotlib.pyplot as plt
import Smoothing
import AutoCorrelation
from glob import glob
import operator

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

#########********** Do the time-domain work here - calculate the ACF/TCFs **********############
#########                                                                           ############

def AutoCorr1d(x):
  return AutoCorrelation.fftCyclicAutocorrelation1D(x)

def FreqAxis(x,dt):
  return numpy.array(numpy.fft.fftfreq(n=len(x), d=dt))/c

def FFT(x):
  x = numpy.array(x)
  N = len(x)
  window = numpy.hanning(N)
  x = x*window
  fft = numpy.fft.fft(x)
  freqs = FreqAxis(x,dt)
  return freqs*freqs*abs(fft)*abs(fft)

# load the column-data file
cdfs = [CDF(i) for i in sys.argv[1:]]
data_x = [i[0] for i in cdfs]   # grab the z-component
data_y = [i[1] for i in cdfs]   # grab the z-component
data_z = [i[2] for i in cdfs]   # grab the z-component
tcfs_x = [numpy.array(AutoCorr1d(d)) for d in data_x]
tcfs_y = [numpy.array(AutoCorr1d(d)) for d in data_y]
tcfs_z = [numpy.array(AutoCorr1d(d)) for d in data_z]

tcfs = [(x+y+z)/3.0 for x,y,z in zip(tcfs_x,tcfs_y,tcfs_z)]
avg_tcf = reduce(operator.add, tcfs)
avg_tcf = avg_tcf[:len(avg_tcf)/2]

# plot the tcf
fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
axs = fig.add_subplot(1,1,1)
axs.set_ylabel(r'ACF', size='xx-large')
axs.set_xlabel(r'Timestep', size='xx-large')

#for t in tcfs:
	#axs.plot(range(len(t)), t)
axs.plot(range(len(avg_tcf)), avg_tcf, linewidth=3.0, color='k')

#########********** Do the frequency-domain work here - calculate the FFT of the TCFs **********############
#########                                                                                       ############

# Set up the figure for the spectral plots
fig = plt.figure(num=None, facecolor='w', edgecolor='w', frameon=True)
axs = fig.add_subplot(1,1,1)
axs.set_ylabel(r'Spectrum', size='xx-large')
axs.set_xlabel(r'Frequency / cm$^{-1}$', size='xx-large')

# use the same axis for all the spectra
freq_axis = FreqAxis(tcfs[0][:len(tcfs[0])/2],dt)

fft_tcfs_x = [FFT(t[:len(t)/2]) for t in tcfs_x]
fft_tcfs_y = [FFT(t[:len(t)/2]) for t in tcfs_y]
fft_tcfs_z = [FFT(t[:len(t)/2]) for t in tcfs_z]
fft_tcfs = [(x+y+z)/3.0 for x,y,z in zip(fft_tcfs_x,fft_tcfs_y,fft_tcfs_z)]
avg_fft = reduce(operator.add, fft_tcfs)/len(fft_tcfs)
#smooth = Smoothing.window_smooth(avg_fft)

#for f in fft_tcfs:
	#axs.plot(freq_axis, f, linewidth=2.0)
axs.plot(freq_axis, avg_fft, linewidth=4.0, linestyle='-', color='k')

plt.xlim(0,6000)
#PlotUtility.ShowLegend(axs)
plt.show()

