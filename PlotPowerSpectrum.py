import numpy, sys, operator, PlotUtility, scipy.signal, Smoothing
from ColumnDataFile import ColumnDataFile as CDF
import matplotlib.pyplot as plt
import Smoothing
import AutoCorrelations

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

def freqAxis():
	#return [float(i)/len(x)/dt/c/2.0 for i in range(len(x))]		# old way
	axis = numpy.array(numpy.fft.fftfreq(n=len(tcf), d=dt))/c
	return axis

#########********** Do the time-domain work here - calculate the ACF/TCFs **********############
#########                                                                           ############

# load the column-data file
cdfs = [CDF(i) for i in sys.argv[1:]]
data = [numpy.array(i[1]) for i in cdfs]
tcfs = [AutoCorrelations.fftCyclicAutocorrelation1D(d) for d in data]

# plot the tcf
fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
axs = fig.add_subplot(1,1,1)
axs.set_ylabel(r'ACF', size='xx-large')
axs.set_xlabel(r'Timestep', size='xx-large')

for t in tcfs:
	axs.plot(range(len(t)), t)

N = len(tcfs[0])
window = numpy.blackman(N)
tcfs = [t*window for t in tcfs]

#########********** Do the frequency-domain work here - calculate the FFT of the TCFs **********############
#########                                                                                       ############

# Set up the figure for the spectral plots
fig = plt.figure(num=None, facecolor='w', edgecolor='w', frameon=True)
axs = fig.add_subplot(1,1,1)
axs.set_ylabel(r'Spectrum', size='xx-large')
axs.set_xlabel(r'Frequency / cm$^{-1}$', size='xx-large')

# use the same axis for all the spectra
freq_axis = numpy.array(numpy.fft.fftfreq(n=N, d=dt))/c
#freq_axis = freq_axis[:len(freq_axis)/2+1]

fft_tcfs = [numpy.fft.fft(t) for t in tcfs]
fft_tcfs = [freq_axis*freq_axis*abs(f)*abs(f) for f in fft_tcfs]

for f in fft_tcfs:
	axs.plot(freq_axis, f, linewidth=2.0)

#fft_smooth = Smoothing.window_smooth(fft_tcf,window_len=10,window='gaussian')
#axs.plot(freq_axis, fft_smooth, linewidth=2.0, color='b')

plt.xlim(0,6000)
#PlotUtility.ShowLegend(axs)
plt.show()

