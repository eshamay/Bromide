import numpy, sys, operator, PlotUtility, scipy.signal, Smoothing
from ColumnDataFile import ColumnDataFile as CDF
import matplotlib.pyplot as plt
import Smoothing
import AutoCorrelation
import operator
import glob

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
  return numpy.array(freqs*freqs*abs(fft)*abs(fft))

def NewCorr(x,y=None):
	ret = []
	if y == None:	# the autocorrelation case
		ret = numpy.array(numpy.correlate(x,x,"full"))[::-1]	# reverse the list with [::-1]
	else:	# cross-correlate x and y
		ret = numpy.array(numpy.correlate(x,y,"full"))[::-1]

	ret = ret[len(ret)/2:]	# only take the 1st half of the result (2nd half is just reversed)
	ret = [n / (len(ret)-tau) for n,tau in zip(ret,range(len(ret)))]	# do the ensemble average over the time lags
	return ret

	''' The original brute-force method -- pretty slow '''
def ManualCorrelate(func,tau,mean,x):
	Nstep = float(len(x))
	Nav = int(Nstep-tau)

	d = 0.0
	for n in range(Nav):
		d = d + func(x[n]-mean, x[tau+n]-mean)	# Autocorrelation
	return d/float(Nav)		# ensemble averaging


def vectorAutocorr(x,tau):
	ManualAutocorr(numpy.dot,x,tau)

def ScalarAutocorr(x,tau):
	ManualAutocorr(operator.mul, tau, x)

def TCFAxis(fig_num=1):
	# plot the tcf
	fig = plt.figure(num=fig_num, facecolor='w', edgecolor='w', frameon=True)
	axs = fig.add_subplot(1,1,1)
	axs.set_ylabel(r'ACF', size='xx-large')
	axs.set_xlabel(r'Timestep', size='xx-large')
	return axs

def PlotTCF(tcf,axs):
	#axs = TCFAxis()
	axs.plot(range(len(tcf)), tcf, linewidth=3.0, color='k')

def AverageTCFs(tcfs):
	avg_tcf = numpy.array(reduce(operator.add,tcfs))/len(tcfs)
	return avg_tcf

def LoadCDFs(files):
	# load the column-data files
	return [CDF(f) for f in files]

def ParseColumnFileData(cdfs,cols):
	data = []
	for cdf in cdfs:
		for col in range(cols):
			data.append (numpy.array(cdf[col]))
	return data

def PowerSpectrumAxis():
	# Set up the figure for the spectral plots
	fig = plt.figure(num=None, facecolor='w', edgecolor='w', frameon=True)
	#axs = fig.add_subplot(1,1,1, autoscaleon=False)
	axs = plt.subplot(1,1,1)
	axs.set_ylabel(r'Spectrum', size='xx-large')
	axs.set_xlabel(r'Frequency / cm$^{-1}$', size='xx-large')
	return axs

def PlotPowerSpectrum(tcf, axs):
	#axs = PowerSpectrumAxis()
	# use the same axis for all the spectra
	freq_axis = FreqAxis(len(tcf),dt)
	fft = FFT(tcf)
	smooth = Smoothing.window_smooth(fft, window='hanning', window_len=15)
	axs.plot(freq_axis, smooth, linewidth=3.5)


def PlotPowerSpectra(tcfs):
	axs = PowerSpectrumAxis()
	# use the same axis for all the spectra
	freq_axis = FreqAxis(len(tcfs[0]),dt)

	# fft each of the individual tcfs
	fft_tcfs = [FFT(t) for t in tcfs]

	# take an average of all the ffts
	#avg_fft = numpy.array(reduce(operator.add, fft_tcfs))/len(fft_tcfs)

	for f in fft_tcfs:
		axs.plot(freq_axis, f, linewidth=2.0)


def LoadColumnFileData(path, num_cols):
	files = glob.glob(path)
	cdfs = LoadCDFs(files)
	data = ParseColumnFileData(cdfs,num_cols)
	return data


def CalcAverageTCFFromFiles(path, num_cols):
	data = LoadColumnFileData(path, num_cols)
	tcfs = [numpy.array(AutoCorr1d(d)) for d in data]
	avg = AverageTCFs(tcfs)
	return avg




'''
#########********** Do the time-domain work here - calculate the TCFs of the data files ********############
#########                                                                                       ############

#bond_avg = CalcAverageTCFFromFiles('[1-5]/so2-bondlengths.dat', 2)
#angle_avg = CalcAverageTCFFromFiles('[1-5]/so2-angles.dat', 1)
#closest_water_oh_avg = CalcAverageTCFFromFiles('[1-5]/closest-water-bondlengths.dat', 6)
closest_water_oh_avg = CalcAverageTCFFromFiles('[1-5]/closest-water-bondlengths.dat', 1)
second_closest_water_oh_avg = CalcAverageTCFFromFiles('[1-5]/2nd-closest-water-bondlengths.dat', 1)
top_water_oh_avg = CalcAverageTCFFromFiles('[1-5]/top-water-bondlengths.dat', 1)
#closest_water_angle_avg = CalcAverageTCFFromFiles('[1-5]/h2o-angle.dat', 3)

axs = TCFAxis()
PlotTCF(closest_water_oh_avg, axs)
#########********** Do the frequency-domain work here - calculate the FFT of the TCFs **********############
#########                                                                                       ############

axs = PowerSpectrumAxis()
#PlotPowerSpectrum(bond_avg, axs)
#PlotPowerSpectrum(angle_avg, axs)
PlotPowerSpectrum(closest_water_oh_avg, axs)
PlotPowerSpectrum(second_closest_water_oh_avg, axs)
PlotPowerSpectrum(top_water_oh_avg, axs)
#PlotPowerSpectrum(closest_water_angle_avg, axs)

plt.xlim(0,6000)
#PlotUtility.ShowLegend(axs)
plt.show()
'''
