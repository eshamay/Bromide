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

def ManualCorrelate(func,tau,x,y=None):
	Nstep = float(len(x))
	Ncorr = tau
	Nshift = 1
	Nav = int((Nstep-Ncorr)/Nshift)

	d = 0.0
	for n in range(Nav):
		if y == None:
			d = d + func(x[n], x[Ncorr+n])	# Autocorrelation
		else:
			d = d + func(x[n], y[Ncorr+n])	# Cross correlation
	return numpy.array(d)/float(Nav)

def vectorAutocorr(x,tau):
	ManualAutocorr(numpy.dot,x,tau)

def ScalarAutocorr(x,tau):
	ManualAutocorr(operator.mul, tau, x)

def TCFAxis():
	# plot the tcf
	fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
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
	axs = fig.add_subplot(1,1,1)
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
