from PlotPowerSpectra import *
from ColumnDataFile import ColumnDataFile as CDF
import sys

import matplotlib.pyplot as plt
from pylab import *

c = 29979245800.0		# speed of light (cm/s)
dt = 0.75e-15	# length of time between each simulation data point
correlation_tau = 7000	# length of the correlation function

def WindowFunction(data):
	return numpy.hanning(len(data))

def SmoothFunction(data):
	return WindowFunction(data) * data

def SmoothSpectrum(data):
	smooth_data = SmoothFunction(data)
	freqs = numpy.array(numpy.fft.fftfreq(n=len(data), d=dt))/c
	fft = numpy.array(numpy.fft.fft(smooth_data))
	fft = fft * freqs
	spectrum = abs(fft) * abs(fft)
	smooth_spectrum = Smoothing.window_smooth(spectrum, window_len=5)
	return (freqs,smooth_spectrum)

class TimeFunction:
	def __init__(self, data):
		self.data = data

		self.tcf = numpy.array(NewCorr(self.data)[:correlation_tau])
		self.time = numpy.array(range(len(self.tcf)))/dt*1.0e15	# in fs
		(self.freqs,self.spectrum) = SmoothSpectrum(self.tcf)

	def TCF(self):
		return self.tcf

	def Time(self):
		return self.time

	def Freqs(self):
		return self.freqs

	def Spectrum(self):
		return self.spectrum


def AverageTCF (tcfs):
	sum_tcf = numpy.array(reduce(operator.add, tcfs))
	return sum_tcf / len(tcfs)
	

files = glob.glob('[6-9]/h2o-bondlengths.normal_modes.dat')
files2 = glob.glob('10/h2o-bondlengths.normal_modes.dat')
files = files + files2
cdfs = [CDF(i) for i in files]
tfs_sym = [TimeFunction(i[0]) for i in cdfs]
tcfs_sym = [tf.TCF() for tf in tfs_sym]
tfs_antisym = [TimeFunction(i[1]) for i in cdfs]
tcfs_antisym = [tf.TCF() for tf in tfs_antisym]

tcf_sym = AverageTCF(tcfs_sym)
tcf_antisym = AverageTCF(tcfs_antisym)

(sym_freqs,sym_spectrum) = SmoothSpectrum(tcf_sym)
(antisym_freqs,antisym_spectrum) = SmoothSpectrum(tcf_antisym)

axs = PowerSpectrumAxis()
axs.plot (sym_freqs[800:], sym_spectrum[800:], linewidth=2.5)
axs.plot (antisym_freqs[800:], antisym_spectrum[800:], linewidth=2.5)

axs.set_xlabel(r'Frequency / cm$^{-1}$', fontsize='64')
axs.set_ylabel('Power Spectrum', fontsize='64')
xticks(fontsize=40)
yticks([])
plt.xlim(2800,4000)

def Autoscale(ax,data):
	#start,stop = plt.xlim()
	#print start, stop
	#d = data[start:stop+1]
	ax.dataLim._points[:,1] = 0.0, max(data)
	ax.autoscale_view(scalex=False,scaley=True)

#Autoscale(axs,antisym_spectrum)
#Autoscale(axs,sym_spectrum)
#print axs.dataLim._points

plt.show()

