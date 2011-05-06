from PlotPowerSpectra import *
from ColumnDataFile import ColumnDataFile as CDF
import sys

import matplotlib.pyplot as plt
import PlotUtility
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
	#smooth_spectrum = spectrum
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
	

def PlotFiles(files, axs, lbl):
	cdfs = [CDF(i) for i in files]

	tfs_sym = [TimeFunction(i[0]) for i in cdfs]
	tcfs_sym = [tf.TCF() for tf in tfs_sym]

	tfs_antisym = [TimeFunction(i[1]) for i in cdfs]
	tcfs_antisym = [tf.TCF() for tf in tfs_antisym]
	
	tcf_sum = AverageTCF(tcfs_sym + tcfs_antisym)
	tcf_sym = AverageTCF(tcfs_sym)
	tcf_antisym = AverageTCF(tcfs_antisym)
	
	(sym_freqs,sym_spectrum) = SmoothSpectrum(tcf_sym)
	(antisym_freqs,antisym_spectrum) = SmoothSpectrum(tcf_antisym)
	(sum_freqs,sum_spectrum) = SmoothSpectrum(tcf_sum)
	
	#axs.plot (sym_freqs, sym_spectrum, linewidth=2.5)
	#axs.plot (antisym_freqs, antisym_spectrum, linewidth=2.5)
	axs.plot (sum_freqs, sum_spectrum, linewidth=2.5, label=lbl)


#filename='h2o-bondlengths.oh.normal_modes.dat'
filename='so2-bond+angles.normal_modes.dat'
#filename='h2o-bondlengths.normal_modes.dat'
files_cold = glob.glob('[1-5]/'+filename)

files_hot = glob.glob('[6-9]/'+filename)
files_hot = files_hot + glob.glob('10/'+filename)

axs = PowerSpectrumAxis()
axs.set_xlabel(r'Frequency / cm$^{-1}$', fontsize='64')
axs.set_ylabel('Power Spectrum', fontsize='64')

PlotFiles (files_cold, axs, 'Cold')
PlotFiles (files_hot, axs, 'Hot')

xticks(fontsize=40)
yticks([])
plt.xlim(2800,4000)

PlotUtility.ShowLegend(axs)
plt.show()
