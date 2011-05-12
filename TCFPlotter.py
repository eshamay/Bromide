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
		self.data = numpy.array(data)
                self.mean = numpy.average(data)
                self.data = self.data - self.mean
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
	

def PlotFiles(files, axs, cols, lbl):
	cdfs = [CDF(i) for i in files]

        # each column gets its own time function
	tfs = [[TimeFunction(c[i]) for c in cdfs] for i in cols]

	# calculate the correlation of each time function
	tcfs = [[t[i].TCF() for t in tfs] for i in cols]
	
	# average the correlations
	avg_tcfs = [AverageTCF(t) for t in tcfs]
        avg_tcf = AverageTCF(avg_tcfs)

	# calculate the spectra
        a = SmoothSpectrum(avg_tcf)
        axs.plot(a[0], a[1], linewidth=3.0)

	spectra = [SmoothSpectrum(t) for t in avg_tcfs]
	for s in spectra:
		axs.plot (s[0], s[1], linewidth=1.5)

filename = 'h2o-bondlengths.normal_modes.dat'
#filename='so2-bond+angles.dat'
files_cold = glob.glob('[1-5]/'+filename)

files_hot = glob.glob('[6-9]/'+filename)
files_hot = files_hot + glob.glob('10/'+filename)

axs = PowerSpectrumAxis()
axs.set_xlabel(r'Frequency / cm$^{-1}$', fontsize='64')
axs.set_ylabel('Power Spectrum', fontsize='64')

PlotFiles (files_cold, axs, [0,1], 'Cold')
PlotFiles (files_cold, axs, [0,1], 'Hot')

xticks(fontsize=40)
yticks([])
plt.xlim(2800,4000)

#PlotUtility.ShowLegend(axs)
plt.show()
