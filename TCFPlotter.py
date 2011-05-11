from PlotPowerSpectra import *
from ColumnDataFile import ColumnDataFile as CDF
import sys
import scipy.stats

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
	#smooth_spectrum = Smoothing.window_smooth(spectrum, window_len=5)
	smooth_spectrum = spectrum
	return (freqs,smooth_spectrum)

class TimeFunction:
	def __init__(self, data):
		self.data = data
		self.tcf = numpy.array(NewCorr(self.data))
		self.time = numpy.array(range(len(self.data)))/dt*1.0e15	# in fs
		self.tcftime = numpy.array(range(len(self.tcf)))/dt*1.0e15	# in fs
		(self.freqs,self.spectrum) = SmoothSpectrum(self.tcf)

	def Data(self):
	  	return self.data

	def TCF(self):
		return self.tcf

	def Time(self):
		return self.time

	def TCFTime(self):
		return self.tcftime

	def Freqs(self):
		return self.freqs

	def Spectrum(self):
		return self.spectrum


def AverageTCF (tcfs):
	sum_tcf = numpy.array(reduce(operator.add, tcfs))
	return sum_tcf / len(tcfs)
	

def PlotFiles(files, axs, cols, labels):
	cdfs = [CDF(i) for i in files]

	# each column gets its own time function
	tfs = [[TimeFunction(cdfs[c][i]) for c in range(len(cdfs))] for i in cols]

	# calculate the correlation of each time function
	tcfs = [[t[i].TCF() for t in tfs] for i in cols]

	# average the correlations
	avg_tcfs = [AverageTCF(t) for t in tcfs]

	# calculate the spectra
	spectra = [SmoothSpectrum(t) for t in avg_tcfs]

	for s,l in zip(spectra,labels):
		axs.plot (s[0], s[1], linewidth=2.5, label=l)
	#axs.plot (antisym_freqs, antisym_spectrum, linewidth=2.5, label='antisym')
	#axs.plot (sum_freqs, sum_spectrum, linewidth=2.5, label=lbl)


#filename='h2o-bondlengths.oh.normal_modes.dat'
filename='so2-bond+angles.normal_modes.dat'
#filename='h2o-bondlengths.normal_modes.dat'
files_cold = glob.glob('[1-5]/'+filename)

#cdf = CDF(sys.argv[1])

# start with a spectrum and work backwards
# first, let's create a simple gaussian distribution of frequencies
gaussian = lambda x: 3*exp(-(1050-x)**2/5.5)
#dw = 1
spectrum = [gaussian(i) for i in range(2000)]

sqrt_spectrum = numpy.array([sqrt(s) for s in spectrum])
sqrt_spectrum = [s/(w+1)*1.0j for s,w in zip(sqrt_spectrum,range(2000))]
#files_hot = glob.glob('[6-9]/'+filename)
#files_hot = files_hot + glob.glob('10/'+filename)
#axs = TCFAxis(1)
#axs.set_xlabel(r'Time / ps', fontsize='64')
#axs.set_ylabel('Bondlength', fontsize='64')
#axs.plot(time, tf)

#PlotFiles (glob.glob(filename), axs, 'Temp')
#PlotFiles (files_cold, axs, [0], ['oh'])
#PlotFiles (files_hot, axs, 'Hot')

#xticks(fontsize=40)
#yticks([])
#plt.xlim(2800,4000)

#PlotUtility.ShowLegend(axs)
plt.show()
