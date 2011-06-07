from PlotPowerSpectra import *
from ColumnDataFile import ColumnDataFile as CDF
import glob

import matplotlib.pyplot as plt

c = 29979245800.0		# speed of light (cm/s)
dt = 0.75e-15	# length of time between each simulation data point
tau = 20000	# length of the correlation function

class TimeFunction:
	def __init__(self, data):
		self.data = numpy.array(data)
                self.mean = numpy.average(data)
                self.covariance = (self.data * self.data).sum()
                self.data = self.data - self.mean

		self.tcf = numpy.array(NewCorr(self.data)[:correlation_tau])/self.covariance
		#self.time = numpy.array(range(len(self.tcf)))/dt*1.0e15	# in fs

		(self.freqs,self.spectrum) = SmoothSpectrum(self.tcf)

	def TCF(self):
		return self.tcf

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
	#avg_tcf = AverageTCF(avg_tcfs)

	# calculate the spectra
	#a = SmoothSpectrum(avg_tcf)
	#axs.plot(a[0], a[1], linewidth=3.0)

	spectra = [SmoothSpectrum(t) for t in avg_tcfs]
	for s in spectra:
		axs.plot (s[0], s[1], linewidth=1.5)

files = glob.glob('oco[1-2].dat')

cdfs = [CDF(f) for f in files]

axs = TCFAxis(1)
axs.set_xlabel(r'Time / ps', fontsize='64')
axs.set_ylabel('Bondlength', fontsize='64')
for cdf in cdfs:
        axs.plot(numpy.array(range(len(cdf[1])))*0.75, cdf[1])

tcfs = [Correlate(i[1]) for i in cdfs]

axs = TCFAxis(2)
axs.set_xlabel(r'Time Lag / ps', fontsize='64')
axs.set_ylabel('TCF', fontsize='64')
for t in tcfs:
        axs.plot(numpy.array(range(len(t)))*0.75,t)

axs = PowerSpectrumAxis(3)
axs.set_xlabel(r'Frequency / cm$^{-1}$', fontsize='64')
axs.set_ylabel('Power Spectrum', fontsize='64')
for t in tcfs:
        freqs_cold,spectrum_cold,smooth_spectrum_cold = PowerSpectrum(t)
        axs.plot(freqs_cold,smooth_spectrum_cold)
        
#axs.set_xlim(2800,4000)

'''
filename = 'h2o-bondlengths.normal_modes.dat'
#filename='so2-bond+angles.dat'
files_cold = glob.glob('[1-5]/'+filename)

files_hot = glob.glob('[6-9]/'+filename)
files_hot = files_hot + glob.glob('10/'+filename)

#filename = 'so2-bond+angles.dat'
filename='h2o-bondlengths.normal_modes.z.dat'
cold = glob.glob('[1-5]/'+filename)
hot = glob.glob('[6-9]/'+filename)
#hot = hot + glob.glob('10/'+filename)

cdfs_cold = [CDF(f) for f in cold]
cdfs_hot = [CDF(f) for f in hot]
tcfs_cold = [Correlate(i[0])[:tau] for i in cdfs_cold]
tcfs_cold = tcfs_cold + [Correlate(i[1])[:tau] for i in cdfs_cold]
tcfs_hot = [Correlate(i[0])[:tau] for i in cdfs_hot]
tcfs_hot = tcfs_hot + [Correlate(i[1])[:tau] for i in cdfs_hot]
time = range(tau)

avg_tcf_cold = AverageTCF(tcfs_cold)
avg_tcf_hot = AverageTCF(tcfs_hot)

axs = TCFAxis(1)
axs.set_xlabel(r'Time / ps', fontsize='64')
axs.set_ylabel('Lag', fontsize='64')
axs.plot(time,avg_tcf_cold)
axs.plot(time,avg_tcf_hot)

freqs_cold,spectrum_cold,smooth_spectrum_cold = PowerSpectrum(avg_tcf_cold)
freqs_hot,spectrum_hot,smooth_spectrum_hot = PowerSpectrum(avg_tcf_hot)

axs = PowerSpectrumAxis(2)
axs.set_xlabel(r'Frequency / cm$^{-1}$', fontsize='64')
axs.set_ylabel('Power Spectrum', fontsize='64')
#xticks(fontsize=36)
#yticks(fontsize=28)
axs.plot(freqs_cold,smooth_spectrum_cold)
axs.plot(freqs_hot,smooth_spectrum_hot)
axs.set_xlim(2800,4000)
'''

plt.show()
