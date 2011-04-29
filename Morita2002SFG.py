#!/usr/bin/python

# This SFG calculator creates both IR and SFG spectra from dipole/polarizability files
# The SFG spectra are all SSP polarized, but that can be changed in the DipPol analyzer

# PrintData prints out 5 columns for frequency and spectral intensity:
#		frequency, IR, SFG_X, SFG_Y, SFG_Z

# PlotData creates 2 figures, one for IR and one for SFG. The SFG figure has 3 axes for X,Y, and Z polarization choices for dipole vector component

import sys
from DipPolAnalyzer import DipolePolarizabilityFile as DPF
from PlotPowerSpectra import *
import Smoothing
import operator

#import matplotlib.pyplot as plt
# the extents of the x-range
xmin = 1000.0
xmax = 5000.0
c = 29979245800.0		# speed of light (cm/s)
dt = 0.75e-15	# length of time between each simulation data point
correlation_tau = 7000	# length of the correlation function

files = glob.glob('*/sfg.cold.dat')
dpfs = [DPF(f) for f in files]

# load the data file
#dpf = DPF(sys.argv[1])
alphas_xx = [dpf.Alpha(0,0) for dpf in dpfs]
alphas_yy = [dpf.Alpha(1,1) for dpf in dpfs]
mus = [dpf.Mu(2) for dpf in dpfs]

# perform the cross correlation of the polarizability with the dipole in the SSP regime
#ccf = numpy.array([ManualCorrelate(operator.mul, tau, alpha, mu) for tau in range(correlation_tau)])
ccfs_xx = [numpy.array(NewCorr(alpha,mu)[:correlation_tau]) for alpha,mu in zip(alphas_xx,mus)]	 # using the new routine
ccfs_yy = [numpy.array(NewCorr(alpha,mu)[:correlation_tau]) for alpha,mu in zip(alphas_yy,mus)]	 # using the new routine
ccfs = ccfs_xx + ccfs_yy
avg_ccf = numpy.array(reduce(operator.add,ccfs))/2.0/len(ccfs)

# set up the time axis and plot the correlation function
axs = TCFAxis()
axs.plot(range(len(avg_ccf)), avg_ccf, linewidth=2.5, color='k')

# apply a smoothing window to the ccf
window = numpy.hanning(len(avg_ccf))
avg_ccf = window * avg_ccf

# fourier transform the smoothed/periodic correlation function
fft = numpy.array(numpy.fft.fft(avg_ccf))	 # this is a complex-valued function

# define the frequency axis
freqs = numpy.array(numpy.fft.fftfreq(n=len(avg_ccf), d=dt))/c

# apply a prefactor
fft = fft * freqs

# now take the mag squared of the function to get the SFG lineshape
chi_2 = abs(fft) * abs(fft)

# smooth out the chi_2
smooth_chi_2 = Smoothing.window_smooth(chi_2, window_len=10)

# set up the frequency axis/figure
axs = PowerSpectrumAxis()
axs.plot (freqs, smooth_chi_2, linewidth=2.5, color='k')

plt.xlim(0,5000)
plt.show()



class MoritaSFG2002:

	def __init__(self, file, dt=0.75-15, temp=300.0):
		# first open up the bond-length trajectory file
		self.datafile = ColumnDataFile(file)
		self.dipoles = apply (zip, [self.datafile[i] for i in range(3)])
		self.polarizabilities = apply(zip, [self.datafile[i] for i in range(3,12)])

		self.dpa = DipPolAnalyzer(self.dipoles,self.polarizabilities)
		self.dt = dt
		self.temp = temp

	def PrintData(self):

		# plotting out IR spectra
		ir = TCFSFGAnalyzer (self.dpa.IR_TCF())
		ir.CalcSFGChi(self.dt,self.temp)
		ir_x = ir.Freq()
		ir_y = ir.FFTNormSquared()

		sfg_data = []

		pol = [0,1,2]	# polarization combos
		for row in range(3):

			# get both of the symmetric SSP data sets together
			# i.e. s1,s1,p & s2,s2,p
			sfg_tcf = self.dpa.SFG_TCF(pol[0],pol[1],pol[2])
			sfg = TCFSFGAnalyzer (sfg_tcf)
			sfg.CalcSFGChi(self.dt,self.temp)

			sfg_data.append(sfg.SFG())

			pol = pol[1:]+pol[:1]	# rotate to the next P polarization

		# limit the output to be within the frequency extents listed above (xmin, xmax)
		freq_min_index = ir_x.index([f for f in ir_x if f > xmin][0])
		freq_max_index = ir_x.index([f for f in ir_x if f < xmax][-1])
		# now get the data all set up
		data = zip(ir_x, ir_y, sfg_data[0], sfg_data[1], sfg_data[2])

		for d in data[freq_min_index:freq_max_index]:
			for i in d:
				print "%12.6e " % (i),
			print



'''
	def PlotData(self):

		# Set up the plot parameters (labels, size, limits, etc)
		fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
		fig2 = plt.figure(num=2, facecolor='w', edgecolor='w', frameon=True)

		# the extents of the x-range
		xmin = 0.0
		xmax = 10000.0

		# plotting out IR spectra
		ir = TCFSFGAnalyzer (self.dpa.IR_TCF())
		ir.CalcSFGChi(self.dt,self.temp)
		axs = fig2.add_subplot(1,1,1)
		axs.plot(ir.Freq(), ir.ChiSquared(), label='xyz IR')
		axs.set_xlim(1000,4000)
		ShowLegend(axs)

		pol = [0,1,2]	# polarization combos
		for row in range(3):

			axs = fig.add_subplot(3,1,row+1)

			# get both of the symmetric SSP data sets together
			# i.e. s1,s1,p & s2,s2,p
			sfg_tcf = self.dpa.SFG_TCF(pol[0],pol[1],pol[2])
			sfg = TCFSFGAnalyzer (sfg_tcf)
			sfg.CalcSFGChi(self.dt,self.temp)

			# now average them

			axs.plot(sfg.Freq(), sfg.SFG(), label="P = "+str(pol[2]))

			axs.set_xlim(1000,4000)
			ShowLegend(axs)

			pol = pol[1:]+pol[:1]	# rotate to the next P polarization
			sys.stdout.flush()
'''


#sfg = MoritaSFG2002(sys.argv[1], 1.0e-15, 300.0)

#sfg.PrintData()
#sfg.PlotData()
#plt.show()
