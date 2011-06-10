#!/usr/bin/python

# This SFG calculator creates both IR and SFG spectra from dipole/polarizability files
# The SFG spectra are all SSP polarized, but that can be changed in the DipPol analyzer

# PrintData prints out 5 columns for frequency and spectral intensity:
#		frequency, IR, SFG_X, SFG_Y, SFG_Z

# PlotData creates 2 figures, one for IR and one for SFG. The SFG figure has 3 axes for X,Y, and Z polarization choices for dipole vector component

import glob
from DipPolAnalyzer import DipolePolarizabilityFile as DPF
from PlotPowerSpectra import *
import PlotUtility
import Smoothing
import operator
from pylab import *

#import matplotlib.pyplot as plt
# the extents of the x-range
tau = 20000	# length of the correlation function


def PlotMorita (files,axs):
	dpfs = [DPF(f) for f in files]
	
	# load the data file
	#dpf = DPF(sys.argv[1])
	alphas_xx = [dpf.Alpha(0,0) for dpf in dpfs]
	alphas_yy = [dpf.Alpha(1,1) for dpf in dpfs]
	alphas_xy = [dpf.Alpha(0,1) for dpf in dpfs]
	alphas_yx = [dpf.Alpha(1,0) for dpf in dpfs]
	mus = [dpf.Mu(2) for dpf in dpfs]
	
	# perform the cross correlation of the polarizability with the dipole in the SSP regime
	ccfs_xx = [numpy.array(Correlate(alpha,mu)) for alpha,mu in zip(alphas_xx,mus)]	 # using the new routine
	ccfs_yy = [numpy.array(Correlate(alpha,mu)) for alpha,mu in zip(alphas_yy,mus)]	 # using the new routine
	ccfs_xy = [numpy.array(Correlate(alpha,mu)) for alpha,mu in zip(alphas_xy,mus)]	 # using the new routine
	ccfs_yx = [numpy.array(Correlate(alpha,mu)) for alpha,mu in zip(alphas_yx,mus)]	 # using the new routine
	#ccfs = ccfs_xx + ccfs_yy + ccfs_xy + ccfs_yx
	ccfs = ccfs_xx + ccfs_yy

	avg_ccf = numpy.array(reduce(operator.add,ccfs))/len(ccfs)
	
	# set up the time axis and plot the correlation function
	#axs = TCFAxis(1)
	#axs.plot(range(len(avg_ccf)), avg_ccf, color='k')

	freqs,spectrum,smooth_spectrum = PowerSpectrum(avg_ccf)
  	smooth_spectrum = smooth_spectrum/smooth_spectrum.max()

	axs.plot (freqs[300:-300], smooth_spectrum[300:-300], linewidth=2.5)
	plt.xlim (500,4000)
	plt.ylim (-0.05,1.2)

	
files_cold = glob.glob('[1-5]/sfg.top10.dat')
files_hot = glob.glob('[6-9]/sfg.top10.dat')

axs = PowerSpectrumAxis(2)
PlotMorita(files_cold,axs)
PlotMorita(files_hot,axs)

xticks(fontsize=48)
yticks([])

plt.show()
