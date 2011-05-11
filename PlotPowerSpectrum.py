#!/usr/bin/python

# This SFG calculator creates both IR and SFG spectra from dipole/polarizability files
# The SFG spectra are all SSP polarized, but that can be changed in the DipPol analyzer

# PrintData prints out 5 columns for frequency and spectral intensity:
#		frequency, IR, SFG_X, SFG_Y, SFG_Z

# PlotData creates 2 figures, one for IR and one for SFG. The SFG figure has 3 axes for X,Y, and Z polarization choices for dipole vector component

import sys
from ColumnDataFile import ColumnDataFile as CDF
from PlotPowerSpectra import *
import Smoothing
import operator
import numpy

#import matplotlib.pyplot as plt
# the extents of the x-range
xmin = 1000.0
xmax = 5000.0
c = 29979245800.0		# speed of light (cm/s)
dt = 0.75e-15	# length of time between each simulation data point
correlation_tau = 3000

files = sys.argv[1:]
cdfs = [CDF(f) for f in files]

tcfs_x = [NewCorr(cdf[0])[:correlation_tau] for cdf in cdfs]
tcfs_y = [NewCorr(cdf[1])[:correlation_tau] for cdf in cdfs]
tcfs_z = [NewCorr(cdf[2])[:correlation_tau] for cdf in cdfs]
tcfs = tcfs_x + tcfs_y + tcfs_z
print len(tcfs)

# set up the time axis and plot the correlation function
axs = TCFAxis()
for t in tcfs:
	axs.plot(range(len(t)), t, linewidth=2.5)

# apply a smoothing window to the tcf
smoothed_tcfs = [numpy.hanning(len(t)) * t for t in tcfs]

# fourier transform the smoothed/periodic correlation function
ffts = [numpy.array(numpy.fft.fft(t)) for t in smoothed_tcfs]	 # this is a complex-valued function

# define the frequency axis
freqs = [numpy.array(numpy.fft.fftfreq(n=len(t), d=dt))/c for t in smoothed_tcfs]

# apply a prefactor (multiply by the frequency)
ffts = [f*w for f,w in zip(ffts,freqs)]

# now take the mag squared of the function to get the SFG lineshape
lineshapes = [abs(f)*abs(f) for f in ffts]

# set up the frequency axis/figure
axs = PowerSpectrumAxis()
for f,i in zip(freqs,lineshapes):
	axs.plot (f, i, linewidth=2.5)

plt.xlim(0,5000)
plt.show()
