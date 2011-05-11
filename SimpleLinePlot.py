import sys
import matplotlib.pyplot as plt
from pylab import *
import PlotUtility
from ColumnDataFile import ColumnDataFile as CDF
import numpy

cdf = CDF(sys.argv[1])

fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
axs = fig.add_subplot(1,1,1)

axs.set_xlabel(r'Time / ns', fontsize='64')
axs.set_ylabel('S-O Bondlength / $\AA$', fontsize='64')
xticks(fontsize=40)
yticks(fontsize=40)


time = numpy.array(range(len(cdf[0])))/0.75/1000.0

axs.plot(time, cdf[0]);

plt.show()
