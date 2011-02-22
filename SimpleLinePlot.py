import sys
import matplotlib.pyplot as plt
import PlotUtility
from ColumnDataFile import ColumnDataFile as CDF
from pylab import *

cdf = CDF(sys.argv[1])

fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
axs = fig.add_subplot(1,1,1)

axs.set_ylabel(r'Binding Water Distance Above Surface / $\AA$', fontsize='100')
axs.set_xlabel('Sulfur dioxide - Water intermolecular distance / $\AA$', fontsize='100')

xticks(fontsize=64)
yticks(fontsize=64)

axs.plot(cdf[0],cdf[1])

plt.show()
