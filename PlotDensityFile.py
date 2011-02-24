import sys
from glob import glob
from ColumnDataFile import ColumnDataFile as CDF
from StatisticMachine import StatisticMachine as SM
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
#from pylab import *
from numpy import array, arange
#from scipy import interpolate
#from operator import itemgetter
import PlotUtility

sm = SM()

# load all of the files for the analysis and parse out the columns of data
file = glob(sys.argv[1])
cdf = CDF(file[0],header=True)

fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)
axs = fig.add_subplot(1,1,1)

dr = 0.25					# thickness of each slice of the system
const = 30.0*30.0*dr*6.02e23*1.0e-24 				# conversion factor to g/mL
mass = {'S':64.0644, 'O':18.002}
axs.plot(cdf['position'], array(cdf['O'])*mass['O']/const, color='k', linestyle=':', label=r'H$_2$O', linewidth=2.0)
S_y = array(cdf['S'])*mass['S']/const 
axs.plot(cdf['position'], S_y, color='r', label=r'SO$_2$', linewidth=3.5)


p0 = array ([1.0, -40.0, 3.0, 3.0, 7.0, 7.0])
residual = sm.FittingFunction(sm.residuals,function=sm.double_gaussian_fit)
plsq = leastsq(residual, p0, args=(S_y,cdf['position']), maxfev=10000)
print plsq[0]
gauss_y = sm.double_gaussian_fit(array(cdf['position']),plsq[0])
axs.plot(cdf['position'], gauss_y, linewidth=1.5, color='k', label='Gaussian Fit')

plt.xlim(-50.0,20.0)
plt.xticks(fontsize=44)
plt.xlabel (r'Distance to Top Surface / $\AA$', fontsize=58)

plt.ylim(0.0,1.2)
plt.yticks(fontsize=44)
plt.ylabel (r'Density / $\frac{g}{mL}$', fontsize=58)

PlotUtility.ShowLegend(axs)
plt.show()
