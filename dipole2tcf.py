import numpy
from ColumnDataFile import ColumnDataFile as CDF

# returns an autocorrelation of the given sample
def autocorr(x):
	result = numpy.correlate(x, x, mode='full')
	return result[result.size/2:]

# performs an ensemble average to account for the number of data points used in the autocorrelation function
def ensemble_avg (x):
	return map(lambda x,tau: float(x)/float(len(x)-tau), x, range(len(x)))



data = CDF('tcf.dat')[0]
acf = autocorr(data)
acf_avg = ensemble_avg(acf)

