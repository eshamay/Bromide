import sys
import glob
from ColumnDataFile import ColumnDataFile as CDF
import matplotlib.pyplot as plt
from pylab import *
import numpy
import operator
import PlotUtility


# translate a numeric value of an SO2 coordination to the equivalent text representation
def CoordinationNumberToName (number):
	name = ''
	for i in range(int(number) / 10):
		name = name + 'S'
	for i in range(int(number) % 10):
		name = name + 'O'

	if number == 0:
		name = 'Unbound'

	return name

def ConvertNumbers(numbers):
	#j = 261
	#print j, 10*(int(j)/100) + (int(j)%10) + (int(j%100)/10)
	return [10*(int(n)/100) + (int(n)%10) + (int(n%100)/10) for n in numbers]

def ExtractCoords(files):
	cdfs = [CDF(f) for f in files]
	data = [c[0] for c in cdfs]

	# convert the abc number to ab
	# a = number of S bonds
	# b,c = numbers of O bonds
	# this converts abc and returns a(b+c)
	data = [ConvertNumbers(d) for d in data]
	data = reduce(operator.add, data)
	# this next line sums the total coordination
	#data = [(int(i) % 10) + (int(i)/10) for i in data] 

	return data
	

def PrintStats (name, data):

	print name
	data = numpy.array(data)

	total = float(data.sum())
	print total
	# Un O OO OOO S SO SOO SOOO SS SSO SSOO SSOOO
	print 'Unbound = ', data[0]/total*100, '%'
	print 'O = ', data[1]/total*100, '%'
	print 'OO = ', data[2]/total*100, '%'
	print 'OOO = ', data[3]/total*100, '%'
	print 'All O = ', (data[1] + data[2] + data[3])/total*100, '%'
	print 'S = ', data[4]/total*100, '%'
	print 'SO = ', data[5]/total*100, '%'
	print 'SOO = ', data[6]/total*100, '%'
	print 'SOOO = ', data[7]/total*100, '%'
	print 'All S = ', (data[4] + data[5] + data[6] + data[7])/total*100, '%'

new_bins = []
num_bins = 4
for i in range(num_bins):
  for j in range(num_bins):
			new_bins.append(i + 10*j)
new_bins.sort()
print new_bins
width = 0.35



#files_cold = glob.glob(sys.argv[1])
files_cold = glob.glob('[1-5]/'+sys.argv[1]+'*')
files_hot = glob.glob('[6-9]/'+sys.argv[1]+'*')
files_hot = files_hot + glob.glob('10/'+sys.argv[1]+'*')

fig = plt.figure(num=1, facecolor='w', edgecolor='w', frameon=True)

data = ExtractCoords(files_cold)

cold_histo, cold_bin_edges = numpy.histogram (data, bins=new_bins)
plt.bar(range(len(new_bins[:-1])), cold_histo, width, color='b', align='center')

PrintStats('Cold', cold_histo)


data = ExtractCoords(files_hot)

hot_histo, hot_bin_edges = numpy.histogram (data, bins=new_bins)
plt.bar([i+width for i in range(len(new_bins[:-1]))], hot_histo, width, color='r', align='center')

PrintStats('Cold', hot_histo)

new_bins = [CoordinationNumberToName(i) for i in new_bins]
xticks(arange(len(new_bins)), new_bins)

plt.show()
