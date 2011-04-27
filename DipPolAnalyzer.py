from ColumnDataFile import ColumnDataFile as CDF
from numpy import mat, array
from Utility import group
import sys

class DipolePolarizabilityFile:
	def __init__(self, filename):
		cdf = CDF(filename)
		self.rho = zip(cdf[0], cdf[1], cdf[2])
		self.rho = [array(r) for r in self.rho]

		self.alpha = zip(cdf[3], cdf[4], cdf[5], cdf[6], cdf[7], cdf[8], cdf[9], cdf[10], cdf[11])
		self.alpha = [apply(mat,[group(a,3)]) for a in self.alpha]
		
	''' The dipoles parsed from the file '''
	def Dipoles(self):
		return self.rho
	''' The polarizability matrices parsed from the file '''
	def Polarizabilities(self):
		return self.alpha

	''' Returns a single element of all the dipoles '''
	def Mu(self,i):
		mu = [r[i] for r in self.rho]
		return mu
	''' returns a particular polarization of the polarizability '''
	def Alpha(self,i,j):
		al = [a[i,j] for a in self.alpha]
		return al

dpf = DipolePolarizabilityFile(sys.argv[1])
		
	

class DipPolAnalyzer:
	def __init__(self,rho,alpha=None):
	  	#print rho[0][0]
		self.rho = [array(r) for r in rho]
		#print self.rho[:10]
		if alpha != None:
			self.alpha = [apply(mat,[group(a,3)]) for a in alpha]

	def Rho(self):
		return self.rho
	def Alpha(self):
		return self.alpha

	def IR_TCF (self):
		ret = []
		try: 
			ret = [dot(i,self.rho[0]) for i in self.rho]
		except ValueError:
			print self.rho[0]

		return ret

	def SFG_TCF(self,s1,s2,p):
		alpha = [(a[s1,s1] + a[s2,s2])/2.0 for a in self.alpha]
		#alpha = [a[s1,s1] for a in self.alpha]
		rho = self.rho[0][p]
		return [a*rho for a in alpha]
		#return [dot(a,self.rho[0][r]) for a in self.alpha]

