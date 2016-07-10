"""
This module computes the property of the star alone
"""
import numpy as np 

class Star:
	def __init__(self,T=None,R=None):
		# define the planet by temperature and radius
		# default to be solar-like star
		if T == None:
			self.T = 5750
		else:
			self.T = T 
		if R == None:
			self.R = 6.955e8
		else:
			self.R = R

	def StarSpectra(self,wavelength,file_path=None):
		"""
		input: wavelength grid in micron
		input: path to the file containing the spectra
		return the luminosity per wavelength
		"""
		pass

	def BlackBodySpectra(self,wavelength):
		"""
		input: wavelength is a vector in micron
		compute blackbody radiation and
		return the luminosity per wavelength (W/m)
		"""
		h = 6.626e-34
		c = 3e+8
		k = 1.38e-23
		# convert from micron to meters
		wavelength = wavelength*1.e-6
		flux_bd = 2*np.pi*h*c*c/(wavelength**5)/(np.exp(h*c/(wavelength*k*self.T))-1.)
		L_bd = flux_bd*4*np.pi*self.R*self.R
		return L_bd

if __name__ == "__main__":
	# wavelegnth in microns
	wavelength = np.logspace(start=-0.5, stop=2, num=500, endpoint=True, base=10.0)
	star = Star()
	L_bd = star.BlackBodySpectra(wavelength)
	import matplotlib.pyplot as plt
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot(wavelength,L_bd)
	ax.set_xlabel('wavelength($\mu$m)')
	ax.set_ylabel('luminosity (W/m)')
	ax.set_yscale('log')
	ax.set_xlim(0.5,20)
	plt.show()
