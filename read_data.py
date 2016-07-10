"""
This module reads the spectra data into the memory
perform preliminary processing to get the source flux
received by the telescope
three components:
star only
star + planet
star obscured by planet
"""
import numpy as np 
import matplotlib.pyplot as plt 

class TransitData:
	def __init__(self,wavelength=None,depth=None):
		self.wavelength = wavelength
		self.depth = depth

	def Read(self,transit_file_path):
		[self.wavelength, self.depth] = np.loadtxt(transit_file_path, \
		comments='#', delimiter=None, \
		skiprows=0, usecols=(0,1), unpack=True)


	def Plot(self,**kwargs):
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.plot(self.wavelength,self.depth,**kwargs)
		return ax

class EmergentData:
	def __init__(self,wavelength=None,flux_per_wavelength=None):
		self.wavelength = wavelength
		self.flux_per_wavelength = flux_per_wavelength

	def Read(self,emergent_file_path):
		[self.wavelength, self.flux_per_wavelength] = np.loadtxt(emergent_file_path, \
		comments='#', delimiter=None, \
		skiprows=0, usecols=(0,1), unpack=True)


	def Plot(self,**kwargs):
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.plot(self.wavelength,self.flux_per_wavelength,**kwargs)
		return ax

if __name__ == '__main__':
	"""
	test the module
	"""
	path = '/Users/wangdong/Documents/research/spectra_project/spectra_data/'+\
	'CNOSP_500K_G_star_1e9/include_all/res_100/spectra-CNOSP500.TRANSIT-IR-Res100'
	transit_data = TransitData()
	transit_data.Read(path)
	ax = transit_data.Plot(lw = 2)
	ax.set_xlabel('wavelength($\mu$m)')
	ax.set_ylabel('transit depth')
	ax.set_xlim(2,20)
	ax.set_xscale('log')
	plt.show()

	path = '/Users/wangdong/Documents/research/spectra_project/spectra_data/'+\
	'CNOSP_500K_G_star_1e9/include_all/res_100/spectra-CNOSP500.IR.Res100.dat'
	emergent_data = EmergentData()
	emergent_data.Read(path)
	ax = emergent_data.Plot(lw = 2)
	ax.set_xlabel('wavelength($\mu$m)')
	ax.set_ylabel('flux per wavelength(W/m^2/m)')
	ax.set_xscale('log')
	ax.set_xlim(2,20)
	plt.show()




