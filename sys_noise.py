import numpy as np 

class SysNoise:
	def __init__(self,wavelength,noise_floor):
		self.wavelength = wavelength
		self.ratio = self.SystematicNoise(self.wavelength,noise_floor)

	def SystematicNoise(self,wavelength,noise_floor):
		# NIRISS: 20 ppm
		# NIRCam: 30 ppm
		# MIRI: 50 ppm
		n = len(wavelength)
		N_S_sys = np.zeros(n)
		
		index_1 = np.nonzero(np.logical_and(wavelength>=1, wavelength < 2.5))
		N_S_sys[index_1] = noise_floor[0]

		index_2 = np.nonzero(np.logical_and(wavelength>=2.5, wavelength < 5.0))
		N_S_sys[index_2] = noise_floor[1]

		index_3 = np.nonzero(np.logical_and(wavelength>=5.0, wavelength < 11.0))
		N_S_sys[index_3] = noise_floor[2]

		return N_S_sys

