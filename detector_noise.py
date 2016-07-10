import numpy as np
class DetectorNoise:
	"""
	compute the total detector noise as a function of wavelength
	"""
	def __init__(self,wavelength,R,n_ints):
		self.wavelength = wavelength
		self.sigma = self.JWST_total_detector_noise(self.wavelength,R,n_ints)

	# utility functions
	def JWST_total_detector_noise(self,wavelength, R, n_ints):
		n = len(wavelength)
		N_d_total = np.zeros(n)
		# For NIRISS
		index_1 = np.nonzero(np.logical_and(wavelength>=1, wavelength < 2.5))
		N_d_NIRISS = 18
		n_pix_NIRISS = 25*2
		R_native_NIRISS = 700
		N_d_total[index_1] = self.detector_noise_eqn(N_d_NIRISS,n_pix_NIRISS,n_ints,R_native_NIRISS,R)

		# for NIRCam
		index_2 = np.nonzero(np.logical_and(wavelength>=2.5, wavelength < 5))
		N_d_NIRCam = 18
		n_pix_NIRCam = 2*2
		R_native_NIRCam = 1700
		N_d_total[index_2] = self.detector_noise_eqn(N_d_NIRCam,n_pix_NIRCam,n_ints,R_native_NIRCam,R)

		# for MIRI
		index_3 = np.nonzero(np.logical_and(wavelength>=5., wavelength < 11.))
		N_d_NIRCam = 28
		n_pix_NIRCam = 2*2
		R_native_NIRCam = 100
		N_d_total[index_3] = self.detector_noise_eqn(N_d_NIRCam,n_pix_NIRCam,n_ints,R_native_NIRCam,R)

		return N_d_total


	def detector_noise_eqn(self,N_d,n_pix,n_ints,R_native,R):
		# N_d: total detector noise of a single integration
		# n_pix: the numeber of spatial x 2 spectra pixels smmed in each R_native 
		#        native 2-pixel resolution element of the selected observing mode
		# n_ints: number of integrations during exposure time t
		# R: final binned spectra resolution
		return N_d*np.sqrt(n_pix*n_ints*R_native/R)