"""
This module computes binned flux from source, background and instruments noise
"""
import numpy as np 
from scipy import interpolate

h = 6.626e-34
c = 3.0e8

class Binning:
	"""
	binning the spectra according to the spectra resolution
	"""
	def __init__(self,wavelength,binned_wavelength_min,binned_wavelength_max,R):
		self.wavelength = wavelength
		N = np.int(np.log(binned_wavelength_max/binned_wavelength_min)*R)
		self.binned_wavelength = np.logspace(np.log10(binned_wavelength_min),np.log10(binned_wavelength_max),N)

	def binning(self,flux):
		f = interpolate.interp1d(self.wavelength[::-1],flux[::-1])
		return f(self.binned_wavelength)

class BinnedTransitPhotonNumber:
	"""
	compute the binned photon number for transit observations
	"""
	def __init__(self,binned_wavelength_min,binned_wavelength_max,transit_flux,R,tau,A_tel,t):
		# convert to meters
		binning_obj = Binning(transit_flux.wavelength,binned_wavelength_min,binned_wavelength_max,R)
		self.binned_wavelength = binning_obj.binned_wavelength
		f_star_bin = binning_obj.binning(transit_flux.f_star)
		f_in_transit_bin = binning_obj.binning(transit_flux.f_in_transit) 
		binned_wavelength = self.binned_wavelength*1e-6
		self.n_star_bin = self.source_photon_number_eqn(binned_wavelength,f_star_bin,R,tau,A_tel)*t
		self.n_in_transit_bin = self.source_photon_number_eqn(binned_wavelength,f_in_transit_bin,R,tau,A_tel)*t

	def source_photon_number_eqn(self,wavelength, F, R, tau, A_tel):
		# Source returns the number of photons in each wavelength bin 
		# F: flux 
		# t: integration time
		# wavelength: wavelength bin grid, should be an array
		# R: spectra resolution
		# tau: total system transmission
		# A_tel: telescope collecting area
		return F*A_tel*wavelength**2/(h*c*R)*tau


class BinnedSecondaryPhotonNumber:
	"""
	compute the binned photon number for secondary eclipse observations
	"""
	def __init__(self,binned_wavelength_min,binned_wavelength_max,secondary_flux,R,tau,A_tel,t):
		# convert to meters
		binning_obj = Binning(secondary_flux.wavelength,binned_wavelength_min,binned_wavelength_max,R)
		self.binned_wavelength = binning_obj.binned_wavelength
		f_star_bin = binning_obj.binning(secondary_flux.f_star)
		f_out_transit_bin = binning_obj.binning(secondary_flux.f_out_transit) 
		binned_wavelength = self.binned_wavelength*1e-6
		self.n_star_bin = self.source_photon_number_eqn(binned_wavelength,f_star_bin,R,tau,A_tel)*t
		self.n_out_transit_bin = self.source_photon_number_eqn(binned_wavelength,f_out_transit_bin,R,tau,A_tel)*t

	def source_photon_number_eqn(self,wavelength, F, R, tau, A_tel):
		# Source returns the number of photons in each wavelength bin 
		# F: flux 
		# t: integration time
		# wavelength: wavelength bin grid, should be an array
		# R: spectra resolution
		# tau: total system transmission
		# A_tel: telescope collecting area
		return F*A_tel*wavelength**2/(h*c*R)*tau


class BinnedJWSTBackgroundPhotonEnergy:
	"""

	"""
	def __init__(self,wavelength,R,t):
		self.wavelength = wavelength
		self.n_background = self.JWST_Bkg_flux(self.wavelength,R)*t 

	# utility functions
	def JWST_Bkg_flux(self,wavelength,R):
		"""
		for a given wavelength, return a background flux
		four modes specifically corresponds to JWST transit observations 

		NIRISS 
		1 - 2.5 micron
		instrument: GR700XD 
		F115W (1-1.3 micron): B = 197.30 e/s/arcsec^2
		F150W (1.3 - 1.7 micron): B = 237.24 e/s/arcsec^2 
		F200W (1.7 - 2.5 micron): B = 193.07 e/s/arcsec^2  
		sum of the values to represent bkg signal for GR700XD

		NIRCam
		2.5 - 3.9 micron
		instrument: F322W2
		F277W (2.5 - 3.1 micron): B = 107.57 e/s/arcsec^2 
		F356W (3.1 - 3.9 micron): B = 96.05 e/s/arcsec^2 
		sum of the values to present bkg for F322W2

		NIRCam
		3.9 - 5.0 micron
		F444W (3.9 - 5.0 micron): B = 308.74 e/s/arcsec^2 

		MIRI 
		5.0 - 11 micron
		instrument: LRS prism
		F560W (5.0 - 6.2 micron): B = 213.69 e/s/arcsec^2 
		F770W (6.2 - 8.8 micron): B = 1970.02 e/s/arcsec^2 
		F1000W (8.8 - 11 micron): B = 2972.29 e/s/arcsec^2  
		sum of these values to represent bkg for LRS prism
		
		wavlength must be a numpy array object
		R is the spectra resolution
		return the background photon flux in each spectra bin
		"""
		n = len(wavelength)
		bkg = np.zeros(n)

		# For NIRISS, bright SOSS mode, GR700XD optics
		index_1 = np.nonzero(np.logical_and(wavelength>=1, wavelength < 2.5))
		B1 = 197.30
		B2 = 237.24
		B3 = 193.07
		B_GR700XD = B1 + B2 + B3
		A_pix_NIRISS = 0.065*0.065
		n_pix_NIRISS = 25*2
		R_native_NIRISS = 700
		bkg[index_1] = self.bkg_flux_eqn(B_GR700XD,A_pix_NIRISS,n_pix_NIRISS,R_native_NIRISS,R)

		# for NIRCam, LW grism mode, F322W2 optics and LW grism, F444W
		index_2 = np.nonzero(np.logical_and(wavelength>=2.5, wavelength < 3.9))
		index_3 = np.nonzero(np.logical_and(wavelength>=3.9, wavelength < 5.0))
		B4 = 107.57
		B5 = 96.05
		B_F322W2 = B4 + B5
		B6 = 308.74
		B_F444W = B6	
		A_pix_NIRCam = 0.064*0.064
		n_pix_NIRCam = 2*2
		R_native_NIRCam = 1700
		bkg[index_2] = self.bkg_flux_eqn(B_F322W2,A_pix_NIRCam,n_pix_NIRCam,R_native_NIRCam,R)
		bkg[index_3] = self.bkg_flux_eqn(B_F444W,A_pix_NIRCam,n_pix_NIRCam,R_native_NIRCam,R)

		# for MIRI, silitless mode, LRS prism optics
		index_4 = np.nonzero(np.logical_and(wavelength>=5.0, wavelength < 11))
		B7 = 213.69
		B8 = 1970.02
		B9 = 2972.29
		B_LRS_prism = B7 + B8 + B9
		A_pix_MIRI = 0.110*0.110
		n_pix_MIRI = 2*2
		R_native_MIRI = 100
		bkg[index_4] = self.bkg_flux_eqn(B_LRS_prism,A_pix_MIRI,n_pix_MIRI,R_native_MIRI,R)

		return bkg

	def bkg_flux_eqn(self,B,A_pix,n_pix,R_native,R):
		# Bkg returns the background signal in each spectra bin
		# B: the backgroud of the intrument mode in electrons per arcsec^2 per second
		# t: total exposure time
		# A_pix: the area subtended by each pixel in arcsec^2
		# n_pix: the numeber of spatial x 2 spectra pixels smmed in each R_native 
		#        native 2-pixel resolution element of the selected observing mode
		# R: final binned spectra resolution
		return B*A_pix*n_pix*R_native/R


