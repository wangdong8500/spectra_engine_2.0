"""
This module computes the flux per wavelength received by the telescope 
"""
import numpy as np 
from star import Star

# 1 parsec to meters
pc = 3.086e+16

class TransitFlux:
	def __init__(self,transit_data,d,L_star):
		"""
		d is the distance from planeary system to the telescope in pc
		L_star: W/m

		return
		wavelength in micron 
		flux in W/m^2/m
		"""
		self.wavelength = transit_data.wavelength
		depth = transit_data.depth
		self.f_star = L_star/(4*np.pi*(d*pc)**2)
		self.f_in_transit = self.f_star*(1-depth)

class SecondaryFlux:
	def __init__(self,emergent_data,d,R_p,L_star):
		"""
		R_p is the planetary radius
		
		return
		wavelength in micron 
		flux in W/m^2/m		"""
		self.wavelength = emergent_data.wavelength
		self.f_star = L_star/(4*np.pi*(d*pc)**2)
		planetary_flux = emergent_data.flux_per_wavelength
		L_p = planetary_flux*(4*np.pi*R_p*R_p)
		self.f_out_transit = (L_p + L_star)/(4*np.pi*(d*pc)**2)

