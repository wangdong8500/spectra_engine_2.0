import matplotlib.pyplot as plt 
import numpy as np 
from read_data import TransitData, EmergentData
from star import Star
from source_flux import TransitFlux, SecondaryFlux
from binned_photon_energy import BinnedTransitPhotonNumber, BinnedSecondaryPhotonNumber, BinnedJWSTBackgroundPhotonEnergy
from detector_noise import DetectorNoise
from sys_noise import SysNoise
from noise_engine import NoiseEngine

class TransitSpectra:
	"""
	compute the transit spectra and store the spectra
	"""
	def __init__(self,binned_wavelength = None, depth = None):
		self.binned_wavelength=binned_wavelength
		self.mean_depth = depth

	def Compute(self,transit_file_path,d,\
		A_tel,tau,R,noise_floor,binned_wavelength_min,binned_wavelength_max,\
		n_ints,t,n_instance):
		"""
		d is the distance in parsec from the planetary system to telescope
		R_p: planetary radius
		tau: transmission efficiency
		R is the spectra resolution
		n_ints: number of integrations
		t is the total integration time
		number of instances for simulating noise
		return the spectra
		"""
		# read transit data
		transit_data = TransitData()
		transit_data.Read(transit_file_path)
		# star luminosity per wavelength
		star = Star()
		L_star = star.BlackBodySpectra(transit_data.wavelength)
		# construct transit flux object
		transit_flux = TransitFlux(transit_data,d,L_star)
		# construct spectra object
		source_obj = BinnedTransitPhotonNumber(binned_wavelength_min,binned_wavelength_max,transit_flux,R,tau,A_tel,t)
		bkg_obj = BinnedJWSTBackgroundPhotonEnergy(source_obj.binned_wavelength,R,t)
		detector_noise = DetectorNoise(source_obj.binned_wavelength,R,n_ints)
		sys_noise = SysNoise(source_obj.binned_wavelength,noise_floor)

		# generate n noisy star signal
		star_noise_engine = NoiseEngine(source_obj.n_star_bin,bkg_obj.n_background,detector_noise,sys_noise)
		noisy_star_signals = star_noise_engine.generate(n_instance)

		# generate n noisy in transit signal
		in_transit_noise_engine = NoiseEngine(source_obj.n_in_transit_bin,bkg_obj.n_background,detector_noise,sys_noise)
		noisy_in_transit_signals = in_transit_noise_engine.generate(n_instance)

		# compute the spectra
		depth = []
		depth_squared = []
		for i in range(n_instance):
			depth_instance = (noisy_star_signals[i] - noisy_in_transit_signals[i])/(noisy_star_signals[i]-bkg_obj.n_background)
			depth.append(depth_instance)
			depth_squared.append(depth_instance*depth_instance)
		self.mean_depth = sum(depth)/len(depth)
		self.sigma = np.sqrt(sum(depth_squared)/len(depth) - self.mean_depth*self.mean_depth)
		self.binned_wavelength = source_obj.binned_wavelength

	def Plot(self,ax,**kwargs):
		ax.errorbar(self.binned_wavelength,self.mean_depth,self.sigma,**kwargs)
		return ax

class SecondarySpectra:
	"""
	compute the transit spectra and store the spectra
	"""
	def __init__(self,binned_wavelength = None, ratio = None):
		self.binned_wavelength=binned_wavelength
		self.mean_ratio = ratio

	def Compute(self,emergent_file_path,d,R_p,\
		A_tel,tau,R,noise_floor,binned_wavelength_min,binned_wavelength_max,\
		n_ints,t,n_instance):
		"""
		d is the distance in parsec from the planetary system to telescope
		R_p: planetary radius
		tau: transmission efficiency
		R is the spectra resolution
		n_ints: number of integrations
		t is the total integration time
		number of instances for simulating noise
		return the spectra
		"""
		# read transit data
		emergent_data = EmergentData()
		emergent_data.Read(emergent_file_path)
		# star luminosity per wavelength
		star = Star()
		L_star = star.BlackBodySpectra(emergent_data.wavelength)
		# construct transit flux object
		secondary_flux = SecondaryFlux(emergent_data,d,R_p,L_star)
		# construct spectra object
		source_obj = BinnedSecondaryPhotonNumber(binned_wavelength_min,binned_wavelength_max,secondary_flux,R,tau,A_tel,t)
		bkg_obj = BinnedJWSTBackgroundPhotonEnergy(source_obj.binned_wavelength,R,t)
		detector_noise = DetectorNoise(source_obj.binned_wavelength,R,n_ints)
		sys_noise = SysNoise(source_obj.binned_wavelength,noise_floor)

		# generate n noisy star signal
		star_noise_engine = NoiseEngine(source_obj.n_star_bin,bkg_obj.n_background,detector_noise,sys_noise)
		noisy_star_signals = star_noise_engine.generate(n_instance)

		# generate n noisy in transit signal
		out_transit_noise_engine = NoiseEngine(source_obj.n_out_transit_bin,bkg_obj.n_background,detector_noise,sys_noise)
		noisy_out_transit_signals = out_transit_noise_engine.generate(n_instance)

		# compute the spectra
		ratio = []
		ratio_squared = []
		for i in range(n_instance):
			ratio_instance = (noisy_out_transit_signals[i] - noisy_star_signals[i])/(noisy_star_signals[i]-bkg_obj.n_background)
			ratio.append(ratio_instance)
			ratio_squared.append(ratio_instance*ratio_instance)
		self.mean_ratio = sum(ratio)/len(ratio)
		self.sigma = np.sqrt(sum(ratio_squared)/len(ratio) - self.mean_ratio*self.mean_ratio)
		self.binned_wavelength = source_obj.binned_wavelength

	def Plot(self,ax,**kwargs):
		ax.errorbar(self.binned_wavelength,self.mean_ratio,self.sigma,**kwargs)
		return ax



