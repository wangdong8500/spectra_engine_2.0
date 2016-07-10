"""
for a given source signal, generate n_instance noisy signals
"""
import numpy as np 

class NoiseEngine:
	def __init__(self,source,bkg,detector_noise,sys_noise):
		"""
		source: a numpy array containing binned number of photons
		bkg: a numpy array containing binned number of photons
		"""
		self.length = source.size
		self.source = source
		self.bkg = bkg
		# total detector noise
		sigma_detector = detector_noise.sigma
		sigma_sys = sys_noise.ratio*source
		self.sigma_total = np.sqrt(source + bkg + sigma_detector**2 + sigma_sys**2)

	def generate(self,n_instance):
		noisy_signal_ensemble = []
		for i in range(n_instance):
			noise_instance = self.sigma_total*np.random.randn(self.length)
			noisy_signal_instance = self.source + self.bkg + noise_instance
			noisy_signal_ensemble.append(noisy_signal_instance)

		return noisy_signal_ensemble


		