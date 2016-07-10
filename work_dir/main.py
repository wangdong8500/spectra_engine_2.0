"""
main program for generating transit or secondary transit spectra
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('/Users/wangdong/Documents/research/spectra_project/noise_engine_2.0')
from spectra import TransitSpectra,SecondarySpectra

# input parameters
# distance in pc
d = 5
# total integration time
t = 11000 # about 3 hrs
# binned spectra resolution
R = 100
# ccd transmission efficiency
tau = 0.5
# radius of Jupiter
R_p = 7e+7
# telescope aperture 
A_tel = 25
# number of integrations
n_ints = 20
# number of instances used in simulating noise
n_instance = 100
# noise floor
noise_floor = [20.e-6,30.e-6,50.e-6]
# range of wavelength 
binned_wavelength_min = 2.
binned_wavelength_max = 11.
# output
save_fig = True
save_data = True

"""
Transit calculations
"""

####################################################################################################################
# input file 1
transit_file_directory_1 = '/Users/wangdong/Documents/research/spectra_project/spectra_data/CNOSP_500K_G_star_1e9/include_all/res_100/'
transit_file_name_1 = 'spectra-CNOSP500.TRANSIT-IR-Res100' 
transit_file_path_1 = transit_file_directory_1 + transit_file_name_1

# input file 2
transit_file_directory_2 = '/Users/wangdong/Documents/research/spectra_project/spectra_data/CNOSP_500K_G_star_1e9/no_PH3/res_100/'
transit_file_name_2 = 'spectra-CNOSP500.TRANSIT-IR-noPH3.Res100'  
transit_file_path_2 = transit_file_directory_2 + transit_file_name_2

# output file name
transit_output_file_name = 'spectra_500K_transit_PH3'

#####################################################################################################################
# compute the transit spectra 

transit_spectra_1 = TransitSpectra()
transit_spectra_1.Compute(transit_file_path_1,d,\
	A_tel,tau,R,noise_floor,binned_wavelength_min,binned_wavelength_max,\
	n_ints,t,n_instance)

transit_spectra_2 = TransitSpectra()
transit_spectra_2.Compute(transit_file_path_2,d,\
	A_tel,tau,R,noise_floor,binned_wavelength_min,binned_wavelength_max,\
	n_ints,t,n_instance)

# plot the transit spectra
fig = plt.figure()
ax = fig.add_subplot(111)
transit_spectra_1.Plot(ax,lw = 2,color = 'b')
transit_spectra_2.Plot(ax,lw = 2,color = 'r')
ax.set_xscale('log')
ax.set_xlim(2,11)
ax.set_xlabel('wavelength($\mu$m)')
ax.set_ylabel('transit depth(%)')
ax.set_xticks([2,3,4,5,6,7,8,9,10,11])
ax.set_xticklabels([2,3,4,5,6,7,8,9,10,11])
if save_fig == True:
	plt.savefig('./spectra/'+transit_output_file_name+'.ps')
	print "transit spectra generated successfully!"
else:
	plt.show()
if save_data == True:
	data = np.column_stack((transit_spectra_1.mean_depth,transit_spectra_1.sigma,transit_spectra_2.mean_depth,transit_spectra_2.sigma))
	np.savetxt('./data/'+transit_output_file_name,data)
	print "transit spectra data saved successfully!"

########################################################################################################################

"""
Secondary eclipse calculations
"""

####################################################################################################################
# input file 1
emergent_file_directory_1 = '/Users/wangdong/Documents/research/spectra_project/spectra_data/CNOSP_500K_G_star_1e9/include_all/res_100/'
emergent_file_name_1 = 'spectra-CNOSP500.IR.Res100.dat' 
emergent_file_path_1 = emergent_file_directory_1 + emergent_file_name_1

# input file 2
emergent_file_directory_2 = '/Users/wangdong/Documents/research/spectra_project/spectra_data/CNOSP_500K_G_star_1e9/no_PH3/res_100/'
emergent_file_name_2 = 'spectra-CNOSP500.IR.noPH3.Res100.dat' 
emergent_file_path_2 = emergent_file_directory_2 + emergent_file_name_2

# output file
secondary_output_file_name = 'spectra_500K_secondary_PH3'

#####################################################################################################################
# compute the secondary eclipse spectra 

secondary_spectra_1 = SecondarySpectra()
secondary_spectra_1.Compute(emergent_file_path_1,d,R_p,\
	A_tel,tau,R,noise_floor,binned_wavelength_min,binned_wavelength_max,\
	n_ints,t,n_instance)

secondary_spectra_2 = SecondarySpectra()
secondary_spectra_2.Compute(emergent_file_path_2,d,R_p,\
	A_tel,tau,R,noise_floor,binned_wavelength_min,binned_wavelength_max,\
	n_ints,t,n_instance)

# plot the transit spectra
fig = plt.figure()
ax = fig.add_subplot(111)
secondary_spectra_1.Plot(ax,lw = 2,color = 'b')
secondary_spectra_2.Plot(ax,lw = 2,color = 'r')
ax.set_xscale('log')
ax.set_xlim(2,11)
ax.set_xlabel('wavelength($\mu$m)')
ax.set_ylabel('f_p/f_s')
ax.set_xticks([2,3,4,5,6,7,8,9,10,11])
ax.set_xticklabels([2,3,4,5,6,7,8,9,10,11])
if save_fig == True:
	plt.savefig('./spectra/'+secondary_output_file_name+'.ps')
	print "secondary eclipse spectra generated successfully!"
else:
	plt.show()
if save_data == True:
	data = np.column_stack((secondary_spectra_1.mean_ratio,secondary_spectra_1.sigma,secondary_spectra_2.mean_ratio,secondary_spectra_2.sigma))
	np.savetxt('./data/'+secondary_output_file_name,data)
	print "secondary eclipse spectra data saved successfully!"

