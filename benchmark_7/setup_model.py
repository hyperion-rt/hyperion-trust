#################################################
## Setup for TRUST-VII filament benchmark test ##
#################################################

#import sys
import numpy as np
import math
from hyperion.model import Model
from hyperion.dust  import SphericalDust
from hyperion.util.constants import pc, lsun, pi, m_h
from cli import command_line_arguments, handle_command_line_arguments, filename


def setup_model(cli):

	#
	# Hyperion setup:
	#
	model = Model()


	if(cli.mode == "temperature"):
		#
		# Dust properties:
		#
		dust_properties = SphericalDust('dust_integrated_full_scattering.hdf5')


		#
		# Write dust properties:
		#
		dust_properties.write('dust_properties.hdf5')
		dust_properties.plot('dust_properties.png')

	
		#
		# Grid setup:
		#
		grid_wmin =  0
		grid_wmax =  5.0*pc # 4.0*pc
		grid_zmin =  0.0*pc
		grid_zmax = 10.0*pc
		grid_pmin =  0
		grid_pmax =  2*pi

		grid_dx = cli.resolution*pc
		grid_dw = grid_dx # uniform resolution
		grid_dz = grid_dx # uniform resolution
		grid_dp = grid_dx # resolution at filament location at r = 1 pc

		grid_Nw   = int((grid_wmax - grid_wmin) / grid_dw)
		grid_Nz   = int((grid_zmax - grid_zmin) / grid_dz)
		grid_Np   = int(2*pi * 1.0*pc / grid_dp)

		if(cli.verbose):
			print("Grid setup:")
			print(" Grid resolution =",cli.resolution, "pc.")
			print(" grid_Nw =",grid_Nw)
			print(" grid_Nz =",grid_Nz)
			print(" grid_Np =",grid_Np)

		#grid_w      = np.logspace(np.log10(grid_wmin), np.log10(grid_wmax), grid_Nw)
		#grid_w      = np.hstack([0., grid_w]) # add innermost cell interface at w=0
		grid_w    = np.linspace(grid_wmin, grid_wmax, grid_Nw+1)
		grid_z    = np.linspace(grid_zmin, grid_zmax, grid_Nz+1)
		grid_p    = np.linspace(grid_pmin, grid_pmax, grid_Np+1)

		model.set_cylindrical_polar_grid(grid_w, grid_z, grid_p)

		#
		# Dust density setup:
		#
		RC  = 0.1*pc
		nC  = 6.6580e+03       # in cm^-3
		nC *= cli.opticaldepth # the optical depth at 1 micron
		nC *= m_h              # in g cm^-3
		nC /= 100.0            # converts from gas to dust density
	
		rho = np.zeros(model.grid.shape)
	
		#
		# n(r) = nC / [ 1.0 + (r/RC)**2.0 ]
		# x = -sin(2.0×pi×t) pc, y = +cos(2.0×pi×t) pc, z = 10.0×t pc, t = [0.0, 1.0]
		#  => t = m.grid.gz / (10*pc)
		#  => phi(t) = mod(360*t+270, 360)
		#
		for k in range(0, grid_Np):
			for j in range(0, grid_Nz):
				for i in range(0, grid_Nw):
				
					t = model.grid.gz[k,j,i] / (10*pc)
				
					if(cli.filament == "linear"):
						filament_center_x  = 0
						filament_center_y  = 0
					elif(cli.filament == "spiraling"):
						filament_center_x  = - math.sin(2*pi*t)*pc
						filament_center_y  = + math.cos(2*pi*t)*pc
				
					spherical_grid_r   = model.grid.gw[k,j,i]
					spherical_grid_phi = model.grid.gp[k,j,i]
				
					cartesian_grid_x   = spherical_grid_r * math.cos(spherical_grid_phi)
					cartesian_grid_y   = spherical_grid_r * math.sin(spherical_grid_phi)
				
					rsquared = (
								(cartesian_grid_x - filament_center_x)**2
								+
								(cartesian_grid_y - filament_center_y)**2
								)
				
					rho[k,j,i] = nC / (1.0 + (rsquared / (RC*RC)))
				
					if rsquared**0.5 > 3*pc:
						rho[k,j,i] = 0

		rho[model.grid.gw > grid_wmax] = 0
		rho[model.grid.gz < grid_zmin] = 0
		rho[model.grid.gz > grid_zmax] = 0

		model.add_density_grid(rho, 'dust_properties.hdf5')


		#
		# Check optical depth through the filament:
		#
		#  (y,z = 0, 2.5 pc goes through the filament center in all setups)
		
		#
		# Determine index of closest grid cell to z = 2.5 pc:
		#
		dz_last = 2*abs(grid_zmax-grid_zmin)
		for j in range(0, grid_Nz):
			dz = abs(model.grid.gz[0,j,0] - 2.5*pc)
			if(dz > dz_last):
				j=j-1
				break
			else:
				dz_last = dz

		#
		# Opacity at 1.0 micron (per gram dust):
		#
		chi = dust_properties.optical_properties.interp_chi_wav(1.0)

		tau_max = 0
		for k in range(0, grid_Np):
			tau = 0
			for i in range(0, grid_Nw):
				dr = model.grid.widths[0,k,j,i]
				dtau = dr * rho[k,j,i] * chi
				tau += dtau
			tau_max = max(tau_max, tau)

		if(cli.filament == "linear"):
			tau_max *= 2

		dev = 100 * abs(cli.opticaldepth - tau_max) / cli.opticaldepth

		if(cli.verbose):
			print("Check:")
			print(" Numerical integration of the optical depth through the filament center yields tau = ", tau_max)
			print(" This corresponds to a deviation to the chosen setup value of", dev, "percent")


		#
		# Source:
		#
		if(cli.sources == "external"):
		
			nu, jnu            = np.loadtxt('bg_intensity_modified.txt', unpack=True)
			source_R           = 5*pc
			source             = model.add_external_spherical_source()
			source.peeloff     = False
			source.position    = (0, 0, 5.0*pc) # in a Cartesian frame
			source.radius      = source_R
			source.spectrum    = (nu, jnu)
			#source_MeanIntensity_J = <integrate bg_intensity.txt>
			#source_Area        = 4.0 * pi * source_R*source_R
			source.luminosity  = 8237.0*lsun #source_Area * pi * source_MeanIntensity_J
		
		elif(cli.sources == "stellar"):

			source             = model.add_point_source()
			source.luminosity  = 3.839e35 # in ergs s^-1
			source.temperature = 10000.0 # in K
			if(cli.filament == "linear"):
				source.position    = (3.0*pc, 0, 5.0*pc)
			elif(cli.filament == "spiraling"):
				source.position    = (0     , 0, 3.0*pc)

		#
		# To compute total photon numbers:
		#
		grid_N = grid_Nw * grid_Nz * grid_Np
		if(cli.verbose):
			print("Radiation setup:")
			print(" photons_temperature / cell =", cli.photons_temperature)
			print(" photons_temperature total  =", grid_N * cli.photons_temperature)

		file = filename(cli, "temperature")
		file += ".rtin"

	else:
		file = filename(cli, "temperature")
		file += ".rtout"
	
		try:
			with open(file):
				if(cli.verbose):
					print("Using the specific energy distribution from file", file)
				model.use_geometry(file)
				model.use_quantities(file, only_initial=False, copy=False)
				model.use_sources(file)

		except IOError:
			print("ERROR: File '", file, "' cannot be found. \nERROR: This file, containing the specific energy density, has to be computed first via calling hyperion.")
			exit(2)

		#
		# To compute total photon numbers:
		#
		grid_Nw = len(model.grid.gw[0,0,:])
		grid_Nz = len(model.grid.gw[0,:,0])
		grid_Np = len(model.grid.gw[:,0,0])
		grid_N = grid_Nw * grid_Nz * grid_Np
		if(cli.verbose):
			print("Grid setup:")
			print(" grid_Nw =",grid_Nw)
			print(" grid_Nz =",grid_Nz)
			print(" grid_Np =",grid_Np)
			print("Radiation setup:")
			print(" photons_temperature / cell =", cli.photons_temperature)
			print(" photons_temperature total  =", grid_N * cli.photons_temperature)
			print(" photons_raytracing / cell  =", cli.photons_raytracing)
			print(" photons_raytracing total   =", grid_N * cli.photons_raytracing)
			print(" photons_imaging / cell     =", cli.photons_imaging)
			print(" photons_imaging total      =", grid_N * cli.photons_imaging)

		file = filename(cli, "")
		file += ".rtin"


	##
	## Temperature, Images, and SEDs:
	##
	if(cli.mode == "temperature"):

		model.set_raytracing(True)
		model.set_n_photons(
						initial            = grid_N * cli.photons_temperature,
						raytracing_sources = grid_N * cli.photons_raytracing,
						raytracing_dust    = grid_N * cli.photons_raytracing,
						imaging            = grid_N * cli.photons_imaging
						)
	
	elif(cli.mode == "images"):
	
		model.set_n_initial_iterations(0)
		model.set_raytracing(True)
		model.set_monochromatic(True, wavelengths=[100.0, 500.0, 0.55, 2.2])
		model.set_n_photons(
						raytracing_sources = grid_N * cli.photons_raytracing,
						raytracing_dust    = grid_N * cli.photons_raytracing,
						imaging_sources    = grid_N * cli.photons_imaging,
						imaging_dust       = grid_N * cli.photons_imaging
						)
	
		# group = 0
		image1x = model.add_peeled_images(sed=False, image=True)
		image1x.set_image_size(300, 300)
		image1x.set_image_limits(-5*pc, +5*pc, 0, 10*pc)
		image1x.set_viewing_angles([90],[0]) # along the x-direction
		image1x.set_uncertainties(True)
		image1x.set_output_bytes(8)
		image1x.set_track_origin('basic')
	
		# group = 1
		image1y = model.add_peeled_images(sed=False, image=True)
		image1y.set_image_size(300, 300)
		image1y.set_image_limits(-5*pc, +5*pc, 0, 10*pc)
		image1y.set_viewing_angles([90],[90]) # along the y-direction
		image1y.set_uncertainties(True)
		image1y.set_output_bytes(8)
		image1y.set_track_origin('basic')
	
		# group = 2
		image1z = model.add_peeled_images(sed=False, image=True)
		image1z.set_image_size(300, 300)
		image1z.set_image_limits(-5*pc, +5*pc, -5*pc, +5*pc)
		image1z.set_viewing_angles([0],[0]) # along the z-direction
		image1z.set_uncertainties(True)
		image1z.set_output_bytes(8)
		image1z.set_track_origin('basic')

	elif(cli.mode == "sed"):
	
		model.set_n_initial_iterations(0)
		model.set_raytracing(True)
		model.set_n_photons(
							raytracing_sources = grid_N * cli.photons_raytracing,
							raytracing_dust    = grid_N * cli.photons_raytracing,
							imaging            = grid_N * cli.photons_imaging
							)
	
		# group = 0
		sed1 = model.add_peeled_images(sed=True, image=False)
		sed1.set_wavelength_range(250, 0.01, 2000.0)
		sed1.set_viewing_angles([90],[0]) # along the x-direction
		sed1.set_peeloff_origin((0, 0, 2.5*pc))
		sed1.set_aperture_range(1, 0.3*pc, 0.3*pc)
		sed1.set_uncertainties(True)
		sed1.set_output_bytes(8)
		sed1.set_track_origin('basic')

		# group = 1
		sed2 = model.add_peeled_images(sed=True, image=False)
		sed2.set_wavelength_range(250, 0.01, 2000.0)
		sed2.set_viewing_angles([90],[0]) # along the x-direction
		sed2.set_peeloff_origin((0, 0, 5.0*pc))
		sed2.set_aperture_range(1, 0.3*pc, 0.3*pc)
		sed2.set_uncertainties(True)
		sed2.set_output_bytes(8)
		sed2.set_track_origin('basic')

		# group = 2
		sed3 = model.add_peeled_images(sed=True, image=False)
		sed3.set_wavelength_range(250, 0.01, 2000.0)
		sed3.set_viewing_angles([90],[0]) # along the x-direction
		sed3.set_peeloff_origin((0, 0, 7.5*pc))
		sed3.set_aperture_range(1, 0.3*pc, 0.3*pc)
		sed3.set_uncertainties(True)
		sed3.set_output_bytes(8)
		sed3.set_track_origin('basic')

	##
	## Write model for hyperion runs:
	##
	model.conf.output.output_density         = 'last'
	model.conf.output.output_specific_energy = 'last'
	model.conf.output.output_n_photons       = 'last'
	model.write(file)
	if(cli.verbose):
		print("The input file for hyperion was written to", file)



if(__name__ == "__main__"):
	cli = command_line_arguments()
	handle_command_line_arguments(cli);
	setup_model(cli);

