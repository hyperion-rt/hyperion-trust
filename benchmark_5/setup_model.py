#############################################
## Setup for TRUST-V Galaxy benchmark test ##
#############################################

import sys
import numpy as np
import math
from hyperion.model import Model
from hyperion.dust  import SphericalDust
from hyperion.util.constants import pc, pi, m_h, msun
from cli import command_line_arguments, handle_command_line_arguments, filename


def setup_model(cli):
	
    lsun_TRUST = 3.839e33
        
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
        # Specify galaxy setup:
        #
        hR                     =  4000.0*pc             # [cm]
        Rmax                   =     5.0*hR             # [cm]
        hz_oldstars            =   350.0*pc             # [cm]
        hz_youngstars          =   200.0*pc             # [cm]
        hz_dust                =   200.0*pc             # [cm]
        zmax_oldstars          =     5.0*hz_oldstars    # [cm]
        zmax_youngstars        =     5.0*hz_youngstars  # [cm]
        zmax_dust              =     5.0*hz_dust        # [cm]
        zmax                   =  zmax_oldstars         # [cm]
        reff                   =  1600.0*pc             # [cm]
        n                      =     3.0
        q                      =     0.6
        bn                     = 2.0*n - 1.0/3.0 + 4.0/405.0/n + 46.0/25515.0/n/n + 131.0/1148175.0/n/n/n
        temperature_oldstars   =  3500.0                # [K]
        temperature_youngstars = 10000.0                # [K]
        temperature_bulge      =  3500.0                # [K]
        luminosity_oldstars    =     4.0e+10*lsun_TRUST # [ergs/s]
        luminosity_youngstars  =     1.0e+10*lsun_TRUST # [ergs/s]
        luminosity_bulge       =     3.0e+10*lsun_TRUST # [ergs/s]
        
        w_oldstars             =     0.25
        w_youngstars           =     0.75
        w_dust                 =     0.75
        phi0_oldstars          =     0.0
        phi0_youngstars        =    20.0 * pi/180.0
        phi0_dust              =    20.0 * pi/180.0
        modes                  =     2
        pitchangle             =    20.0 * pi/180.0
        
        
        
        #
        # Grid setup:
        #
        grid_wmin =  0.0
        grid_wmax =  Rmax
        grid_zmin = -zmax
        grid_zmax = +zmax
        grid_pmin =  0.0
        grid_pmax =  2.0*pi
        
        grid_dx = cli.resolution*pc
        grid_dw = grid_dx # uniform resolution
        grid_dz = grid_dx # uniform resolution
        grid_dp = grid_dx # resolution at characteristic radial disk spatial scale hR = 4000.0 pc
        
        grid_Nw   = int((grid_wmax - grid_wmin) / grid_dw) + 1
        grid_Nz   = int((grid_zmax - grid_zmin) / grid_dz) + 1
        if(cli.case == 1):
            grid_Np = 1
        if(cli.case == 2):
            grid_Np = int((grid_pmax - grid_pmin) * hR / grid_dp)
        
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
        # Dust density and sources setup:
        #
        rho_oldstars   = np.zeros(model.grid.shape)
        rho_youngstars = np.zeros(model.grid.shape)
        rho_bulge      = np.zeros(model.grid.shape)
        rho_dust       = np.zeros(model.grid.shape)
        
        for k in range(0, grid_Np):
            for j in range(0, grid_Nz):
                for i in range(0, grid_Nw):
                    
                    R = model.grid.gw[k,j,i]
                    z = model.grid.gz[k,j,i]
                    m = math.sqrt(R*R + z*z/q/q)
                    
                    rho_dust[k,j,i]       = math.exp(- R/hR -abs(z)/hz_dust      )
                    rho_oldstars[k,j,i]   = math.exp(- R/hR -abs(z)/hz_oldstars  )
                    rho_youngstars[k,j,i] = math.exp(- R/hR -abs(z)/hz_youngstars)
                    rho_bulge[k,j,i]      = math.pow(m/reff, 0.5/n - 1.0) * math.exp(- bn * math.pow(m/reff, 1.0/n))
                    
                    if(cli.case == 2):
                        phi = model.grid.gp[k,j,i]
                        perturb = math.sin(modes * (math.log(R/hR) / math.tan(pitchangle) - (phi - phi0_dust)))
                        rho_dust[k,j,i]       *= (1.0 + w_dust       * perturb)
                        perturb = math.sin(modes * (math.log(R/hR) / math.tan(pitchangle) - (phi - phi0_oldstars)))
                        rho_oldstars[k,j,i]   *= (1.0 + w_oldstars   * perturb)
                        perturb = math.sin(modes * (math.log(R/hR) / math.tan(pitchangle) - (phi - phi0_youngstars)))
                        rho_youngstars[k,j,i] *= (1.0 + w_youngstars * perturb)
        
        rho_dust[model.grid.gw > grid_wmax] = 0
        rho_dust[model.grid.gz < grid_zmin] = 0
        rho_dust[model.grid.gz > grid_zmax] = 0
        
        kappa_ref     = dust_properties.optical_properties.interp_chi_wav(0.55693)
        rho0          = cli.opticaldepth / (2.0 * hz_dust * kappa_ref)
        rho_dust[:]  *= rho0
        model.add_density_grid(rho_dust, 'dust_properties.hdf5')
        
        source_oldstars                = model.add_map_source()
        source_oldstars.luminosity     = luminosity_oldstars
        source_oldstars.temperature    = temperature_oldstars
        source_oldstars.map            = rho_oldstars
        
        source_youngstars              = model.add_map_source()
        source_youngstars.luminosity   = luminosity_youngstars
        source_youngstars.temperature  = temperature_youngstars
        source_youngstars.map          = rho_youngstars
        
        source_bulge                   = model.add_map_source()
        source_bulge.luminosity        = luminosity_bulge
        source_bulge.temperature       = temperature_bulge
        source_bulge.map               = rho_bulge
        
        
        #
        # Check face-on optical depth at 1.0 micron (per gram dust) through the dust disk:
        #
        tau   = 0
        
        k = 0
        i = 0
        for j in range(0, grid_Nz):
            #print(model.grid.gz[k,j,i]/pc, rho_dust[k,j,i])
            dz   = model.grid.widths[1,k,j,i]
            dtau = dz * rho_dust[k,j,i] * kappa_ref
            tau += dtau
        
        deviation = 100.0 * abs(cli.opticaldepth - tau) / cli.opticaldepth
        
        if(cli.verbose):
            print("Check optical depth of dust density setup:")
            print(" kappa(0.55693 micron) = ", kappa_ref, "cm^2 g^-1")
            print(" Numerical integration of the face-on optical depth at 0.55693 micron through the central dust disk yields tau = ", tau)
            print(" This corresponds to a deviation to the chosen setup value of", deviation, "percent")
    
        #
        # Check central dust density:
        #
        rho_max = np.max(rho_dust)
        if(cli.opticaldepth < 1.0):
            rho_setup = 1.04366e-4 * msun/pc/pc/pc
        if(cli.opticaldepth < 3.0):
            rho_setup = 5.21829e-4 * msun/pc/pc/pc
        else:
            rho_setup = 2.60915e-3 * msun/pc/pc/pc

        deviation = 100.0 * abs(rho_setup - rho_max) / rho_setup

        if(cli.verbose):
            print("Check value of central dust density:")
            print(" rho_max = ", rho_max, "g cm^-3")
            print(" This corresponds to a deviation to the chosen setup value of", deviation, "percent")

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
        # old setup: model.set_monochromatic(True, wavelengths=[0.4, 1.0, 10.0, 100.0, 500.0])
        model.set_monochromatic(True, wavelengths=[0.45483, 1.2520, 26.114, 242.29])
        model.set_n_photons(
            raytracing_sources = grid_N * cli.photons_raytracing,
            raytracing_dust    = grid_N * cli.photons_raytracing,
            imaging_sources    = grid_N * cli.photons_imaging,
            imaging_dust       = grid_N * cli.photons_imaging
        )
    
        # group = 0
        image1 = model.add_peeled_images(sed=False, image=True)
        image1.set_image_size(501, 501)
        image1.set_image_limits(-12500.0*pc, +12500.0*pc, -12500.0*pc, +12500.0*pc)
        image1.set_viewing_angles([30],[0])
        image1.set_uncertainties(True)
        image1.set_output_bytes(8)
        image1.set_track_origin('basic')
    
        # group = 1
        image2 = model.add_peeled_images(sed=False, image=True)
        image2.set_image_size(501, 501)
        image2.set_image_limits(-12500.0*pc, +12500.0*pc, -12500.0*pc, +12500.0*pc)
        image2.set_viewing_angles([80],[90])
        image2.set_uncertainties(True)
        image2.set_output_bytes(8)
        image2.set_track_origin('basic')
    
        # group = 2
        image3 = model.add_peeled_images(sed=False, image=True)
        image3.set_image_size(501, 501)
        image3.set_image_limits(-12500.0*pc, +12500.0*pc, -12500.0*pc, +12500.0*pc)
        image3.set_viewing_angles([88],[0]) # mostly edge-on
        image3.set_uncertainties(True)
        image3.set_output_bytes(8)
        image3.set_track_origin('basic')

    elif(cli.mode == "seds"):
        
        model.set_n_initial_iterations(0)
        model.set_raytracing(True)
        model.set_n_photons(
            raytracing_sources = grid_N * cli.photons_raytracing,
            raytracing_dust    = grid_N * cli.photons_raytracing,
            imaging            = grid_N * cli.photons_imaging
        )
    
        # group = 0
        sed1 = model.add_peeled_images(sed=True, image=False)
        sed1.set_wavelength_range(47, 0.081333, 1106.56)
        sed1.set_viewing_angles([30],[0])
        sed1.set_peeloff_origin((0, 0, 0))
        sed1.set_aperture_range(1, 25000.0*pc, 25000.0*pc)
        sed1.set_uncertainties(True)
        sed1.set_output_bytes(8)
        sed1.set_track_origin('basic')
        
        # group = 1
        sed2 = model.add_peeled_images(sed=True, image=False)
        sed2.set_wavelength_range(47, 0.081333, 1106.56)
        sed2.set_viewing_angles([80],[0])
        sed2.set_peeloff_origin((0, 0, 0))
        sed2.set_aperture_range(1, 25000.0*pc, 25000.0*pc)
        sed2.set_uncertainties(True)
        sed2.set_output_bytes(8)
        sed2.set_track_origin('basic')
    
        # group = 2
        sed3 = model.add_peeled_images(sed=True, image=False)
        sed3.set_wavelength_range(47, 0.081333, 1106.56)
        sed3.set_viewing_angles([88],[0])
        sed3.set_peeloff_origin((0, 0, 0))
        sed3.set_aperture_range(1, 25000.0*pc, 25000.0*pc)
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

