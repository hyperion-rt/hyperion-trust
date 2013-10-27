import numpy as np
from hyperion.model import Model
from hyperion.dust import SphericalDust
from hyperion.util.constants import pc

for tau_v in [0.1, 1.0, 20.0]:

    m = Model()

    # Global geometry:
    #
    # * slab
    # * system size = 10x10x10 pc
    # * system coordinates (x,y,z min/max) = -5 to +5 pc
    # * slab z extent = -2 to -5 pc
    # * slab xy extend = -5 pc to 5 pc
    # * z optical depth @0.55um in slab = 0.1, 1, 20
    # * optical depth outside slab = 0

    x = np.linspace(-5 * pc, 5 * pc, 100)
    y = np.linspace(-5 * pc, 5 * pc, 100)
    z = np.hstack([np.linspace(-5 * pc, -2 * pc, 100), 5 * pc])

    m.set_cartesian_grid(x, y, z)

    # Grain Properties:

    d = SphericalDust('integrated_hg_scattering.hdf5')
    chi_v = d.optical_properties.interp_chi_wav(0.55)

    # Determine density in slab
    rho0 = tau_v / (3 * pc * chi_v)

    # Set up density grid
    density = np.ones(m.grid.shape) * rho0
    density[-1,:,:] = 0.

    m.add_density_grid(density, d)

    # Set up illuminating source:
    s = m.add_point_source()
    s.position = (0., 0., 4 * pc)
    s.temperature = 10000.0
    s.luminosity = 3.845e38

    # TODO: Specrtum is provided, so we could use that

    i = m.add_peeled_images()
    i.set_viewing_angles([0., 30., 60., 90., 120., 180.], [0., 0., 0., 0., 0., 0.])
    i.set_image_limits(-7.5 * pc, 7.5 * pc, -7.5 * pc, 7.5 * pc)
    i.set_image_size(300, 300)

    # Read in wavelength grid
    wav = np.loadtxt('wave_grid_bm1_res5.dat')

    # Set up monochromatic mode
    m.set_monochromatic(True, wavelengths=wav)

    # Use raytracing
    m.set_raytracing(True)

    # Set up number of photons
    m.set_n_photons(initial=1e8, imaging_sources=1e8, imaging_dust=1e8,
                    raytracing_sources=1e8, raytracing_dust=1e8)

    # Write out and run
    m.write('bm1_slab_eff_tau{0:05.2f}.rtin'.format(tau_v), overwrite=True)
    m.run('bm1_slab_eff_tau{0:05.2f}.rtout'.format(tau_v), mpi=True, n_processes=4, overwrite=True)
