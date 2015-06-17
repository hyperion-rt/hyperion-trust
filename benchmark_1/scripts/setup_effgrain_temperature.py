import os

import numpy as np

from hyperion.model import Model
from hyperion.dust import SphericalDust
from hyperion.util.constants import pc, c

NPHOTONS = 1e9
NITER_MAX = 20

if not os.path.exists('models'):
    os.mkdir('models')


TAU_LABEL = {}
TAU_LABEL[0.01] = "1e-2"
TAU_LABEL[0.1] = "1e-1"
TAU_LABEL[1] = "1e+0"
TAU_LABEL[10] = "1e+1"
TAU_LABEL[100] = "1e+2"


for tau_v in [0.01, 0.1, 1, 10, 100]:

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
    zfine = -2 - np.logspace(-3, np.log10(3), 200)
    z = np.hstack([zfine[::-1], -2, 5]) * pc

    m.set_cartesian_grid(x, y, z)

    # Grain Properties:

    d = SphericalDust('../dust/integrated_hg_scattering_mid.hdf5')
    chi_v = d.optical_properties.interp_chi_wav(0.55)

    # Determine density in slab
    rho0 = tau_v / (3 * pc * chi_v)

    # Set up density grid
    density = np.ones(m.grid.shape) * rho0
    density[-1,:,:] = 0.

    m.add_density_grid(density, d)

    # Set up illuminating source:

    wav, fnu = np.loadtxt('data/BB_T10000_L100000.dat', usecols=[0,1], unpack=True)
    nu = c / (wav * 1.e-4)
    nu = nu[::-1]
    fnu = fnu[::-1]

    s = m.add_point_source()
    s.position = (0., 0., 4 * pc)
    s.luminosity = 3.839e38
    s.spectrum = (nu, fnu)

    # Set up number of photons
    m.set_n_photons(initial=NPHOTONS, imaging=0)

    m.conf.output.output_specific_energy = 'all'
    m.set_n_initial_iterations(NITER_MAX)

    m.set_convergence(True, percentile=99.9, absolute=2., relative=1.01)

    # Write out and run
    model_name = 'models/hyper_slab_eff_t{0}_temperature'.format(TAU_LABEL[tau_v])
    m.write(model_name + '.rtin', overwrite=True)
    # m.run(model_name + '.rtout', mpi=True, overwrite=True, n_processes=12)
