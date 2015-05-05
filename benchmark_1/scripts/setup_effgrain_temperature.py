import os

import numpy as np

from hyperion.model import Model
from hyperion.dust import SphericalDust
from hyperion.util.constants import pc

if not os.path.exists('models'):
    os.mkdir('models')


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

    d = SphericalDust('../dust/integrated_hg_scattering.hdf5')
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
    s.luminosity = 3.839e38

    # Set up number of photons
    m.set_n_photons(initial=1e9, imaging=0)

    # Write out and run
    m.write('models/bm1_slab_eff_tau{0:05.2f}_temperature.rtin'.format(tau_v), overwrite=True)
