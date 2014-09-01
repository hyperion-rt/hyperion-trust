import numpy as np
from hyperion.model import Model
from hyperion.dust import SphericalDust
from hyperion.util.constants import pc, au, sigma, pi, rsun

# TODO: remove dust around source

m = Model()

x = np.linspace(0., 60. * au, 256)
y = np.linspace(0., 60. * au, 256)
z = np.linspace(0., 60. * au, 256)

x = np.hstack([-10 * au, x])
y = np.hstack([-10 * au, y])
z = np.hstack([-10 * au, z])

m.set_cartesian_grid(x, y, z)

# Grain Properties:

d = SphericalDust('integrated_hg_scattering.hdf5')
chi_v = d.optical_properties.interp_chi_wav(0.55)

# Determine density in ambient medium
tau_0 = 0.1
size = 60 * au
rho_0 = tau_0 / (size * chi_v)

# Determine density in sphere 1
tau_1 = 1000
r_1 = 5 * au
rho_1 = tau_1 / (2 * r_1 * chi_v)

# Determine density in sphere 2
tau_2 = 100
r_2 = 20 * au
rho_2 = tau_2 / (2 * r_2 * chi_v)

# Set up ambient medium
density = np.ones(m.grid.shape) * rho_0

# Set up sphere 1
in_sphere_1 = (m.grid.gx - 10 * au) ** 2 + (m.grid.gy - 15 * au) ** 2 + (m.grid.gz - 20 * au) ** 2 < r_1 ** 2
density[in_sphere_1] = rho_1

# Set up sphere 2
in_sphere_2 = (m.grid.gx - 26.666667 * au) ** 2 + (m.grid.gy - 31.666667 * au) ** 2 + (m.grid.gz - 28.333333 * au) ** 2 < r_2 ** 2
density[in_sphere_2] = rho_2

# Remove dust outside original cube
density[0,:,:] = 0
density[:,0,:] = 0
density[:,:,0] = 0

m.add_density_grid(density, d)

# Set up illuminating source:
s = m.add_spherical_source()
s.radius = 6.6 * rsun
s.temperature = 33000.
s.luminosity = 4 * pi * s.radius ** 2 * sigma * s.temperature ** 4

# Set up number of photons
m.set_n_photons(initial=1e7, imaging=0)

# Write out and run
m.write('bm2_eff_temperature.rtin', overwrite=True)
m.run('bm2_eff_temperature.rtout', mpi=True, n_processes=12)
