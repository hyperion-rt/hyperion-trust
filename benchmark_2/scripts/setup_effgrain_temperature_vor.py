import os

import numpy as np

from hyperion.model import Model
from hyperion.dust import SphericalDust
from hyperion.util.constants import pc, au, sigma, pi, rsun

NPHOTONS = 1e7

if not os.path.exists('models'):
    os.mkdir('models')

n0 = 100000
n1 = 100000
n2 = 100000
n3 = 100000

RSUB = 15 * au

def random_in_sphere(rmin, rmax, n):
    phi = np.random.uniform(0., 2 * np.pi, n)
    mu = np.random.uniform(-1., 1., n)
    sin_theta = (1 - mu*mu)**0.5
    r = (np.random.random(n) * (rmax ** 3 - rmin ** 3) + rmin ** 3) ** (1. / 3.)
    x = r * np.cos(phi) * sin_theta
    y = r * np.sin(phi) * sin_theta
    z = r * mu
    return x, y, z

# TODO: remove dust around source

m = Model()

# Grain Properties:

d = SphericalDust('../dust/integrated_hg_scattering.hdf5')
chi_v = d.optical_properties.interp_chi_wav(0.55)

# Determine density and density in ambient medium
tau_0 = 0.1
size = 60 * au
rho_0 = tau_0 / (size * chi_v)
mass_0 = rho_0 * size ** 3

# Determine density and mass in sphere 1
tau_1 = 1000
r_1 = 5 * au
rho_1 = tau_1 / (2 * r_1 * chi_v)
mass_1 = (rho_1 - rho_0) * 4./3. * pi * r_1 ** 3

# Determine density and mass in sphere 2
tau_2 = 100
r_2 = 20 * au
rho_2 = tau_2 / (2 * r_2 * chi_v)
mass_2 = (rho_2 - rho_0) * 4./3. * pi * r_2 ** 3

# Ambient medium
x0 = np.random.uniform(0., 60. * au, n0)
y0 = np.random.uniform(0., 60. * au, n0)
z0 = np.random.uniform(0., 60. * au, n0)

# Sphere 1
x1, y1, z1 = random_in_sphere(r_1 * 0.95, r_1 * 1.05, n1)
x1 += 10 * au
y1 += 15 * au
z1 += 20 * au

# Sphere 2
x2, y2, z2 = random_in_sphere(r_2 * 0.95, r_2 * 1.05, n2)
x2 += 26.666667 * au
y2 += 31.666667 * au
z2 += 28.333333 * au

# Inner edge
x3, y3, z3 = random_in_sphere(RSUB * 0.95, RSUB * 1.05, n3)
x3 = np.abs(x3)
y3 = np.abs(y3)
z3 = np.abs(z3)

# Concatenate
x = np.hstack([x0, x1, x2, x3])
y = np.hstack([y0, y1, y2, y3])
z = np.hstack([z0, z1, z2, z3])

# np.savetxt('points.txt', list(zip(x, y, z)))

m.set_voronoi_grid(x, y, z,
                   xmin=-10 * au, xmax=60 * au,
                   ymin=-10 * au, ymax=60 * au,
                   zmin=-10 * au, zmax=60 * au)

density = np.zeros(len(x))

# Ambient medium
in_box = (x > 0.) & (x < 60 * au) & (y > 0.) & (y < 60 * au) & (z > 0.) & (z < 60 * au)
density[in_box] = rho_0

# Set up sphere 1
in_sphere_1 = (x - 10 * au) ** 2 + (y - 15 * au) ** 2 + (z - 20 * au) ** 2 < r_1 ** 2
density[in_sphere_1] = rho_1

# Set up sphere 2
in_sphere_2 = (x - 26.666667 * au) ** 2 + (y - 31.666667 * au) ** 2 + (z - 28.333333 * au) ** 2 < r_2 ** 2
density[in_sphere_2] = rho_2

# Remove dust close to source
in_rsub = np.sqrt(x * x + y * y + z * z) < RSUB
density[in_rsub] = 0.

m.add_density_grid(density, d)

# m.set_propagation_check_frequency(1.0)

# Set up illuminating source:
s = m.add_spherical_source()
s.radius = 6.6 * rsun
s.temperature = 33000.
s.luminosity = 4 * pi * s.radius ** 2 * sigma * s.temperature ** 4

# Set up number of photons
m.set_n_photons(initial=NPHOTONS, imaging=0)

# Write out and run
m.write(os.path.join('models', 'bm2_eff_vor_temperature.rtin'), overwrite=True)
m.run(os.path.join('models', 'bm2_eff_vor_temperature.rtout'), mpi=True, overwrite=True)
