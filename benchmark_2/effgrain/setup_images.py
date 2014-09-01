import numpy as np
from hyperion.model import Model
from hyperion.dust import SphericalDust
from hyperion.util.constants import au

m = Model.read('bm2_eff_vor_temperature.rtout', only_initial=False)

m.set_n_initial_iterations(0)

del m.n_photons['initial']
del m.n_photons['last']

i = m.add_peeled_images()
i.set_viewing_angles([0., 90., 90., 90., 90., 180.], [0., 0., 90., 180., 270., 0.])
i.set_image_limits(-60 * au, 60 * au, -60 * au, 60 * au)
i.set_image_size(300, 300)

# Set up monochromatic mode
m.set_monochromatic(True, wavelengths=[0.10019, 0.55165, 2.00293, 10.03850, 101.15800])

# Use raytracing
m.set_raytracing(True)

# Set up number of photons
m.set_n_photons(imaging_sources=1e7, imaging_dust=1e7,
                raytracing_sources=1, raytracing_dust=1e7)

# Write out and run
m.write('bm2_eff_images.rtin', overwrite=True)
m.run('bm2_eff_images.rtout', mpi=True)
