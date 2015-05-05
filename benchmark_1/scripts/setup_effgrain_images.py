import numpy as np
from hyperion.model import Model
from hyperion.dust import SphericalDust
from hyperion.util.constants import pc

for tau_v in [0.1, 1.0, 20.0]:

    m = Model.read('bm1_slab_eff_tau{0:05.2f}_temperature.rtout'.format(tau_v), only_initial=False)

    m.set_n_initial_iterations(0)

    del m.n_photons['initial']
    del m.n_photons['last']

    i = m.add_peeled_images()
    i.set_viewing_angles([0., 30., 60., 90., 120., 150., 180.], [0., 0., 0., 0., 0., 0., 0.])
    i.set_image_limits(-7.5 * pc, 7.5 * pc, -7.5 * pc, 7.5 * pc)
    i.set_image_size(300, 300)

    # Set up monochromatic mode
    m.set_monochromatic(True, wavelengths=[0.165, 0.570, 21.3, 161.6])

    # Use raytracing
    m.set_raytracing(True)

    # Set up number of photons
    m.set_n_photons(imaging_sources=1e9, imaging_dust=1e9,
                    raytracing_sources=1, raytracing_dust=1e9)

    # Write out and run
    m.write('bm1_slab_effgrain_tau_{0:05.2f}_images.rtin'.format(tau_v), overwrite=True)
