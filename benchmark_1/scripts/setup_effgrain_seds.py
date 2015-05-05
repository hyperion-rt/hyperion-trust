import numpy as np
from hyperion.model import Model
from hyperion.dust import SphericalDust
from hyperion.util.constants import pc

for tau_v in [0.1, 1.0, 20.0]:

    m = Model.read('bm1_slab_eff_tau{0:05.2f}_temperature.rtout'.format(tau_v), only_initial=False)

    m.set_n_initial_iterations(0)

    del m.n_photons['initial']
    del m.n_photons['last']

    i = m.add_peeled_images(sed=True, image=False)
    i.set_viewing_angles([0., 30., 60., 90., 120., 150., 180.], [0., 0., 0., 0., 0., 0., 0.])

    # Set up monochromatic mode
    wavelengths = np.loadtxt('wave_grid_bm1_res5.dat')
    m.set_monochromatic(True, wavelengths=wavelengths)

    # Use raytracing
    m.set_raytracing(True)

    # Set up number of photons
    m.set_n_photons(imaging_sources=1e6, imaging_dust=1e6,
                    raytracing_sources=1, raytracing_dust=1e6)

    # Write out and run
    m.write('bm1_slab_effgrain_tau_{0:05.2f}_seds.rtin'.format(tau_v), overwrite=True)
