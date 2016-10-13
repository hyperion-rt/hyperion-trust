import os
import glob

import numpy as np

from hyperion.model import Model
from hyperion.dust import SphericalDust
from hyperion.util.constants import pc

NPHOTONS = 1e6
NPHOTONS_RAY = 1e8

WAV = np.logspace(-1, 3, 45)

for model_path in glob.glob(os.path.join('models', '*_temperature.rtout')):

    m = Model.read(model_path, only_initial=False)

    m.set_n_initial_iterations(0)

    del m.n_photons['initial']
    del m.n_photons['last']

    i = m.add_peeled_images(sed=True, image=False)
    i.set_viewing_angles([0., 30., 60., 90., 120., 150., 180.],
                         [0., 0., 0., 0., 0., 0., 0.])
    i.set_track_origin('basic')

    i = m.add_peeled_images(sed=True, image=False)
    i.set_viewing_angles([0., 30., 60., 90., 120., 150., 180.],
                         [0., 0., 0., 0., 0., 0., 0.])
    i.set_track_origin('basic')
    i.set_ignore_optical_depth(True)

    # Set up monochromatic mode
    m.set_monochromatic(True, wavelengths=WAV, energy_threshold=1e-120)

    # Use raytracing
    m.set_raytracing(True)

    # Set up number of photons
    m.set_n_photons(imaging_sources=NPHOTONS, imaging_dust=NPHOTONS,
                    raytracing_sources=1, raytracing_dust=NPHOTONS_RAY)

    m.set_copy_input(False)
    m.conf.output.output_specific_energy = 'none'

    # Set up forced first scattering, which will ensure that we get the correct
    # scattered flux at high optical depths in the slab.
    m.set_forced_first_interaction(True, algorithm='baes16', baes16_xi=0.2)

    # Write out and run
    model_name = os.path.join('models', os.path.basename(model_path).replace('temperature.rtout', 'seds'))
    m.write(model_name + '.rtin', overwrite=True)
    #m.run(model_name + '.rtout', mpi=True, overwrite=True, n_processes=4)

