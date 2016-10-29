import os
import glob

import numpy as np

from hyperion.model import Model
from hyperion.dust import SphericalDust
from hyperion.util.constants import pc

NPHOTONS = 1e8
NPHOTONS_RAY = 1e8

WAV = np.logspace(-1, 3, 45)
WAV_IM = WAV[np.array([2, 8, 28, 35])]

for model_path in glob.glob(os.path.join('models', '*_temperature.rtout')):

    m = Model.read(model_path, only_initial=False)

    m.set_n_initial_iterations(0)

    del m.n_photons['initial']
    del m.n_photons['last']

    i = m.add_peeled_images()
    i.set_viewing_angles([0., 30., 60., 90., 120., 150., 180.], np.repeat(270, 7))
    i.set_image_limits(-7.5 * pc, 7.5 * pc, -7.5 * pc, 7.5 * pc)
    i.set_track_origin('basic')

    if 't1e+1' in model_path:
        size = 600
    else:
        size = 300

    i.set_image_size(size, size)

    # Set up monochromatic mode
    m.set_monochromatic(True, wavelengths=WAV_IM, energy_threshold=1e-120)

    # Use raytracing
    m.set_raytracing(True)

    # Set up number of photons
    m.set_n_photons(imaging_sources=NPHOTONS, imaging_dust=NPHOTONS,
                    raytracing_sources=1, raytracing_dust=NPHOTONS_RAY)

    # Set up forced first scattering, which will ensure that we get the correct
    # scattered flux at high optical depths in the slab.
    m.set_forced_first_interaction(True, algorithm='baes16', baes16_xi=0.2)

    # When we read in the model above, the link to the dust file is lost, so we
    # replace the dust object by the filename to avoid making the size of
    # input/output files too large.
    m.dust = ['../dust/effgrain_1.0.hdf5']

    # Don't copy input into output
    m.set_copy_input(False)

    # Write out and run
    model_name = os.path.join('models', os.path.basename(model_path).replace('temperature.rtout', 'images'))
    m.write(model_name + '.rtin', overwrite=True, copy=False)
    #m.run(model_name + '.rtout', overwrite=True, mpi=True)
