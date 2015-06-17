import os
import glob

import numpy as np

from hyperion.model import Model
from hyperion.dust import SphericalDust
from hyperion.util.constants import pc

NPHOTONS = 1e8
NPHOTONS_RAY = 1e8

WAV = np.logspace(-1, 3, 45)
WAV_IM = WAV[np.array([2, 8, 26, 35])]

for model_path in glob.glob(os.path.join('models', '*_temperature.rtout')):

    m = Model.read(model_path, only_initial=False)

    m.set_n_initial_iterations(0)

    del m.n_photons['initial']
    del m.n_photons['last']

    i = m.add_peeled_images()
    i.set_viewing_angles([0., 30., 60., 90., 120., 150., 180.], np.repeat(270, 7))
    i.set_image_limits(-7.5 * pc, 7.5 * pc, -7.5 * pc, 7.5 * pc)

    if 't1e+2' in model_path:
        size = 3000
    elif 't1e+1' in model_path:
        size = 600
    else:
        size = 300

    i.set_image_size(size, size)

    # Set up monochromatic mode
    m.set_monochromatic(True, wavelengths=WAV_IM)

    # Use raytracing
    m.set_raytracing(True)

    # Set up number of photons
    m.set_n_photons(imaging_sources=NPHOTONS, imaging_dust=NPHOTONS,
                    raytracing_sources=1, raytracing_dust=NPHOTONS_RAY)

    # Write out and run
    model_name = os.path.join('models', os.path.basename(model_path).replace('temperature.rtout', 'images'))
    m.write(model_name + '.rtin', overwrite=True)
    #m.run(model_name + '.rtout', overwrite=True, mpi=True)
