import os
import glob

import numpy as np

from hyperion.model import Model
from hyperion.dust import SphericalDust
from hyperion.util.constants import pc

NPHOTONS = 1e6

for model_path in glob.glob(os.path.join('models', '*_temperature.rtout')):

    m = Model.read(model_path, only_initial=False)

    m.set_n_initial_iterations(0)

    del m.n_photons['initial']
    del m.n_photons['last']

    i = m.add_peeled_images(sed=True, image=False)
    i.set_viewing_angles([0., 30., 60., 90., 120., 150., 180.],
                         [0., 0., 0., 0., 0., 0., 0.])

    # Set up monochromatic mode
    wavelengths = np.loadtxt('data/wave_grid_bm1_res5.dat')
    m.set_monochromatic(True, wavelengths=wavelengths)

    # Use raytracing
    m.set_raytracing(True)

    # Set up number of photons
    m.set_n_photons(imaging_sources=NPHOTONS, imaging_dust=NPHOTONS,
                    raytracing_sources=1, raytracing_dust=NPHOTONS)

    # Write out and run
    model_name = os.path.basename(model_path).rsplit('.')[0].replace('temperature', 'seds')
    m.write(model_name + '.rtin', overwrite=True)
    m.run(model_name + '.rtout', overwrite=True)
