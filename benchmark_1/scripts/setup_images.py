import os
import glob

from hyperion.model import Model
from hyperion.dust import SphericalDust
from hyperion.util.constants import pc

NPHOTONS = 1e9

for model_path in glob.glob(os.path.join('models', '*_temperature.rtout')):

    m = Model.read(model_path, only_initial=False)

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
    m.set_n_photons(imaging_sources=NPHOTONS, imaging_dust=NPHOTONS,
                    raytracing_sources=1, raytracing_dust=NPHOTONS)

    # Write out and run
    model_name = os.path.basename(model_path).rsplit('.')[0].replace('temperature', 'images')
    m.write(model_name + '.rtin', overwrite=True)
    m.run(model_name + '.rtout', overwrite=True)
