import os
import glob

import numpy as np
from astropy.io import fits

from hyperion.model import ModelOutput
from hyperion.util.constants import kpc

if not os.path.exists('seds'):
    os.mkdir('seds')

for model_path in glob.glob(os.path.join('models', '*_seds.rtout')):

    m = ModelOutput(model_path)

    model_name = os.path.basename(model_path).replace('_seds.rtout', '')

    for iincl, theta in enumerate([0, 30, 60, 90, 120, 150, 180]):

        sed = m.get_sed(inclination=iincl, units='Jy', distance=10. * kpc, aperture=-1)

        output_file = 'seds/{name}_theta_{theta:03d}_sed.dat'.format(name=model_name, theta=theta)

        np.savetxt(output_file, list(zip(sed.wav, sed.val)), fmt="%10.4e")
