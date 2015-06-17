import os
import glob

import numpy as np

from astropy.io import fits

from hyperion.model import ModelOutput
from hyperion.util.constants import kpc

WAV = np.logspace(-1, 3, 45)
WAV_IM = WAV[np.array([2, 8, 26, 35])]

if not os.path.exists('images'):
    os.mkdir('images')

for model_path in glob.glob(os.path.join('models', '*_images.rtout')):

    m = ModelOutput(model_path)

    model_name = os.path.basename(model_path).replace('_images.rtout', '')

    for iincl, theta in enumerate([0, 30, 60, 90, 120, 150, 180]):

        image = m.get_image(inclination=iincl, units='MJy/sr', distance=10. * kpc)

        # output_file = 'images/{name}_theta_{theta:03d}.fits'.format(name=model_name, theta=theta)
        # fits.writeto(output_file, np.rot90(image.val).transpose().swapaxes(1,2), clobber=True)
        # fits.writeto(output_file, image.val.transpose()[:,:,::-1], clobber=True)

        for iwav, wav in enumerate(WAV_IM):
        
            output_file = 'images/{name}_i{theta:03d}a000_w{wav:06.2f}.fits'.format(name=model_name, theta=theta, wav=wav)
        
            fits.writeto(output_file, image.val[:, :, iwav], clobber=True)
