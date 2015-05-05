import os
import glob

from astropy.io import fits

from hyperion.model import ModelOutput
from hyperion.util.constants import kpc

if not os.path.exists('images'):
    os.mkdir('images')

for model_path in glob.glob(os.path.join('models', '*_images.rtout')):

    m = ModelOutput(model_path)

    model_name = os.path.basename(model_path).replace('_images.rtout', '')

    for iincl, theta in enumerate([0, 30, 60, 90, 120, 150, 180]):

        image = m.get_image(inclination=iincl, units='MJy/sr', distance=10. * kpc)

        for iwav, wav in enumerate([0.165, 0.570, 21.3, 161.6]):

            output_file = 'images/{name}_theta_{theta:03d}_wave_{wav:07.3f}.fits'.format(name=model_name, theta=theta, wav=wav)

            fits.writeto(output_file, image.val[:, :, iwav], clobber=True)
