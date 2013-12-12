from hyperion.model import ModelOutput
from hyperion.util.constants import kpc
from astropy.io import fits


for tau in [0.1, 1.0, 20.]:

	input_file = 'bm1_slab_eff_tau{tau:05.2f}_images.rtout'.format(tau=tau)

	m = ModelOutput(input_file)

	for iincl, theta in enumerate([0, 30, 60, 90, 120, 180]):

		image = m.get_image(inclination=iincl, units='MJy/sr', distance=10. * kpc)

		for iwav, wav in enumerate([0.165, 0.570, 21.3, 161.6]):

			output_file = 'bm1_hyperion/bm1_slab_effgrain_theta_{theta:03d}_wave_{wav:05.1f}.fits'.format(theta=theta, wav=wav)

			fits.writeto(output_file, image.val[:, :, iwav], clobber=True)