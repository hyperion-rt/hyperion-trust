from hyperion.model import ModelOutput
from hyperion.util.constants import kpc
from astropy.io import fits


for tau in [0.1, 1.0, 20.]:

    input_file = 'bm1_slab_effgrain_tau_{tau:05.2f}_seds.rtout'.format(tau=tau)

    import h5py
    f = h5py.File(input_file, 'r')
    print f['Peeled'].keys()

    m = ModelOutput(input_file)

    for iincl, theta in enumerate([0, 30, 60, 90, 120, 150, 180]):

        sed = m.get_sed(inclination=iincl, units='Jy', distance=10. * kpc, aperture=-1)

        output_file = 'seds/bm1_slab_effgrain_tau_{tau:05.2f}_theta_{theta:03d}_sed.txt'.format(tau=tau, theta=theta)

        np.savetxt(output_file, zip(sed.wav, sed.val), fmt="%10.4e")
