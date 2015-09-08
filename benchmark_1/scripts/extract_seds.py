import os
import glob

import numpy as np
from astropy.io import fits

from hyperion.model import ModelOutput
from hyperion.util.constants import kpc

import matplotlib.pyplot as plt

if not os.path.exists('seds'):
    os.mkdir('seds')

for model_path in glob.glob(os.path.join('models', '*_seds.rtout')):

    m = ModelOutput(model_path)

    model_name = os.path.basename(model_path).replace('_seds.rtout', '')

    for iincl, theta in enumerate([0, 30, 60, 90, 120, 150, 180]):

        sed_total = m.get_sed(inclination=iincl, units='Jy', distance=10. * kpc, aperture=-1, group=0, component='total')
        sed_semit = m.get_sed(inclination=iincl, units='Jy', distance=10. * kpc, aperture=-1, group=0, component='source_emit')
        sed_sscat = m.get_sed(inclination=iincl, units='Jy', distance=10. * kpc, aperture=-1, group=0, component='source_scat')
        sed_demit = m.get_sed(inclination=iincl, units='Jy', distance=10. * kpc, aperture=-1, group=0, component='dust_emit')
        sed_dscat = m.get_sed(inclination=iincl, units='Jy', distance=10. * kpc, aperture=-1, group=0, component='dust_scat')

        sed_trans = m.get_sed(inclination=iincl, units='Jy', distance=10. * kpc, aperture=-1, group=1, component='source_emit')

        output_file = 'seds/{name}_i{theta:03d}a000.sed'.format(name=model_name, theta=theta)

        # TODO: remove me once fixed upstream
        output_file = output_file.replace('t1e+0', 't1e0').replace('t1e+1', 't1e1').replace('t1e+2', 't1e2')

        np.savetxt(output_file, list(zip(sed_total.wav,
        								 sed_total.val,
        								 sed_semit.val,
        								 sed_sscat.val,
        								 sed_demit.val,
        								 sed_dscat.val,
                                         sed_trans.val)), fmt="%12.5e")

        # fig = plt.figure()
        # ax = fig.add_subplot(1,1,1)
        # ax.plot(sed_total.wav, sed_total.val, color='k', label='total')
        # ax.plot(sed_total.wav, sed_semit.val, color='blue', label='source (emit)')
        # ax.plot(sed_total.wav, sed_sscat.val, color='purple', label='source (scat)')
        # ax.plot(sed_total.wav, sed_demit.val, color='red', label='dust (emit)')
        # ax.plot(sed_total.wav, sed_dscat.val, color='orange', label='dust (scat)')
        # ax.plot(sed_total.wav, sed_trans.val, color='0.5', label='transparent')
        # ax.set_xscale('log')
        # ax.set_yscale('log')
        # ax.set_ylim(0.000001, 1000)
        # ax.legend(fontsize=7, loc=2)
        # fig.savefig(output_file + '.png')
        # plt.close(fig)