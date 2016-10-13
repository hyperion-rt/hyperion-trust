import os

import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
from hyperion.util.constants import kpc

plt.rcParams['font.family'] = 'Arial'


fig1 = plt.figure(figsize=(10,7))
fig2 = plt.figure(figsize=(10,7))

for i, tau_v in enumerate(['1e-2', '1e-1', '1e0', '1e1']):

    ax1 = fig1.add_subplot(2,2,i+1)
    ax2 = fig2.add_subplot(2,2,i+1)

    for theta in [0, 30, 60, 90, 120, 150, 180]:

        reference = os.path.join('reference', 'dirty', 'dirty_slab_eff_t{0}_i{1:03d}a000.sed'.format(tau_v, theta))
        hyperion = os.path.join('hyper_slab_eff_0.9.2_v4', 'hyper_slab_eff_t{0}_i{1:03d}a000.sed'.format(tau_v, theta))

        wav_ref, fnu_ref = np.loadtxt(reference, usecols=[0,1], unpack=True)

        # Temporarily compare to input spectrum
        #wav_ref, fnu_ref = np.loadtxt('data/BB_T10000_L100000.dat', usecols=[0,1], unpack=True)

        #fnu_ref /= (4 * np.pi * (10*kpc)**2)
        #fnu_ref *= 1e23

        wav_hyp, fnu_hyp = np.loadtxt(hyperion, usecols=[0,1], unpack=True)

        ax1.loglog(wav_hyp, fnu_hyp, 'o', color='0.5', zorder=1, mec='none')
        ax1.loglog(wav_ref, fnu_ref, 'k-', zorder=2)

        fnu_ref_hyp = np.interp(wav_hyp, wav_ref, fnu_ref)

        ax2.plot(wav_hyp, 100 * (fnu_hyp - fnu_ref_hyp) / fnu_ref_hyp)
        
    ax1.text(0.05, 0.95, "tau = {0}".format(tau_v), ha='left', va='top', transform=ax1.transAxes)
    ax2.text(0.05, 0.95, "tau = {0}".format(tau_v), ha='left', va='top', transform=ax2.transAxes)

    if i == 2:
        ax1.set_xlabel("Wavelength (microns)")
        ax1.set_ylabel("Flux (Jy)")
        ax2.set_xlabel("Wavelength (microns)")
        ax2.set_ylabel("Flux Difference (%)")

    ax1.set_ylim(1.e-6, 200.)
    ax2.set_ylim(-5, 5)
    ax2.set_xscale('log')

fig1.savefig('comparison.png')
fig2.savefig('ratio.png')
