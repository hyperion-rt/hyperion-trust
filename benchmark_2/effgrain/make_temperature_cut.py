import numpy as np
from hyperion.model import ModelOutput
from hyperion.util.constants import au, lsun

RES = 256

mo = ModelOutput('bm2_eff_vor_temperature.rtout')

g = mo.get_quantities()

from scipy.spatial import cKDTree

sites = np.array([g.x, g.y, g.z]).transpose()

tree = cKDTree(sites)

ymin, ymax = 0 * au, 60 * au
zmin, zmax = 0 * au, 60 * au

y = np.linspace(ymin, ymax, RES)
z = np.linspace(zmin, zmax, RES)

Y, Z = np.meshgrid(y, z)
YR = Y.ravel()
ZR = Z.ravel()

for x_cut in [10 * au, 26.666667 * au]:

    XR = np.ones(YR.shape) * x_cut

    map_sites = np.array([XR, YR, ZR]).transpose()

    d, i = tree.query(map_sites, 1)

    print(i.dtype)

    temperature = g['temperature'].array[0][i]

    temperature = temperature.reshape(Y.shape)

    from matplotlib import pyplot as plt

    plt.rc('axes', titlesize=12)
    plt.rc('axes', labelsize=12)
    plt.rc('xtick', labelsize=9)
    plt.rc('ytick', labelsize=9)
    plt.rc('xtick.major', size=2)
    plt.rc('ytick.major', size=2)
    plt.rc('xtick.minor', size=1)
    plt.rc('ytick.minor', size=1)
    plt.rc('font', family='serif')
    plt.rc('axes', linewidth=0.5)
    plt.rc('patch', linewidth=0.5)

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    im = ax.imshow(np.log10(temperature), extent=[ymin / au, ymax / au, zmin / au, zmax / au], interpolation='nearest')
    fig.colorbar(im, label='Log[Temperature]')
    ax.set_xlabel('x (au)')
    ax.set_ylabel('z (au)')
    ax.patch.set_facecolor(plt.cm.jet(0.001))
    fig.savefig('temperature_{0:06.3f}au.png'.format(x_cut / au))