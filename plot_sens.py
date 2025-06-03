#! /usr/bin/env python
#
import matplotlib.pyplot as plt
import numpy as np

DATA_DIR = 'data'
USE_LOG = True
fig, ax = plt.subplots(2, 1, figsize=(8, 6), constrained_layout=True)
#ax =[ax]

ax[0].set_title('Effective Area / System Temperature')
hf = np.loadtxt(f'{DATA_DIR}/AeTsys_HF.dat', delimiter=',')
ax[0].fill_between([8, 3000], [.25, 0.25], color='lightgray')
ax[0].fill_between(hf[:, 0], hf[:, 1], hf[:, 2], color='mediumslateblue')

vl = np.loadtxt(f'{DATA_DIR}/AeTsys_VL.dat', delimiter=',')
ax[0].fill_between(vl[:, 0], vl[:, 1], vl[:, 2], color='mediumslateblue')

vh = np.loadtxt(f'{DATA_DIR}/AeTsys_VH.dat', delimiter=',')
ax[0].fill_between(vh[:, 0], vh[:, 1], vh[:, 2], color='mediumslateblue')

ul = np.loadtxt(f'{DATA_DIR}/AeTsys_UL.dat', delimiter=',')
ax[0].fill_between(ul[:, 0], ul[:, 1], ul[:, 2], color='mediumslateblue')

uh = np.loadtxt(f'{DATA_DIR}/AeTsys_UH.dat', delimiter=',')
ax[0].fill_between(uh[:, 0], uh[:, 1], uh[:, 2], color='mediumslateblue')

if USE_LOG:
    ax[0].set_xscale('log')
    ax[0].set_yscale('log')
    ax[0].axis(xmin=8, xmax=3000, ymin=6E-6, ymax=0.3)
else:
    ax[0].axis(ymin=0.0, ymax=.18)
ax[0].set_xlabel('Freq [MHz]')
ax[0].set_ylabel(r'Ae/Tsys [$m^2 / K$]')
ax[0].grid()

ax[1].set_title('Source Equivalent Flux Density')
shf = np.loadtxt(f'{DATA_DIR}/sefd_hf.dat', delimiter=',')
ax[1].fill_between([8, 3000], [4e8, 4e8], color='lightgray')
ax[1].fill_between(shf[:, 0], shf[:, 1], shf[:, 2], color='mediumslateblue')

svl = np.loadtxt(f'{DATA_DIR}/sefd_vl.dat', delimiter=',')
ax[1].fill_between(svl[:, 0], svl[:, 1], svl[:, 2], color='mediumslateblue')

svh = np.loadtxt(f'{DATA_DIR}/sefd_vh.dat', delimiter=',')
ax[1].fill_between(svh[:, 0], svh[:, 1], svh[:, 2], color='mediumslateblue')

sul = np.loadtxt(f'{DATA_DIR}/sefd_ul.dat', delimiter=',')
ax[1].fill_between(sul[:, 0], sul[:, 1], sul[:, 2], color='mediumslateblue')

suh = np.loadtxt(f'{DATA_DIR}/sefd_uh.dat', delimiter=',')
ax[1].fill_between(suh[:, 0], suh[:, 1], suh[:, 2], color='mediumslateblue')

if USE_LOG:
    ax[1].set_xscale('log')
    ax[1].set_yscale('log')
    ax[1].axis(xmin=8, xmax=3000, ymin=1E4, ymax=4E8)
else:
    ax[1].axis(ymin=0.0, ymax=1.25e6)
ax[1].set_xlabel('Freq [MHz]')
ax[1].set_ylabel('SEFD [Jy]')
ax[1].grid()

plt.show()
plt.savefig('freqbands.png', dpi=300)