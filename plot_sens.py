import matplotlib.pyplot as plt
import numpy as np

DATA_DIR = 'data'
fig, ax = plt.subplots(2, 1, figsize=(8, 6), constrained_layout=True)

ax[0].set_title('Effective Area / System Temperature')
lo = np.loadtxt(f'{DATA_DIR}/AeTsys_lo.dat', delimiter=',')
ax[0].fill_between(lo[:, 0], lo[:, 1], lo[:, 2], color='mediumslateblue')

fm = np.loadtxt(f'{DATA_DIR}/AeTsys_fm.dat', delimiter=',')
ax[0].fill_between(fm[:, 0], fm[:, 1], fm[:, 2], color='mediumslateblue')

mid = np.loadtxt(f'{DATA_DIR}/AeTsys_mid.dat', delimiter=',')
ax[0].fill_between(mid[:, 0], mid[:, 1], mid[:, 2], color='mediumslateblue')

hi = np.loadtxt(f'{DATA_DIR}/AeTsys_hi.dat', delimiter=',')
ax[0].fill_between(hi[:, 0], hi[:, 1], hi[:, 2], color='mediumslateblue')

ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[0].set_xlabel('Freq [MHz]')
ax[0].set_ylabel(r'Ae/Tsys [$m^2 / K$]')
ax[0].grid()

ax[1].set_title('Source Equivalent Flux Density')
lo = np.loadtxt(f'{DATA_DIR}/sefd_lo.dat', delimiter=',')
ax[1].fill_between(lo[:, 0], lo[:, 1], lo[:, 2], color='mediumslateblue')

fm = np.loadtxt(f'{DATA_DIR}/sefd_fm.dat', delimiter=',')
ax[1].fill_between(fm[:, 0], fm[:, 1], fm[:, 2], color='mediumslateblue')

mid = np.loadtxt(f'{DATA_DIR}/sefd_mid.dat', delimiter=',')
ax[1].fill_between(mid[:, 0], mid[:, 1], mid[:, 2], color='mediumslateblue')

hi = np.loadtxt(f'{DATA_DIR}/sefd_hi.dat', delimiter=',')
ax[1].fill_between(hi[:, 0], hi[:, 1], hi[:, 2], color='mediumslateblue')

ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[1].set_xlabel('Freq [MHz]')
ax[1].set_ylabel('SEFD [Jy]')
ax[1].grid()