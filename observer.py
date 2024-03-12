import numpy as np
import matplotlib.pyplot as plt
from copy import copy
import lunar_sys
import lunar_obs


NSIDE = 512


color_palette = [
    (0.12156862745098039, 0.4666666666666667, 0.7058823529411765, 1.0),
    (1.0, 0.4980392156862745, 0.054901960784313725, 1.0),
    (0.17254901960784313, 0.6274509803921569, 0.17254901960784313, 1.0),
    (0.8392156862745098, 0.15294117647058825, 0.1568627450980392, 1.0),
    (0.5803921568627451, 0.403921568627451, 0.7411764705882353, 1.0),
    (0.5490196078431373, 0.33725490196078434, 0.29411764705882354, 1.0)
        ]


class Observe:
    def __init__(self, N=50, deck_diameter=3.0, element_low=400.0, array_low=None, start=300.0, stop=900, step=10):
        self.LB = [1, 50]
        self.FM = [80, 110]
        self.MA = [start, stop]
        self.HB = [start*2, stop*2]
        self.system = lunar_sys.System(N=N, deck_diameter=deck_diameter, element_low=element_low, array_low=array_low, start=start, stop=stop, step=step)
        self.system.gen_Trcvr()
        self.system.gen_Aeff('vivaldi')
        self.system.gen_FWHM()
        self.system.check()
        self.system.show_sys()
        self.Tsys = None
        self.galaxy = None

    def get_galaxy(self, pointing='moon_ptg'):
        self.galaxy = lunar_obs.Galaxy(freqs=self.system.freqs, fwhm=self.system.fwhm, idisplay=self.system.idisplay)
        self.galaxy.pointing_type = pointing
        self.galaxy.gen_map_cube()
        if pointing == 'moon_ptg':
            self.galaxy.gen_pointings()
        else:
            self.galaxy.gen_locs(112000)  # randomish
        self.galaxy.view()

    def get_sky_Tsys(self, pointing='moon_ptg'):
        if self.galaxy is None:
            self.get_galaxy(pointing=pointing)
        self.Gal = []
        self.Tsys = []
        for i in range(len(self.galaxy.locations)):
            self.Gal.append(self.galaxy.map_cube[:, self.galaxy.locations[i]])
            self.Tsys.append(self.galaxy.map_cube[:, self.galaxy.locations[i]] + self.system.Trcvr)

    def get_minmax(self):
        vmin, vmax = 1E6, 0
        for i, T in enumerate(self.Tsys):
            v = self.system.Aeff[0] / T[0]
            if v > vmax:
                self.imax = copy(i)
                vmax = copy(v)
            if v < vmin:
                self.imin = copy(i)
                vmin = copy(v)

    def run(self):
        self.get_sky_Tsys()
        # for i in range(len(self.galaxy.locations)):
            # lon, lat = hp.pixelfunc.pix2ang(NSIDE, i, lonlat=True)
            # hp.visufunc.projplot(lon, lat, 'c*', lonlat=True)
        # Find median galaxy and use

        hp.visufunc.projplot(self.galaxy.coord.lowview.galactic.l.to_value(), self.galaxy.coord.lowview.galactic.b.to_value(), 'k', lonlat=True)
        hp.visufunc.projplot(self.galaxy.coord.hiview.galactic.l.to_value(), self.galaxy.coord.hiview.galactic.b.to_value(), 'k', lonlat=True)

        plt.figure('System')
        # plt.semilogy(self.system.freqs, self.system.Trcvr, linewidth=3, label='Trcvr')
        # plt.semilogy(self.system.freqs, self.system.Aeff, linewidth=3, label='Aeff')
        for G in self.Gal:
            plt.semilogy(self.system.freqs, G)
        plt.grid()
        plt.legend()
        plt.xlabel('Freq [MHz]')
        #plt.ylabel(r'K/m$^2$')
        plt.ylabel('K')

        plt.figure("Sensitivity")
        plt.title(r'$A_e/T_{sys}$')
        for T in self.Tsys:
            plt.plot(self.system.freqs, self.system.Aeff / T)
        plt.grid()
        plt.xlabel('Freq [MHz]')
        plt.ylabel(f'$m^2/K$')
        #peak, bw = peak_bandwidth(system.Aeff / T)

    def plot_bands(self):
        if self.Tsys is None:
            self.get_sky_Tsys()
        self.get_minmax()
        plt.figure('Bands')
        for i, band in enumerate([self.LB, self.FM, self.HB]):
            hgt = 0.04
            plt.barh(0.0, band[1]-band[0], left=band[0], height=hgt, align='center', color=color_palette[i], alpha=0.75)
        ylo = 4.0 * self.system.Aeff / self.Tsys[self.imin]
        yhi = 4.0 * self.system.Aeff / self.Tsys[self.imax]
        plt.fill_between(self.system.freqs, ylo, yhi, color=color_palette[-1])
        locsy, labelsy = plt.yticks([0, 0.35], ['Antenna', 'Array'])
        #plt.setp(labelsy, fontsize=14)
        #plt.axis([0, 2100, 0.5, 3.0])
        plt.xlabel('MHz')

    def peak_bandwidth(self, spectra):
        """
        Returns array of the indices of the peak value and a pair of indices for the 3dB bandwidth
        """
        peak, bw = [], []
        for this_spectrum in spectra:
            smax = max(this_spectrum)
            imax = np.where(this_spectrum > 0.999 * smax)[0]