import numpy as np
from pygdsm import GlobalSkyModel
import matplotlib.pyplot as plt
import healpy as hp
from copy import copy
import moon_ptg


directivity = {'vivaldi': 3.1, 'dipole': 1.6}
NSIDE = 512

# see Danny Price https://github.com/telegraphic/pygdsm/tree/master

color_palette = [
    (0.12156862745098039, 0.4666666666666667, 0.7058823529411765, 1.0),
    (1.0, 0.4980392156862745, 0.054901960784313725, 1.0),
    (0.17254901960784313, 0.6274509803921569, 0.17254901960784313, 1.0),
    (0.8392156862745098, 0.15294117647058825, 0.1568627450980392, 1.0),
    (0.5803921568627451, 0.403921568627451, 0.7411764705882353, 1.0),
    (0.5490196078431373, 0.33725490196078434, 0.29411764705882354, 1.0)
        ]

class System:
    def __init__(self, N, deck_diameter, element_low, array_low, start, stop, step):
        """
        Parameters
        ----------
        N : int
            number of elements
        deck_diameter : float
            diameter of deck in [m]
        element_low : float
            low frequency element "cut-off" MHz
        array_low : float
            critically sampled at MHz
        start : float
            frequency to start MHz
        stop : float
            frequency to stop MHz
        step : float
            step in MHz

        """
        self.N = N
        self.deck_diameter = deck_diameter
        self.element_low = element_low  # Low frequency size scale (MHz)
        if array_low is None:
            self.array_low = element_low
        else:
            self.array_low = array_low  # Critical spacing frequency (MHz)
        if step > 0.0:
            self.idisplay = int((self.element_low - start) / step)
            self.freqs = np.arange(start, stop, step)
        else:
            self.freqs = np.logspace(np.log10(start), np.log10(stop), int(abs(step)))
            adf = abs(self.freqs - self.element_low)
            self.idisplay = np.where(adf < 1.001 * min(adf))[0][0]
        self.fwhm = None

        # Derived:        
        self.deck_area = np.pi * np.power(self.deck_diameter/2.0, 2.0)
        self.f_crit = 150 * np.sqrt(N) / self.deck_diameter  # Frequency where deck/freq go critical
        self.element_width = 0.5 * 300.0 / self.element_low
        self.element_spacing = 300.0 / self.array_low
        self.array_width = np.sqrt(self.N) * self.element_spacing
        self.min_width = np.sqrt(self.N) * self.element_width

    def gen_Trcvr(self):
        m = (45 - 30) / (2500 - 100)
        self.Trcvr = 30 + m * (self.freqs - 100.0)
        plt.figure('Trcvr')
        plt.plot(self.freqs, self.Trcvr)
        plt.xlabel('MHz')
        plt.ylabel(r'T$_R$ [K]')

    def gen_Aeff(self, antenna_type='vivaldi'):
        self.Aeff = self.N * np.power(300.0 / self.freqs, 2.0) * directivity[antenna_type] / (4.0 * np.pi)
        plt.figure('Aeff')
        # plt.plot(self.freqs, self.Aeff)
        # Now accommodate critical spacing
        maxscale = min(self.array_width, self.deck_diameter)
        lambda_crit = 300.0 / self.f_crit
        wl = 300.0 / self.freqs
        icrit = np.where(self.freqs < self.f_crit)
        skirt = np.sqrt(wl[icrit] * (1.0 - (wl[icrit] - wl[0]) / (lambda_crit - wl[0])))
        # plt.plot(self.freqs[icrit], skirt)
        self.Aeff[icrit] = self.Aeff[icrit[0][-1]+1] + skirt
        plt.plot(self.freqs, self.Aeff)
        plt.xlabel('MHz')
        plt.ylabel(f'A$_e$ [m$^2$]')

    def gen_FWHM(self):
        self.fwhm = (300 / self.freqs) / self.deck_diameter

    def check(self):
        if self.f_crit > self.array_low:
            print("\t ---> however your design won't fit with these parameters.")
        if self.element_width > self.element_spacing:
            print(f"Can't have the antennas {self.element_width} be bigger than the spacing {self.element_spacing}")
        if self.min_width > self.deck_diameter:
            print(f"Can't have the array {self.array_width} be bigger than the deck {self.deck_diameter}")

    def show_sys(self):
        print(f"N = {self.N}")
        print(f"deck = {self.deck_diameter} m")
        print(f"element low design freq = {self.element_low} MHz")
        print(f"array low (critical) design freq = {self.array_low} MHz")
        print(f"critical frequency based on deck size (lowest can still have critical) = {self.f_crit} MHz")

class Galaxy:
    def __init__(self, system):
        self.gsm = GlobalSkyModel(include_cmb=True)  # This uses GSM08 (?)
        self.system = system

    def gen_map_cube(self):
        map_raw = self.gsm.generate(self.system.freqs)
        if self.system.fwhm is None:
            self.map_cube = map_raw
        else:
            self.map_cube = []
            for i in range(len(self.system.freqs)):
                self.map_cube.append(hp.sphtfunc.smoothing(map_raw[i], self.system.fwhm[i]))
        self.map_cube = np.array(self.map_cube)

    def gen_pointings(self):
        self.coord = moon_ptg.Pointing()
        self.locations = hp.pixelfunc.ang2pix(NSIDE, self.coord.midview.galactic.l.to_value(), self.coord.midview.galactic.b.to_value(), lonlat=True)

    def gen_locs(self, step=112000):
        minval = min(self.map_cube[0])
        imin = np.where(self.map_cube[0] < minval + 0.05)[0][0]
        maxval = max(self.map_cube[0])
        imax = np.where(self.map_cube[0] > maxval - 0.05)[0][0]
        self.locations = [imin, imax] + list(range(0, self.map_cube.shape[1], step))

    def view(self, view_fullres=True, logged=True):
        plt.figure('Pointings')
        if view_fullres:
            self.gsm.view(self.system.idisplay, logged=True)
        else:
            if logged:
                hp.visufunc.mollview(np.log10(self.map_cube[self.system.idisplay]))
            else:
                hp.visufunc.mollview(self.map_cube[self.system.idisplay])


class Observe:
    def __init__(self, N=50, deck_diameter=3.0, element_low=400.0, array_low=None, start=300.0, stop=900, step=10):
        self.LB = [1, 50]
        self.FM = [80, 110]
        self.MA = [start, stop]
        self.HB = [start*2, stop*2]
        self.system = System(N=N, deck_diameter=deck_diameter, element_low=element_low, array_low=array_low, start=start, stop=stop, step=step)
        self.system.gen_Trcvr()
        self.system.gen_Aeff('vivaldi')
        self.system.gen_FWHM()
        self.system.check()
        self.system.show_sys()
        self.Tsys = None
        self.galaxy = None

    def get_galaxy(self, pointing='moon_ptg'):
        self.galaxy = Galaxy(self.system)
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