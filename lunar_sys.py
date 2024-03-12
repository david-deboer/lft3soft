import matplotlib.pyplot as plt
import numpy as np

directivity = {'vivaldi': 3.1, 'dipole': 1.6}


class System:
    def __init__(self, N=50, deck_diameter=3.0, element_low=400.0, array_low=400.0, start=300.0, stop=900.0, step=20.0):
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
