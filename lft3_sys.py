import matplotlib.pyplot as plt
import numpy as np


class Base:
    def plot_beam(self):
        plt.figure(self.name)
        for f in np.linspace(self.start, self.stop, 50):
            self.beam_1D(f=f, plot=True)

    def get_fwhm(self):
        self.fwhm = []
        for f in self.freqs:
            self.beam_1D(f=f)
            a = np.where(self.beam > 0.5)
            self.fwhm.append((self.theta[a[0][-1]] - self.theta[a[0][0]]) * np.pi / 180.0)
        self.obs_fwhm = []
        for idx in self.nearest_indices:
            self.obs_fwhm.append(self.fwhm[idx])
        self.fwhm = np.array(self.fwhm)

    def proc_obs_freqs(self, freqs):
        """
        For each value in arr1, find the index of the closest value in arr2.

        Parameters:
        arr1 (array-like): Array of values to find nearest neighbors for.
        arr2 (array-like): Array of values to search in.

        Returns:
        list: Indices of the closest values in arr2 for each value in arr1.
        """
        if freqs is None:
            self.obs_freqs = []
            self.nearest_indices = []
            return
        if np.isscalar(freqs):
            obs_freqs = np.array([freqs])
        else:
            obs_freqs = np.asarray(freqs)
        ind = np.where((obs_freqs >= self.start) & (obs_freqs <= self.stop))
        self.obs_freqs = obs_freqs[ind]

        self.nearest_indices = []
        for f in self.obs_freqs:
            # Compute absolute differences between val and all elements in arr2
            differences = np.abs(self.freqs - f)
            # Get index of the minimum difference
            idx = np.argmin(differences)
            self.nearest_indices.append(idx)


class HF(Base):
    def __init__(self, obs_freqs=None):
        self.name = 'HF'
        self.start = 1.0
        self.stop = 50.0
        self.N = 1
        self.directivity = 1.6
        self.freqs = np.arange(self.start, self.stop+0.1, 1.0)
        self.idisplay = len(self.freqs) // 2  # Which frequency index to display
        self.TR()
        self.AE()
        self.proc_obs_freqs(obs_freqs)

    def beam_1D(self, f=10.0, plot=False):
        self.theta = np.arange(0, 180, 1.0) + 1E-6
        th = self.theta * np.pi / 180.0
        beam_small = np.sin(th)
        beam_half = np.cos((np.pi / 2.0) * np.cos(th)) / np.sin(th)
        beam_full = (np.cos(np.pi * np.cos(th)) + 1.0) / (2.0 * np.sin(th))
        if f < 10.0:
            self.beam = beam_small
        elif f < 25.0:
            mc = (1.0 / 15.0) * (f - 10.0)
            self.beam = beam_small * (1.0 - mc) + beam_half * mc
        else:
            mc = (1.0 / 25.0) * (f - 25.0)
            self.beam = beam_half * (1.0 - mc) + beam_full * mc
        self.beam = self.beam ** 2.0
        if plot:
            plt.plot(self.theta, self.beam)
    
    def TR(self):
        self.Trcvr = 300.0 * np.ones(len(self.freqs))

    def AE(self):
        a = 0.48 / 2499
        c = 0.02 - a
        backlobe = a * self.freqs * self.freqs + c
        self.Aeff = (self.N * np.power(300.0 / self.freqs, 2.0) * self.directivity / (4.0 * np.pi)) * backlobe
        # plt.plot(self.freqs, backlobe)
    
    def plot(self, ptype='T', ax=None, xscale='linear', yscale='linear', lw=4):
        if ax is None:
            fig, ax = plt.subplots()
        if ptype == 'T':
            y = self.Trcvr
        elif ptype == 'L':
            y = 1.0 * np.ones(len(self.freqs))
        elif ptype == 'A':
            y = self.Aeff
        elif ptype == 'S':
            y = self.Aeff / self.Trcvr
        ax.plot(self.freqs, y, label='HF', lw=lw)
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)
            

class VHF_LO(Base):
    def __init__(self, obs_freqs=None):
        self.name = 'VHF-LO'
        self.start = 60.0
        self.stop = 110.0
        self.N = 1
        self.directivity = 1.6
        self.freqs = np.arange(self.start, self.stop+0.1, 1.0)
        self.idisplay = len(self.freqs) // 2  # Which frequency index to display
        self.TR()
        self.AE()
        self.proc_obs_freqs(obs_freqs)

    def beam_1D(self, f=100.0, plot=False):
        self.theta = np.arange(0, 180, 1.0) + 1E-6
        th = self.theta * np.pi / 180.0
        self.beam = (np.cos((np.pi/ 2.0) * np.cos(th)) / np.sin(th)) ** 2.0
        if plot:
            plt.plot(self.theta, self.beam)

    def TR(self):
        self.Trcvr = 100.0 * np.ones(len(self.freqs))

    def AE(self):
        self.Aeff = self.N * np.power(300.0 / self.freqs, 2.0) * self.directivity / (4.0 * np.pi)

    def plot(self, ptype='T', ax=None, xscale='linear', yscale='linear', lw=4):
        if ax is None:
            fig, ax = plt.subplots()
        if ptype == 'T':
            y = self.Trcvr
        elif ptype == 'L':
            y = 1.0 * np.ones(len(self.freqs))
        elif ptype == 'A':
            y = self.Aeff
        elif ptype == 'S':
            y = self.Aeff / self.Trcvr
        ax.plot(self.freqs, y, label='VHF-LO (VL)', lw=lw)
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)


class VHF_HI(Base):
    def __init__(self, obs_freqs=None):
        self.name = 'VHF-HI'
        self.start = 150.0
        self.stop = 250.0
        self.N = 1
        self.directivity = 1.6
        self.freqs = np.arange(self.start, self.stop+0.1, 1.0)
        self.idisplay = len(self.freqs) // 2  # Which frequency index to display
        self.TR()
        self.AE()
        self.proc_obs_freqs(obs_freqs)

    def beam_1D(self, f=200.0, plot=False):
        self.theta = np.arange(0, 180, 1.0) + 1E-6
        th = self.theta * np.pi / 180.0
        self.beam = (np.cos((np.pi/ 2.0) * np.cos(th)) / np.sin(th)) ** 2.0
        if plot:
            plt.plot(self.theta, self.beam)
    
    def TR(self):
        self.Trcvr = 100.0 * np.ones(len(self.freqs))

    def AE(self):
        self.Aeff = self.N * np.power(300.0 / self.freqs, 2.0) * self.directivity / (4.0 * np.pi)

    def plot(self, ptype='T', ax=None, xscale='linear', yscale='linear', lw=4):
        if ax is None:
            fig, ax = plt.subplots()
        if ptype == 'T':
            y = self.Trcvr
        elif ptype == 'L':
            y = 1.0 * np.ones(len(self.freqs))
        elif ptype == 'A':
            y = self.Aeff
        elif ptype == 'S':
            y = self.Aeff / self.Trcvr
        ax.plot(self.freqs, y, label='VHF-HI (VH)', lw=lw)
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)


class UHF_LO(Base):
    def __init__(self, obs_freqs=None):
        self.name = 'UHF-LO'
        self.start = 300.0
        self.stop = 900.0
        self.N = 48
        self.D = 3.5
        self.directivity = 3.1
        self.freqs = np.arange(self.start, self.stop+0.1, 2.0)
        self.idisplay = len(self.freqs) // 2  # Which frequency index to display
        self.TR()
        self.AE()
        self.proc_obs_freqs(obs_freqs)

    def beam_1D(self, f=600.0, plot=False):
        self.theta = np.arange(0, 180, 1.0) + 1E-6
        th = (self.theta - 90.0) * np.pi / 180.0
        from scipy.special import j1
        x = (2.0 * np.pi / (300.0 / f)) * (self.D / 2.0) * np.sin(th)
        self.beam = (2.0 * j1(x) / x) ** 2
        if plot:
            plt.plot(self.theta, self.beam)

    def TR(self):
        self.Trcvr = 30.0 * np.ones(len(self.freqs))

    def AE(self):
        self.Aeff = self.N * np.power(300.0 / self.freqs, 2.0) * self.directivity / (4.0 * np.pi)

    def plot(self, ptype='T', ax=None, xscale='linear', yscale='linear', lw=4):
        if ax is None:
            fig, ax = plt.subplots()
        if ptype == 'T':
            y = self.Trcvr
        elif ptype == 'L':
            y = 1.0 * np.ones(len(self.freqs))
        elif ptype == 'A':
            y = self.Aeff
        elif ptype == 'S':
            y = self.Aeff / self.Trcvr
        ax.plot(self.freqs, y, label='UHF-LO (UL)', lw=lw)
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)


class UHF_HI(Base):
    def __init__(self, obs_freqs=None):
        self.name = 'UHF-HI'
        self.start = 900.0
        self.stop = 2700.0
        self.N = 8
        self.directivity = 3.1
        self.freqs = np.arange(self.start, self.stop+0.1, 10.0)
        self.idisplay = len(self.freqs) // 2  # Which frequency index to display
        self.TR()
        self.AE()
        self.proc_obs_freqs(obs_freqs)

    def beam_1D(self, f=1800.0, plot=False):
        self.theta = np.arange(0, 180, 1.0) + 1E-6
        th = self.theta * np.pi / 180.0
        self.beam = (np.cos((np.pi/ 2.0) * np.cos(th)) / np.sin(th)) ** 2.0
        if plot:
            plt.plot(self.theta, self.beam)
    
    def TR(self):
        self.Trcvr = 35.0 * np.ones(len(self.freqs))

    def AE(self):
        self.Aeff = self.N * np.power(300.0 / self.freqs, 2.0) * self.directivity / (4.0 * np.pi)

    def plot(self, ptype='T', ax=None, xscale='linear', yscale='linear', lw=4):
        if ax is None:
            fig, ax = plt.subplots()
        if ptype == 'T':
            y = self.Trcvr
        elif ptype == 'L':
            y = 1.0 * np.ones(len(self.freqs))
        elif ptype == 'A':
            y = self.Aeff
        elif ptype == 'S':
            y = self.Aeff / self.Trcvr
        ax.plot(self.freqs, y, label='UHF-HI (UH)', lw=lw)
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)


class System:
    def __init__(self, obs_freqs=None):
        self.obs_freqs = obs_freqs
        self.obs_fwhm = []
        self.HF = HF(obs_freqs=obs_freqs)
        self.HF.get_fwhm()
        for fw in self.HF.obs_fwhm:
            self.obs_fwhm.append(fw)
        self.VHF_LO = VHF_LO(obs_freqs=obs_freqs)
        self.VHF_LO.get_fwhm()
        for fw in self.VHF_LO.obs_fwhm:
            self.obs_fwhm.append(fw)
        self.VHF_HI = VHF_HI(obs_freqs=obs_freqs)
        self.VHF_HI.get_fwhm()
        for fw in self.VHF_HI.obs_fwhm:
            self.obs_fwhm.append(fw)
        self.UHF_LO = UHF_LO(obs_freqs=obs_freqs)
        self.UHF_LO.get_fwhm()
        for fw in self.UHF_LO.obs_fwhm:
            self.obs_fwhm.append(fw)
        self.UHF_HI = UHF_HI(obs_freqs=obs_freqs)
        self.UHF_HI.get_fwhm()
        for fw in self.UHF_HI.obs_fwhm:
            self.obs_fwhm.append(fw)
        self.obs_fwhm = np.array(self.obs_fwhm)
        
    def plot_freq_plan(self, ptype='T', xscale='linear', yscale='linear'):
        fig, ax = plt.subplots()
        self.HF.plot(ptype=ptype, ax=ax, xscale=xscale, yscale=yscale)
        self.VHF_LO.plot(ptype=ptype, ax=ax, xscale=xscale, yscale=yscale)
        self.VHF_HI.plot(ptype=ptype, ax=ax, xscale=xscale, yscale=yscale)
        self.UHF_LO.plot(ptype=ptype, ax=ax, xscale=xscale, yscale=yscale)
        self.UHF_HI.plot(ptype=ptype, ax=ax, xscale=xscale, yscale=yscale)
        ax.legend()
        ax.set_xlabel('Frequency (MHz)')
        ax.set_ylabel(ptype)

    def plot_beams(self):
        fig, ax = plt.subplots()
        self.HF.plot_beam()
        self.VHF_LO.plot_beam()
        self.VHF_HI.plot_beam()
        self.UHF_LO.plot_beam()
        self.UHF_HI.plot_beam()