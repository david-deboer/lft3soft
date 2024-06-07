import sigs
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import correlate
import argparse
import rf_seti

#https://www.dsprelated.com/showarticle/1004.php  calc power from DFT



### INPUTS
INPUTS = {
    'N_site': 2,
    'site': [rf_seti.site1, rf_seti.site2],
    'BW': 3.0e6,  # full bandwidth
    'Tsky':  100.0,
    'oversample_pc': 1.0,  #percent to oversample from Nyquist
    'time_per_single': 0.4,  # s
    'time_integration': 18.0,  # s
    'time_start': 0.0,  # offset for start time
    'signal_start_freq': 0.8e6,  # Of signal
    'signal_drift_rate': 0.0,
    'distance': 100.0,  # in LY
    'EIRP': 1.0E12,  # in W
    'SNR': 10.0
    }


class SiteResponse:
    def __init__(self):
        self.chan_noise_v = None
        self.chan_noise_f = None
        self.chan_sig_v = None
        self.chan_sig_f = None
        self.snr_detected = None

class Observing:
    def __init__(self, sys):
        self.sys = sys
        # Per site
        for i in range(self.sys.N_site):
            self.sys.site[i].ant = {}
            self.noise = {}
            self.rx = {}
            self.rfi = {}
            self.sig = {}
        # Per baseline
        self.cross = {}
        self.response = []
        for i in self.sys.N_site:
            self.response.append(SiteResponse())

    def make_noise(self):
        print("Making sky and system noise.")
        self.sky = sigs.BandLimitedWhiteNoise(self.sys, None, {'T': 'Tsky'})
        self.sky.band_limited_white_noise()
        for i in range(self.sys.N_site):
            self.ant[i] = sigs.BandLimitedWhiteNoise(self.sys, i, {'T': 'Tsys'})
            self.ant[i].band_limited_white_noise()
            self.noise[i] = sigs.CombinedSignal(self.sys, self.ant[i], self.sky)
            self.noise[i].power_spectrum()
            self.noise[i].power_from_v()
            self.noise[i].power_from_f()
            self.response[i].chan_noise_v = self.noise[i].Iv2 * self.noise[i].sys.resolution_BW
            self.response[i].chan_noise_f = self.noise[i].If * self.noise[i].sys.resolution_BW

    def make_cw(self):
        print("Making technosignature.")
        for i in range(self.sys.N_site):
            self.sig[i] = sigs.DriftingCW(self.sys, i)
            self.sig[i].cw()
            self.sig[i].power_spectrum()
            self.sig[i].power_from_v()
            self.sig[i].power_from_f()
            self.response[i].chan_sig_v = self.sig[i].Iv2 * self.sys.N
            self.response[i].chan_sig_f = self.sig[i].If * self.sys.N

    def make_rfi(self):
        print("Making RFI (fm) per site")
        self.rfi[0] = sigs.FM(self.sys)
        self.rfi[0].band_limited_uniform()
        self.rfi[1] = sigs.FM(self.sys)
        self.rfi[1].silent()

    def auto_observe(self):
        print("Making observation and computing autos.")
        # self.noise[0].integrate()
        for i in range(self.sys.N_site):
            self.rx[i] = sigs.CombinedSignal(self.sys, self.sig[i], self.rfi[i], self.noise[i])
            self.rx[i].power_spectrum()
            self.rx[i].power_from_v()
            self.rx[i].power_from_f()
            self.response[i].snr_detected = self.sig[i].channel_power / self.noise[i].channel_power
        self.f = self.rx[0].f  # pick one

    def cross_observe(self, i=0, j=1):
        if self.sys.N_site > 1:
            cind = f"{i}{j}"
            print(f"Computing cross power {i},{j}.")
            self.f_cross, self.cross[cind] = sigs.cross_power(self.rx[i].signal, self.rx[j].signal, self.sys.fs, self.sys.N)
        else:
            print("Not enough receivers!")
    
    def info(self, site=0):
        print(f"Noise = {self.noise[site].dB('channel_power')}  dB[W]")
        print(f"Signal = {self.sig.dB('channel_power')} dB[W]")
        print(f"SNR(threshhold) = {self.sys.SNR}")
        print(f"SNR(detected) = {self.response[site].snr_detected}")
        print(f"Freq resolution = {self.sys.resolution_BW} Hz")
        print(f"Band = {self.sys.BW} Hz")
        print("Integrating voltage^2 ...")
        print(f"  self.noise[{site}] = {sigs.to_dB(self.response[site].chan_noise_v)}")
        print(f"  sig = {sigs.to_dB(self.chan_sig_v)}")
        #print(f"  rx[0] = {self.rx[ant].dB('Iv2')}")
        print("Integrating power spectrum ...")
        print(f"  self.noise[{site}] = {sigs.to_dB(self.response[site].chan_noise_f)}")
        print(f"  sig = {sigs.to_dB(self.chan_sig_f)}")
        #print(f"  rx[0] = {self.rx[ant].dB('If')}")
        print("DIFF")
        print(self.noise[site].dB('If') - self.noise[site].dB('Iv2'))
        print(self.sig.dB('If') - self.sig.dB('Iv2'))

    def time_plot(self, plot_span=500):
        print("Ignoring time plot.")
        return
        t_plot = self.ant[0].sys.t[:plot_span]
        figt, axt = plt.subplots()
        axt.plot(t_plot, self.noise[0].signal[:plot_span], 'b')
        axt.plot(t_plot, self.sig.signal[:plot_span], 'g')
        axt.plot(t_plot, self.rx[0].signal[:plot_span], 'r')
        axt.plot([t_plot[0], t_plot[-1]], [np.sqrt(self.noise[0].channel_power), np.sqrt(self.noise[0].channel_power)], 'k--')
        axt.plot([t_plot[0], t_plot[-1]], [np.sqrt(self.sig.W / 2.0), np.sqrt(self.sig.W / 2.0)], 'g--')

    def freq_plot(self):
        self.detection = sigs.to_dB(self.noise[0].channel_power * self.sys.SNR)
        figf, axf = plt.subplots()
        for i in range(self.sys.N_site):
            axf.plot(self.f, self.rx[i].dB('S'))
        axf.plot([self.f[0], self.f[-1]], [self.detection, self.detection])
        # axf.plot([self.f[0], self.f[-1]], [self.sig.dB('channel_power'), self.sig.dB('channel_power')])
        for bline, cspec in self.cross.items():
            axf.plot(self.f, sigs.to_dB(cspec[:(self.sys.N//2)]), label=bline)
        if self.sys.N_site > 10:  # ignore for now
            plt.figure("Cross")
            plt.plot(self.f_cross, self.cross['01'].real)
            plt.plot(self.f_cross, self.cross['01'].imag)
        axf.set_xlim(left=0, right=self.sys.BW)
        #axf.set_ylim(bottom=-220)


sys = sigs.System(**INPUTS)
# obs = Observing(sys)
# obs.make_noise()
# obs.make_cw()
# obs.make_rfi()
# obs.auto_observe()
# obs.cross_observe()
# obs.info()
# obs.time_plot()
# obs.freq_plot()




#a = correlate(x, x, mode='same')
#print(type(a))
#print(a.shape)
#axf.plot(a)

# u1 = sigs.float_to_uint8(x1)
# u2 = sigs.float_to_uint8(x2)

# sigs.write_numpy_array_to_binary(u1, 'obs1.bin')
# sigs.write_numpy_array_to_binary(u2, 'obs2.bin')