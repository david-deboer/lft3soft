import sigs
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import correlate
import argparse

#https://www.dsprelated.com/showarticle/1004.php  calc power from DFT

### INPUTS
INPUTS = {
    'N_rcvr': 1,
    'BW': 3.0e6,  # full bandwidth
    'oversample_pc': 1.0,  #percent to oversample from Nyquist
    'time_obs_0': 0.4,  # s
    'time_integration': 18.0,  # s
    'Tsys':  30.0,  # K
    'Tsky':  0.0,
    'fm_freq': 1.1e6,
    'fm_mod': 5.e5,
    'fm_power': (10.0 * 1000.0) / (4.0 * np.pi * 1000.0**2),
    'signal_start_freq': 1.9e6,  # Of signal
    'signal_phase':  0.0,
    'signal_drift_rate': 0.0,
    'distance': 1000.0,  # in LY
    'EIRP': 1.0E12,  # in W
    'Ae': 50 * 50.0,
    'SNR': 10.0
    }


class Observing:
    def __init__(self, sys):
        self.ant = {}
        self.noise = {}
        self.rx = {}
        self.per_rcvr = argparse.Namespace(chan_noise_v={}, chan_noise_f={}, snr_detected={})
        self.sys = sys

    def set_noise(self):
        print("Making sky and system noise.")
        self.sky = sigs.BandLimitedWhiteNoise(self.sys, {'T': 'Tsky'})
        self.sky.band_limited_white_noise()
        for i in range(self.sys.N_rcvr):
            self.ant[i] = sigs.BandLimitedWhiteNoise(self.sys, {'T': 'Tsys'})
            self.ant[i].band_limited_white_noise()
            self.noise[i] = sigs.CombinedSignal(self.sys, self.ant[0], self.sky)
            self.noise[i].power_spectrum()
            self.noise[i].power_from_v()
            self.noise[i].power_from_f()
            self.per_rcvr.chan_noise_v[i] = self.noise[i].Iv2 * self.noise[i].sys.resolution_BW
            self.per_rcvr.chan_noise_f[i] = self.noise[i].If * self.noise[i].sys.resolution_BW

    def set_cw(self):
        print("Making technosignature.")
        self.sig = sigs.DriftingCW(self.sys)
        self.sig.cw()
        self.sig.power_spectrum()
        self.sig.power_from_v()
        self.sig.power_from_f()
        self.chan_sig_v = self.sig.Iv2 * self.sys.N
        self.chan_sig_f = self.sig.If * self.sys.N

    def set_rfi(self):
        print("Making RFI (fm)")
        self.fm = sigs.FM(self.sys)
        self.fm.band_limited_uniform()

    def auto_observe(self):
        print("Making observation.")
        self.noise[0].integrate()
        for i in range(self.sys.N_rcvr):
            self.rx[i] = sigs.CombinedSignal(self.sys, self.sig, self.fm, self.noise[0])
            self.rx[i].power_spectrum()
            self.rx[i].power_from_v()
            self.rx[i].power_from_f()
            self.per_rcvr.snr_detected[i] = self.sig.channel_power/ self.noise[i].channel_power
            self.per_rcvr.chan_noise_v[i] = self.per_rcvr.chan_noise_v[i] #update for integration
    
    def info(self):
        print(f"Noise = {self.noise[0].dB('channel_power')}  dB[W]")
        print(f"Signal = {self.sig.dB('channel_power')} dB[W]")
        print(f"SNR(threshhold) = {self.sys.SNR}")
        print(f"SNR(detected) = {self.per_rcvr.snr_detected[0]}")
        print(f"Freq resolution = {self.sys.resolution_BW} Hz")
        print(f"Band = {self.sys.BW} Hz")
        print("Integrating voltage^2 ...")
        print(f"  self.noise[0] = {sigs.to_dB(self.per_rcvr.chan_noise_v[0])}")
        print(f"  sig = {sigs.to_dB(self.chan_sig_v)}")
        #print(f"  rx[0] = {self.rx[0].dB('Iv2')}")
        print("Integrating power spectrum ...")
        print(f"  self.noise[0] = {sigs.to_dB(self.per_rcvr.chan_noise_f[0])}")
        print(f"  sig = {sigs.to_dB(self.chan_sig_f)}")
        #print(f"  rx[0] = {self.rx[0].dB('If')}")
        print("DIFF")
        print(self.noise[0].dB('If') - self.noise[0].dB('Iv2'))
        print(self.sig.dB('If') - self.sig.dB('Iv2'))

    def time_plot(self, plot_span=500):
        t_plot = self.ant[0].sys.t[:plot_span]
        figt, axt = plt.subplots()
        axt.plot(t_plot, self.noise[0].signal[:plot_span], 'b')
        axt.plot(t_plot, self.sig.signal[:plot_span], 'g')
        axt.plot(t_plot, self.rx[0].signal[:plot_span], 'r')
        axt.plot([t_plot[0], t_plot[-1]], [np.sqrt(self.noise[0].channel_power), np.sqrt(self.noise[0].channel_power)], 'k--')
        axt.plot([t_plot[0], t_plot[-1]], [np.sqrt(self.sig.W / 2.0), np.sqrt(self.sig.W / 2.0)], 'g--')

    def freq_plot(self):
        figf, axf = plt.subplots()
        axf.plot(self.rx[0].f, self.rx[0].dB('S'), 'b')
        axf.plot([self.rx[0].f[0], self.rx[0].f[-1]], [self.noise[0].dB('channel_power') + sigs.to_dB(self.sys.SNR), self.noise[0].dB('channel_power') + sigs.to_dB(self.sys.SNR)])
        axf.plot([self.rx[0].f[0], self.rx[0].f[-1]], [self.sig.dB('channel_power'), self.sig.dB('channel_power')])
        axf.set_xlim(left=0, right=self.sys.BW)
        axf.set_ylim(bottom=-220)

    def cross_observe(self):
        print("YEP")

sys = sigs.System(**INPUTS)
fm = sigs.FM(sys)
fm.band_limited_uniform()
fm.power_spectrum()
plt.plot(fm.f, fm.dB('S'))
# plt.plot(fm.signal[:1000])
# obs = Observing(sys)
# obs.set_noise()
# obs.set_cw()
# obs.auto_observe()
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