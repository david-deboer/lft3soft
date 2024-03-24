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
    'oversample': 1.0,
    'sample_time': 0.4,  # s
    'Tsys':  30.0,  # K
    'Tsky':  0.0,
    'freq': 1.9e5,  # Of signal
    'phase':  0.0,
    'drift_rate': 0.0,
    'distance': 10000.0,  # in LY
    'EIRP': 1.0E12,  # in W
    'Ae': 50 * 50.0,
    'SNR': 10.0,
    'integration_time': 18.0
    }

class System:
    def __init__(self, **kwargs):
        self.ant = {}
        self.noise = {}
        self.rx = {}
        self.per_rcvr = argparse.Namespace(chan_noise_v={}, chan_noise_f={},
                                           snr_detected={})
        for p, v in kwargs.items():
            setattr(self, p, v)

    def set_noise(self):
        self.sky = sigs.BandLimitedWhiteNoise(self.Tsky, self.BW, self.sample_time, self.oversample)
        self.sky.band_limited_white_noise()
        for i in range(self.N_rcvr):
            self.ant[i] = sigs.BandLimitedWhiteNoise(self.Tsys, self.BW, self.sample_time, self.oversample)
            self.ant[i].band_limited_white_noise()
            self.noise[i] = sigs.CombinedSignal(self.ant[0], self.sky)
            self.noise[i].power_spectrum()
            self.noise[i].power_from_v()
            self.noise[i].power_from_f()
            self.per_rcvr.chan_noise_v[i] = self.noise[i].Iv2 * self.noise[i].freq_resolution
            self.per_rcvr.chan_noise_f[i] = self.noise[i].If * self.noise[i].freq_resolution

    def set_cw(self):
        self.sig = sigs.DriftingCW(self.ant[0].fs, self.sample_time, starting_freq=self.freq, drift_rate=self.drift_rate, EIRP=self.EIRP, distance=self.distance)
        self.sig.cw(self.Ae, self.phase)
        self.sig.power_spectrum()
        self.sig.power_from_v()
        self.sig.power_from_f()
        self.chan_sig_v = self.sig.Iv2 * self.sig.N
        self.chan_sig_f = self.sig.If * self.sig.N

    def auto_observe(self):
        self.noise[0].integrate(self.integration_time)
        self.freq_resolution = self.noise[0].freq_resolution  # just copy one over
        self.sampled_BW = self.noise[0].fs / 2.0
        for i in range(self.N_rcvr):
            self.rx[i] = sigs.CombinedSignal(self.sig, self.noise[0])
            self.rx[i].power_spectrum()
            self.rx[i].power_from_v()
            self.rx[i].power_from_f()
            self.per_rcvr.snr_detected[i] = self.sig.channel_power/ self.noise[i].channel_power
            self.per_rcvr.chan_noise_v[i] = self.per_rcvr.chan_noise_v[i] #update for integration

    def info(self):
        print(self.noise[0])
        print(f"Integration time = {self.integration_time}")
    
    def auto_info(self):
        print(f"Noise = {self.noise[0].dB('channel_power')}  dB[W]")
        print(f"Signal = {self.sig.dB('channel_power')} dB[W]")
        print(f"SNR(threshhold) = {self.SNR}")
        print(f"SNR(detected) = {self.per_rcvr.snr_detected[0]}")
        print(f"Freq resolution = {self.freq_resolution} Hz")
        print(f"Band = {self.sampled_BW} Hz")
        print("Integrating voltage^2 ...")
        print(f"  self.noise[0] = {sigs.to_dB(self.chan_sig_v)}")
        print(f"  sig = {sigs.to_dB(self.chan_sig_v)}")
        #print(f"  rx[0] = {self.rx[0].dB('Iv2')}")
        print("Integrating power spectrum ...")
        print(f"  self.noise[0] = {self.per_rcvr.chan_noise_f[0]}")
        print(f"  sig = {self.chan_sig_f}")
        #print(f"  rx[0] = {self.rx[0].dB('If')}")
        print("DIFF")
        print(self.noise[0].dB('If') - self.noise[0].dB('Iv2'))
        print(self.sig.dB('If') - self.sig.dB('Iv2'))

    def time_plot(self, plot_span=500):
        t_plot = self.ant[0].t[:plot_span]
        figt, axt = plt.subplots()
        axt.plot(t_plot, self.noise[0].signal[:plot_span], 'b')
        axt.plot(t_plot, self.sig.signal[:plot_span], 'g')
        axt.plot(t_plot, self.rx[0].signal[:plot_span], 'r')
        axt.plot([t_plot[0], t_plot[-1]], [np.sqrt(self.noise[0].channel_power), np.sqrt(self.noise[0].channel_power)], 'k--')
        axt.plot([t_plot[0], t_plot[-1]], [np.sqrt(self.sig.W / 2.0), np.sqrt(self.sig.W / 2.0)], 'g--')

    def freq_plot(self):
        figf, axf = plt.subplots()
        axf.plot(self.rx[0].f, self.rx[0].dB('S'), 'b')
        axf.plot([self.rx[0].f[0], self.rx[0].f[-1]], [self.noise[0].dB('channel_power') + sigs.to_dB(self.SNR), self.noise[0].dB('channel_power') + sigs.to_dB(self.SNR)])
        axf.plot([self.rx[0].f[0], self.rx[0].f[-1]], [self.sig.dB('channel_power'), self.sig.dB('channel_power')])
        axf.set_xlim(left=0, right=self.BW)
        axf.set_ylim(bottom=-220)

    def cross_observe(self):
        print("YEP")

obs = System(**INPUTS)
obs.set_noise()
obs.set_cw()
obs.auto_observe()
obs.info()
obs.auto_info()
obs.time_plot()
obs.freq_plot()
#a = correlate(x, x, mode='same')
#print(type(a))
#print(a.shape)
#axf.plot(a)

# u1 = sigs.float_to_uint8(x1)
# u2 = sigs.float_to_uint8(x2)

# sigs.write_numpy_array_to_binary(u1, 'obs1.bin')
# sigs.write_numpy_array_to_binary(u2, 'obs2.bin')