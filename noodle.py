import sigs
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import correlate

#https://www.dsprelated.com/showarticle/1004.php  calc power from DFT

### INPUTS
BW = 3.0e6  # full bandwidth
oversample=1.0
total_time = 0.4  # s
Tsys = 30.0  # K
Tsky = 0.0
freq = 1.9e5  # Of signal
phase = 0.0
drift_rate = 0.0
distance =10000.  # in LY
EIRP = 1.0E12  # in W
Ae = 50 * 50.
SNR = 10.0
integration_time = 18.0

### NOISE
sky = sigs.BandLimitedWhiteNoise(Tsky, BW, total_time, oversample)
sky.band_limited_white_noise(integration_time=integration_time)
sys1 = sigs.BandLimitedWhiteNoise(Tsys, BW, total_time, oversample)
sys1.band_limited_white_noise(integration_time=integration_time)
noise1 = sigs.CombinedSignal(sys1, sky)
noise1.power_spectrum()
noise1.integrate_v()
noise1.integrate_f()

## ET
sig = sigs.DriftingCW(sys1.fs, total_time, starting_freq=freq, drift_rate=drift_rate, EIRP=EIRP, distance=distance)
sig.cw(Ae, phase)
sig.power_spectrum()
sig.integrate_v()
sig.integrate_f()

## RECEIVED
rx1 = sigs.CombinedSignal(sig, noise1)
rx1.power_spectrum()
rx1.integrate_v()
rx1.integrate_f()

print(noise1)
print(f"Integration time = {integration_time}")
print(f"Noise = {noise1.dB('channel_power')}  dB[W]")
print(f"Signal = {sig.dB('channel_power')} dB[W]")
print(f"SNR(threshhold) = {SNR}")
print(f"SNR(detected) = {sig.channel_power/ noise1.channel_power}")
print(f"Freq resolution = {sys1.fs / sys1.N} Hz")
print(f"Band = {sys1.fs / 2.0} Hz")
print("Integrating voltage^2 ...")
print(f"  noise1 = {noise1.dB('Iv2') + sigs.to_dB(noise1.freq_resolution)}")
print(f"  sig = {sig.dB('Iv2') + sigs.to_dB(sig.N)}")
print(f"  rx1 = {rx1.dB('Iv2')}")
print("Integrating power spectrum ...")
print(f"  noise1 = {noise1.dB('If') + sigs.to_dB(noise1.freq_resolution)}")
print(f"  sig = {sig.dB('If') + sigs.to_dB(sig.N)}")
print(f"  rx1 = {rx1.dB('If')}")
print("DIFF")
print(noise1.dB('If') - noise1.dB('Iv2'))
print(sig.dB('If') - sig.dB('Iv2'))

### TIME PLOT
plot_span = 500
t_plot = sys1.t[:plot_span]
figt, axt = plt.subplots()
axt.plot(t_plot, noise1.signal[:plot_span], 'b')
axt.plot(t_plot, sig.signal[:plot_span], 'g')
axt.plot(t_plot, rx1.signal[:plot_span], 'r')
axt.plot([t_plot[0], t_plot[-1]], [np.sqrt(noise1.channel_power), np.sqrt(noise1.channel_power)], 'k--')
axt.plot([t_plot[0], t_plot[-1]], [np.sqrt(sig.W / 2.0), np.sqrt(sig.W / 2.0)], 'g--')

### FREQUENCY PLOT
figf, axf = plt.subplots()
axf.plot(rx1.f, rx1.dB('S'), 'b')
axf.plot([rx1.f[0], rx1.f[-1]], [noise1.dB('channel_power') + sigs.to_dB(SNR), noise1.dB('channel_power') + sigs.to_dB(SNR)])
axf.plot([rx1.f[0], rx1.f[-1]], [sig.dB('channel_power'), sig.dB('channel_power')])
axf.set_xlim(left=0, right=BW)
axf.set_ylim(bottom=-220)


#a = correlate(x, x, mode='same')
#print(type(a))
#print(a.shape)
#axf.plot(a)

# u1 = sigs.float_to_uint8(x1)
# u2 = sigs.float_to_uint8(x2)

# sigs.write_numpy_array_to_binary(u1, 'obs1.bin')
# sigs.write_numpy_array_to_binary(u2, 'obs2.bin')