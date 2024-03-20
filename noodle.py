import sigs

BW = 50.0e6  # full bandwidth
B = 0.1e6  # channel bandwidth
tau = 2.0  # s
oversample=10
total_time = 2  # s
Tsys = 50.0  # K
Tsky = 100.0
SNR = 0.000001
freq = 23e6
phase = 0.0


A = sigs.SignalProc(Tsys, BW, B, tau, total_time, oversample=10.0, mu=0)
s = A.signal(SNR, freq, phase)
n_sky = A.band_limited_white_noise(Tsky)

n_sys1 = A.band_limited_white_noise(Tsys)
n_sys2 = A.band_limited_white_noise(Tsys)

x1 = s + n_sky + n_sys1
x2 = s + n_sky + n_sys2

u1 = sigs.float_to_uint8(x1)
u2 = sigs.float_to_uint8(x2)

sigs.write_numpy_array_to_binary(u1, 'obs1.bin')
sigs.write_numpy_array_to_binary(u2, 'obs2.bin')