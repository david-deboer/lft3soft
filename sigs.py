import numpy as np
import matplotlib.pyplot as plt


kB = 1.380649e-23
BW = 50.0  # MHz, full bandwidth
B = 1.0  # MHz, channel bandwidth
tau = 2.0  # s
sr = BW * 1E6 * 2.1  # /s
total_time = 10  # s
N = int(total_time * sr)
Tsys = 50.0  # K, see the full analysis...



# AWGN https://stackoverflow.com/questions/14058340/adding-noise-to-a-signal-in-python
def white_noise(Tsys, B, tau, sr, n, mu=0):
    """
    Parameters
    ----------
    Tsys : float
        System temperature in K
    B : float
        Channel bandwidth in MHz
    tau : float
        Integration time in sec
    sr : float
        Sample rate in MHz
    n : int
        Number of samples
    mu : float
        Mean

    """
    # f = np.fft.fftfreq(N, 1 / sr)
    B *= 1e6
    sr *= 1e6
    rho = kB * Tsys / np.sqrt(B * tau)
    sigma = rho * np.sqrt(sr/2.0)
    noise = np.random.normal(mu, sigma, n)
    lpnoise = butter_lowpass_filter(noise, 10.0e6, sr)
    plt.figure('time')
    plt.plot(noise)
    plt.plot(lpnoise)
    plt.figure('Freq')
    plt.plot(np.fft.fft(noise).real)
    plt.plot(np.fft.fft(lpnoise).real)
    #return noise

def butter_lowpass_filter(data, cutoff, fs, order=8):
    from scipy.signal import butter, filtfilt, freqz
    # print("Cutoff freq " + str(cutoff))
    nyq = 0.5 * fs # Nyquist Frequency
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    w, ff = freqz(b, a, fs=fs)
    #plt.semilogy(w, np.abs(ff), label='Filter')
    y = filtfilt(b, a, data)
    return y

def band_limited_noise(min_freq, max_freq, N):
    min_freq = min_freq * 1E6  # Convert MHz to Hz
    max_freq = max_freq * 1E6  # Convert MHz to Hz
    samplerate = 2.1 * max_freq
    freqs = np.abs(np.fft.fftfreq(N, 1/samplerate))
    f = np.zeros(N)
    idx = np.where(np.logical_and(freqs>=min_freq, freqs<=max_freq))[0]
    f[idx] = 1
    plt.plot(freqs, f)
    return fftnoise(f)

def fftnoise(f):
    f = np.array(f, dtype='complex')
    Np = (len(f) - 1) // 2
    phases = np.random.rand(Np) * 2 * np.pi
    phases = np.cos(phases) + 1j * np.sin(phases)
    f[1:Np+1] *= phases
    f[-1:-1-Np:-1] = np.conj(f[1:Np+1])
    return np.fft.ifft(f).real