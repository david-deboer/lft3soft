import numpy as np
# import matplotlib.pyplot as plt

rng = np.random.default_rng()


kB = 1.380649e-23  # J/K
LY = 9.461e+15  # m/ly


def to_dB(v):
    return 10.0 * np.log10(v)

def butter_lowpass_filter(data, fs, BW, order=8):
    from scipy.signal import butter, filtfilt, freqz
    nyq = 0.5 * fs # Nyquist Frequency
    normal_cutoff = BW / nyq
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    # w, ff = freqz(b, a, fs=fs)
    # plt.semilogy(w, np.abs(ff), label='Filter')
    y = filtfilt(b, a, data)
    return y

class Signal:
    def __repr__(self):
        return f"T={self.T}\nBW={self.BW}\nfs={self.fs}\nN={self.N}"

    def power_spectrum(self):
        f, sv = fft(self.fs, self.signal)
        S = 2.0 * (np.abs(sv))**2 / self.N
        pn = len(f) // 2
        self.f = f[:pn]
        self.S = S[:pn]

    def dB(self, a):
        v = getattr(self, a)
        return to_dB(v)

    def delay(self, t_delay=0.0):
        """
        Parameter
        ---------
        t_delay : float
            start delay in microseconds
        """
        print("Delay the signal.")

    def integrate(self, integration_time):
        self.integration_time = integration_time
        self.bTau = np.sqrt(self.integration_time * self.freq_resolution)
        print("NOW NEED TO APPLY bTau AS APPROPRIATE")

    def power_from_v(self):
        self.Iv2 =  np.trapz(self.signal**2, self.t) / self.total_time
        
    def power_from_f(self):
        self.If =  np.trapz(self.S, self.f) * self.total_time / self.N

class FM(Signal):
    def __init__(self, fs, BW, N, station=1e6, modu=2e3):
        self.fs = fs
        self.BW = BW
        self.N = N
        self.station = station
        self.modu = modu
        
    def band_limited_uniform(self):
        self.signal = butter_lowpass_filter(rng.uniform(-1.0, 1.0, self.N), self.fs, self.modu)


class CombinedSignal(Signal):
    def __init__(self, signal1, signal2):
        self.signal = signal1.signal + signal2.signal
        self.channel_power = signal1.channel_power + signal2.channel_power
        self.fs = signal1.fs
        self.N = signal1.N
        self.t = signal1.t
        self.freq_resolution = signal1.freq_resolution
        self.t = signal1.t
        self.BW = signal1.BW
        self.total_time = signal1.total_time
        self.T = 0.0
        try:
            self.T += signal1.T
        except AttributeError:
            pass
        try:
            self.T += signal2.T
        except AttributeError:
            pass

class BandLimitedWhiteNoise(Signal):
    def __init__(self, T=50.0, BW=50.0e6, total_time=1.0, oversample=2.0):
        """
        Parameters
        ----------
        T : float
            System temperature in K
        BW : float
            Bandwidth in Hz
        total_time : float
            Total time to "record" voltages in s
        oversample : float
            Percentage to oversample (primarily to see the roll-off)

        """
        self.T = T
        self.BW = BW
        self.total_time = total_time
        self.mu = 0.0  # Could be a bias someday...?
        self.fs = np.floor(self.BW * 2.0 * (1.0 + oversample/100.0))  # /s
        self.N = int(self.total_time * self.fs)
        self.freq_resolution = self.fs / self.N
        self.t = np.linspace(0.0, self.total_time, self.N)
        self.channel_power = kB * self.T * self.freq_resolution
        self.total_power = kB * self.T * self.BW

    # AWGN https://stackoverflow.com/questions/14058340/adding-noise-to-a-signal-in-python
    def band_limited_white_noise(self):
        var = kB * self.T
        sigma = np.sqrt(var)
        self.signal = butter_lowpass_filter(rng.normal(self.mu, sigma, self.N), self.fs, self.BW, order=4)


class DriftingCW(Signal):
    def __init__(self, fs, total_time, starting_freq, drift_rate=0.0, EIRP=1.0E12, distance=10.0):
        """
        Parameters
        ----------
        fs : float
            sample rate in Hz
        total_time : float
            total_time of recording in s
        starting_freq : float
            starting frequency in Hz
        drift_rate : float
            linear drift rate in Hz/s
        EIRP : float
            effective isotropic radiated power in W
        distance : float
            distance to signal in ly

        """
        self.fs = fs
        self.total_time = total_time
        self.N = int(self.total_time * self.fs)
        self.freq_resolution = self.fs / self.N
        self.BW = self.fs / 2.0
        self.starting_freq = starting_freq
        self.drift_rate = drift_rate
        self.EIRP = EIRP
        self.distance = distance * LY
        self.Wm2 = self.EIRP / (4.0 * np.pi * self.distance**2)
        self.t = np.linspace(0, self.total_time, self.N)

    def cw(self, Ae=1.0, phase=0.0):
        self.A = Ae
        self.W = self.Wm2 * Ae
        self.phase = phase
        Vm = np.sqrt(self.W)
        self.freq = np.linspace(self.starting_freq, self.starting_freq + self.total_time * self.drift_rate, self.N)
        self.signal = Vm *  np.sin(2.0 * np.pi * self.freq * self.t + phase)
        self.channel_power = self.W * self.N / 2.0


def fft(fs, data):
    f = np.fft.fftfreq(len(data), 1 / fs)
    sigf = np.fft.fft(data)
    return f, sigf

def Smin(self):
    """
    This just implements Eq 1/2 in https://iopscience.iop.org/article/10.3847/1538-3881/acfc1e/pdf
    for canonical case.

    """
    SNR = 5.0  # Factor of 5 in SNR
    dnu = 1.0  # 1 Hz
    tobs = 5.0 * 60.0  # Five minutes
    npol = 2.0
    T = 25.0
    Ae = 0.85 * np.pi * 50.0 * 50.0  # GBT
    Smin = SNR * 2.0 * kB / (Ae / T) * np.sqrt(dnu / (npol*tobs))
    d = 4.2 * 9.461e+15  # m to nearest star
    EIRPmin = 4.0 * np.pi * d * d * Smin
    print(Smin, EIRPmin / 1e9)

import struct
def float_to_uint8(array):
    min_val = array.min()
    max_val = array.max()
    scaled_array = (array - min_val) / (max_val - min_val) * (2**8 - 1)
    uint_array = np.clip(scaled_array, 0, 2**8 - 1).astype(np.uint8)
    return uint_array

def write_numpy_array_to_binary(array, filename, bits_per_element=8):
    """
    Writes a NumPy array to a binary file in a fixed number of bits per element.

    Parameters:
        array (numpy.ndarray): The NumPy array to be written.
        filename (str): The name of the binary file to write to.
        bits_per_element (int): The number of bits per element for encoding.

    Returns:
        None
    """
    dtype = array.dtype
    if dtype not in [np.uint8, np.uint16, np.uint32, np.uint64]:
        raise ValueError("Unsupported data type. Only uint8, uint16, uint32, uint64 are supported.")
    with open(filename, 'wb') as f:
        array.tofile(f)
    
    # # Calculate the number of bytes required per element
    # bytes_per_element = (bits_per_element + 7) // 8
    # print(bytes_per_element)
    
    # with open(filename, 'wb') as f:
    #     f.write(struct.pack('<I', bits_per_element))  # Write bits per element as 4-byte unsigned int
    #     f.write(struct.pack('<I', len(array)))   # Write number of rows as 4-byte unsigned int
    #     # f.write(struct.pack('<I', array.shape[1]))   # Write number of columns as 4-byte unsigned int
        
    #     # for row in array:
    #     for value in array:
    #         # Convert value to binary string with specified number of bits
    #         bin_value = format(value, f'0{bits_per_element}b')
            
    #         # Pack binary string into bytes and write to file
    #         for x in map(int, bin_value):
    #             f.write(struct.pack(f'<{bytes_per_element}B', x))

# Example usage:
# arr = np.array([1, 2, 3], dtype=np.uint8)
# write_numpy_array_to_binary(arr, 'array.bin', bits_per_element=4)

def read_numpy_array_from_binary(filename):
    """
    Reads a NumPy array from a binary file.

    Parameters:
        filename (str): The name of the binary file to read from.

    Returns:
        numpy.ndarray: The NumPy array read from the binary file.
    """
    with open(filename, 'rb') as f:
        # Read bits per element as 4-byte unsigned int
        bits_per_element = struct.unpack('<I', f.read(4))[0]
        
        # Read number of rows as 4-byte unsigned int
        rows = struct.unpack('<I', f.read(4))[0]
        print(bits_per_element, rows)
        # Read number of columns as 4-byte unsigned int
        # cols = struct.unpack('<I', f.read(4))[0]
        cols = 1
        
        # Calculate bytes per element and number of bits to read
        bytes_per_element = (bits_per_element + 7) // 8
        bits_to_read = rows * cols * bits_per_element
        
        # Read binary data from file and convert to NumPy array
        binary_data = f.read(bits_to_read)
        int_values = struct.unpack(f'<{rows * cols}{bytes_per_element}B', binary_data)
        
        # Convert integers to binary strings, then concatenate and reshape to original shape
        binary_strings = [''.join(format(i, f'08b') for i in int_values)]
        array = np.array([int(binary_strings[0][i:i+bits_per_element], 2) for i in range(0, len(binary_strings[0]), bits_per_element)], dtype=np.uint8)
        # array = array.reshape((rows, cols))
        
        return array

# Example usage:
# arr_read = read_numpy_array_from_binary('array.bin')
# print(arr_read)

    # def psd(self, sig):
    #     from scipy import signal
    #     return signal.periodogram(sig, self.fs)  #returns f, psd

    # def band_limited_noise(min_freq, max_freq, N):
    #     min_freq = min_freq * 1E6  # Convert MHz to Hz
    #     max_freq = max_freq * 1E6  # Convert MHz to Hz
    #     samplerate = 2.1 * max_freq
    #     freqs = np.abs(np.fft.fftfreq(N, 1/samplerate))
    #     f = np.zeros(N)
    #     idx = np.where(np.logical_and(freqs>=min_freq, freqs<=max_freq))[0]
    #     f[idx] = 1
    #     plt.plot(freqs, f)
    #     return fftnoise(f)

    # def fftnoise(f):
    #     f = np.array(f, dtype='complex')
    #     Np = (len(f) - 1) // 2
    #     phases = np.random.rand(Np) * 2 * np.pi
    #     phases = np.cos(phases) + 1j * np.sin(phases)
    #     f[1:Np+1] *= phases
    #     f[-1:-1-Np:-1] = np.conj(f[1:Np+1])
    #     return np.fft.ifft(f).real