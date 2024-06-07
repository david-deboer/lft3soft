import numpy as np
# import matplotlib.pyplot as plt
from tabulate import tabulate

rng = np.random.default_rng()


kB = 1.380649e-23  # J/K
LY = 9.461e+15  # m/ly


def to_dB(v):
    return 10.0 * np.log10(np.abs(v))

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

def cross_power(obs1, obs2, fs, N):
    f, sv1 = fft(fs, obs1)
    f, sv2 = fft(fs, obs2)
    S = (sv1 * np.conjugate(sv2)) / N
    pn = len(f) // 2
    return f, S

class System:
    input_constrained = ['fs', 'BW', 'time_per_single', 'N']
    input_to_modify = ['distance']
    def __init__(self, **kwargs):
        M = 2.0 * (1.0 + kwargs['oversample_pc']/100.0)
        if 'fs' in kwargs:
            self.fs = kwargs['fs']
            if 'BW' in kwargs:
                raise ValueError("Can't specify fs and BW")
            self.BW = self.fs / M
        elif 'BW' in kwargs:
            self.BW = kwargs['BW']
            if 'fs' in kwargs:
                raise ValueError("Can't specify BW and fs")
            self.fs = M * self.BW
        else:
            raise ValueError("Need fs or BW")
        if 'time_per_single' in kwargs:
            self.time_per_single = kwargs['time_per_single']
            if 'N' in kwargs:
                raise ValueError("Can't specify N and time_per_single")
            self.N = int(self.time_per_single * self.fs)
        elif 'N' in kwargs:
            self.N = kwargs['N']
            if 'time_per_single' in kwargs:
                raise ValueError("Can't specify time_per_single and N")
            self.time_per_single = self.N / self.fs
        else:
            raise ValueError("Need time_per_single or N")
        self.resolution_BW = self.fs / self.N
        self.t = np.linspace(0.0, self.time_per_single, self.N)
        table_data = [['resolution_BW', self.resolution_BW, 'Hz']]
        for p in self.input_constrained:
            table_data.append([p, getattr(self, p), ''])
        for p, v in kwargs.items():
            if p in self.input_to_modify:  # This is very ugly
                if p == 'distance':
                    v = v * LY
            if p not in self.input_constrained:
                setattr(self, p, v)
                table_data.append([p, v, ''])
        self.parameter_table = tabulate(table_data, headers=['Parameter', 'Value', 'Unit'])
        print(self)

    def __repr__(self):
        return self.parameter_table

class Signal:
    def power_spectrum(self):
        f, sv = fft(self.sys.fs, self.signal)
        S = 2.0 * (np.abs(sv))**2 / self.sys.N
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

    def integrate(self):
        self.bTau = np.sqrt(self.sys.time_integration * self.sys.resolution_BW)
        print("NOW NEED TO APPLY bTau AS APPROPRIATE")

    def power_from_v(self):
        self.Iv2 =  np.trapz(self.signal**2, self.sys.t) / self.sys.time_per_single
        
    def power_from_f(self):
        self.If =  np.trapz(self.S, self.f) * self.sys.time_per_single / self.sys.N

class FM(Signal):
    name = "FM"
    def __init__(self, sys, n):
        self.sys = sys
        self.system_number = n
        if sys is not None:
            self.mod = self.sys.site[n].fm.mod
            self.freq = self.sys.site[n].fm.freq
            self.pwr = self.sys.site[n].fm.Tx
        
    def band_limited_uniform(self):
        if self.sys is None:
            self.silent()
        else:
            mod_f = butter_lowpass_filter(rng.uniform(-1.0, 1.0, self.sys.N), self.sys.fs, self.mod, order=4)
            mod_f /= max(mod_f[100:-100])
            omega = 2.0 * np.pi * (self.freq + mod_f) 
            self.signal = self.pwr * np.cos(omega * self.sys.t)
            self.channel_power = 0.0  # Ignore this for now

    def silent(self):
        self.signal = np.zeros(self.sys.N)
        self.channel_power = 0.0

class CombinedSignal(Signal):
    name = "CombinedSignal"
    def __init__(self, sys, *args):
        self.sys = sys
        self.signal = np.zeros(self.sys.N)
        self.channel_power = 0.0
        self.T = 0.0
        print("Combining:", end=' ')
        for this_signal in args:
            print(this_signal.name, end=' ')
            self.signal += this_signal.signal
            self.channel_power += this_signal.channel_power
            try:
                self.T += this_signal.T
            except AttributeError:
                pass
            try:
                self.T += this_signal.T
            except AttributeError:
                pass
        print()
class BandLimitedWhiteNoise(Signal):
    name = "BandLimitedWhiteNoise"
    def __init__(self, sys, n, parmap):
        """
        Parameters
        ----------
        sys : System instance
            System

        """
        self.sys = sys
        self.system_number = n
        for p_self, p_sys in parmap.items():
            if n is None:
                setattr(self, p_self, getattr(self.sys.site, p_sys))
            else:
                setattr(self, p_self, getattr(self.sys.site, p_sys)[n])
        self.mu = 0.0
        self.channel_power = kB * self.T * self.sys.resolution_BW
        self.total_power = kB * self.T * self.sys.BW

    # AWGN https://stackoverflow.com/questions/14058340/adding-noise-to-a-signal-in-python
    def band_limited_white_noise(self):
        var = kB * self.T
        sigma = np.sqrt(var)
        self.signal = butter_lowpass_filter(rng.normal(self.mu, sigma, self.sys.N), self.sys.fs, self.sys.BW, order=4)


class DriftingCW(Signal):
    name = "DriftingCW"
    def __init__(self, sys):
        """
        Parameters
        ----------
        fs : float
            sample rate in Hz
        time_per_single : float
            time_per_single of recording in s
        starting_freq : float
            starting frequency in Hz
        drift_rate : float
            linear drift rate in Hz/s
        EIRP : float
            effective isotropic radiated power in W
        distance : float
            distance to signal in ly

        """
        self.sys = sys
        self.Wm2 = self.sys.EIRP / (4.0 * np.pi * self.sys.distance**2)

    def cw(self):
        self.W = self.Wm2 * self.sys.Ae
        Vm = np.sqrt(self.W)
        tstop = self.sys.signal_start_freq + self.sys.time_per_single * self.sys.signal_drift_rate
        self.freq_drifted = np.linspace(self.sys.signal_start_freq, tstop, self.sys.N)
        self.signal = Vm *  np.sin(2.0 * np.pi * self.freq_drifted * self.sys.t + self.sys.signal_phase)
        self.channel_power = self.W * self.sys.N / 2.0

##########################################################################################################
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