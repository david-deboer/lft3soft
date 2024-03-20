import numpy as np
import matplotlib.pyplot as plt


kB = 1.380649e-23


class SignalProc:

    def __init__(self, Tsys=None, BW=None, B=None, tau=None, total_time=None, oversample=10.0, mu=0):
        """
        Parameters
        ----------
        Tsys : float
            System temperature in K
        B : float
            Channel bandwidth in Hz
        tau : float
            Integration time in sec
        sr : float
            Sample rate in Hz
        N : int
            Number of samples
        mu : float
            Mean

        """
        self.Tsys = Tsys
        self.BW = BW
        self.total_time = total_time
        self.mu = mu
        self.B = B
        self.tau = tau
        if self.BW is not None:
            self.fs = np.floor(self.BW * 2.0 * (1.0 + oversample/100.0))  # /s
            self.N = N = int(self.total_time * self.fs)
            print(f"Tsys={Tsys}\nBW={BW}\nB={B}\ntau={tau}\nsr={self.fs}\nN={N}")

    def signal(self, SNR, freq, phase=0.0):
        self.freq = freq
        self.phase = phase
        self.t = np.linspace(0, self.total_time, self.N)
        return np.sqrt(SNR * kB * self.Tsys * self.B) *  np.sin(2.0 * np.pi * freq * self.t + phase)

    # AWGN https://stackoverflow.com/questions/14058340/adding-noise-to-a-signal-in-python
    def band_limited_white_noise(self, T):
        # f = np.fft.fftfreq(N, 1 / sr)
        rho2 = kB * T
        sigma = np.sqrt(rho2) # * sr / 2.0)
        return self.butter_lowpass_filter(np.random.normal(self.mu, sigma, self.N), order=8)
        # plt.figure('time')
        # plt.plot(noise)
        # plt.plot(lpnoise)
        # plt.figure('Freq')
        # plt.plot(np.fft.fft(noise).real)
        # plt.plot(np.fft.fft(lpnoise).real)

    def butter_lowpass_filter(self, data, order=8):
        from scipy.signal import butter, filtfilt, freqz
        # print("Cutoff freq " + str(cutoff))
        nyq = 0.5 * self.fs # Nyquist Frequency
        normal_cutoff = self.BW / nyq
        # Get the filter coefficients 
        b, a = butter(order, normal_cutoff, btype='low', analog=False)
        w, ff = freqz(b, a, fs=self.fs)
        # plt.semilogy(w, np.abs(ff), label='Filter')
        y = filtfilt(b, a, data)
        return y

    def fft(self, data):
        return np.fft.fft(data).real

    def Smin(self):
        """
        This just implements Eq 1/2 in https://iopscience.iop.org/article/10.3847/1538-3881/acfc1e/pdf
        for canonical case.

        """
        SNR = 5.0  # Factor of 5 in SNR
        dnu = 1.0  # 1 Hz
        tobs = 5.0 * 60.0  # Five minutes
        npol = 2.0
        Tsys = 25.0
        Ae = 0.85 * np.pi * 50.0 * 50.0  # GBT
        Smin = SNR * 2.0 * kB / (Ae / Tsys) * np.sqrt(dnu / (npol*tobs))
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