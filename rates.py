from argparse import Namespace
import numpy as np
from tabulate import tabulate

bytes_per_sample = 1
sample_rate = 100E6  # Msps
total_bw = 50E6
pol = 2
bands = {'LB': {'beams': 1, 'tunings': 1, 'pols': 2},
         'FM': {'beams': 1, 'tunings': 1, 'pols': 2},
         'MB': {'beams': 6, 'tunings': 12, 'pols': 2},
         'HB': {'beams': 1, 'tunings': 12, 'pols': 2}}

concurrent = 6 + 1 + 1 + 1
multi50 = 12 + 12 + 1 + 1

full_baseband = 0
for band in bands:
    for beam in range(bands[band]['beams']):
        for tuning in range(bands[band]['tunings']):
            for pol in range(bands[band]['pols']):
                full_baseband += 1
BB = full_baseband * sample_rate * bytes_per_sample
print(f"Baseband {BB / 1E9} GBps")

Iwf50 = Namespace(B = [1, 100, 1000], tau=2,  Ni=100)
print(f"tau: {Iwf50.tau} s")
print(f"Ni: {Iwf50.Ni}")

table_data = []
for B in Iwf50.B:
    table_row = [B]
    Nch = total_bw / B
    table_row.append(Nch / 1E6)
    Nwf = Nch * Iwf50.Ni
    table_row.append(Nwf / 1E6)
    T = Iwf50.Ni * Iwf50.tau
    table_row.append(T)
    Nsamp = Nwf * multi50
    table_row.append(Nsamp / 1E9)
    Ttime = T * multi50
    table_row.append(Ttime)
    rate = Nwf / T
    table_row.append(rate)
    table_data.append(table_row)

print(tabulate(table_data, ['B', 'Nch [M]', 'Nwf [M]', 'T [s]', 'Nsamp [G]', 'Ttime [s]', 'rate [sps]']))
