import numpy as np
import matplotlib.pyplot as plt
from astroquery.jplhorizons import Horizons


obs_lon = 182.0
obs_lat = -23.0

AU = 1.496e+11

Rearth = 6378.
Rmoon = 1737.

start_time = '2024-01-01T0:00:00'
stop_time = '2024-12-31T23:59:59'
step = '1d'

mobj = Horizons(id=301, location=500, epochs={'start': start_time, 'stop': stop_time, 'step': step})
meph = mobj.ephemerides()
w = np.where(meph['PDObsLon'] > 180.0)
meph['PDObsLon'][w] = meph['PDObsLon'][w] - 360.0


plt.figure('sublunar')
plt.plot(meph['PDObsLon'], meph['PDObsLat'])

def draw_circle(r, xo=0.0, yo=0.0):
    th = np.arange(0, 2.0 * np.pi, 0.01)
    xc = r * np.cos(th) + xo
    yc = r * np.sin(th) + yo
    return xc, yc

plt.figure('system')
# Draw the earth
xe, ye = draw_circle(Rearth, 0.0, 0.0)
plt.plot(xe, ye)
# Draw the moon per step
for i, dist in enumerate(meph['delta']):
    D = dist * AU / 1000.0
    dec_moon = meph['DEC'][i] * np.pi / 180.0
    xm, ym = draw_circle(Rmoon, D * np.cos(dec_moon), D * np.sin(dec_moon))
    plt.plot(xm, ym)
plt.axis('square')
plt.grid()

plt.figure('ra_v_dec')
ra = meph['RA'] + meph['PDObsLon'] + (obs_lon - 180.0)
dec = meph['DEC'] + meph['PDObsLat'] + obs_lat
plt.plot(ra, dec, '.')

plt.figure('point_v_jd')
plt.plot(meph['datetime_jd'], ra)
plt.plot(meph['datetime_jd'], dec)