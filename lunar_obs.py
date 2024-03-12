from astroquery.jplhorizons import Horizons
from pygdsm import GlobalSkyModel
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import healpy as hp


# see Danny Price https://github.com/telegraphic/pygdsm/tree/master
NSIDE = 512

class Galaxy:
    def __init__(self, freqs, fwhm, idisplay):
        self.gsm = GlobalSkyModel(include_cmb=True)  # This uses GSM08 (?)
        self.freqs = freqs
        self.fwhm = fwhm
        self.idisplay = idisplay

    def gen_map_cube(self):
        map_raw = self.gsm.generate(self.freqs)
        if self.fwhm is None:
            self.map_cube = map_raw
        else:
            self.map_cube = []
            for i in range(len(self.freqs)):
                self.map_cube.append(hp.sphtfunc.smoothing(map_raw[i], self.fwhm[i]))
        self.map_cube = np.array(self.map_cube)

    def gen_pointings(self):
        self.coord = Pointing()
        self.locations = hp.pixelfunc.ang2pix(NSIDE, self.coord.midview.galactic.l.to_value(), self.coord.midview.galactic.b.to_value(), lonlat=True)

    def gen_locs(self, step=112000):
        minval = min(self.map_cube[0])
        imin = np.where(self.map_cube[0] < minval + 0.05)[0][0]
        maxval = max(self.map_cube[0])
        imax = np.where(self.map_cube[0] > maxval - 0.05)[0][0]
        self.locations = [imin, imax] + list(range(0, self.map_cube.shape[1], step))

    def view(self, view_fullres=True, logged=True):
        plt.figure('Pointings')
        if view_fullres:
            self.gsm.view(self.idisplay, logged=True)
        else:
            if logged:
                hp.visufunc.mollview(np.log10(self.map_cube[self.idisplay]))
            else:
                hp.visufunc.mollview(self.map_cube[self.idisplay])

class Pointing:
    def __init__(self, observer='Apollo-11 LRRR @ 301', target='399', flip=True, view=30.0,
                 start='2027-07-01T00:00', stop='2027-07-30T23:59', step='240m'):
        #observer = 'g: 180.0, -5.0, 0 @ 301'
        self.target = target
        self.epoch = {'start': start, 'stop': stop, 'step': step}
        self.observer = observer
        apollo_class = Horizons(id=target, location=observer, epochs=self.epoch)
        apollo = apollo_class.ephemerides()
        dt = apollo['datetime_jd'] - apollo['datetime_jd'][0]

        if flip:
            RA = apollo['RA'] + 180.0
            i360 = np.where(RA >= 360)
            if len(i360):
                RA[i360] = RA[i360] - 360
            DEC = -1.0 * apollo['DEC']
        else:
            RA = apollo['RA']
            DEC = apollo['DEC']
        plt.plot(dt, RA / 15.0, '.')
        plt.plot(dt, DEC, '.')
        plt.figure()
        plt.plot(RA / 15.0, DEC, '.')
        plt.xlabel('RA [hr]')
        plt.ylabel('Dec [deg]')

        self.lowview = SkyCoord(ra=np.array(RA)*u.degree, dec=(np.array(DEC) - view/2.0)*u.degree, frame='icrs')
        self.midview = SkyCoord(ra=np.array(RA)*u.degree, dec=np.array(DEC)*u.degree, frame='icrs')
        self.hiview = SkyCoord(ra=np.array(RA)*u.degree, dec=(np.array(DEC) + view/2.0)*u.degree, frame='icrs')

#observer = '500@301'
#luna_class = Horizons(id=target, location=observer, epochs=epoch)
#luna = luna_class.ephemerides()

# plt.figure()
# plt.plot(dt, apollo['RA'])
# plt.plot(dt, apollo['DEC'])
# plt.plot(dt, luna['RA'])
# plt.plot(dt, luna['DEC'])

# plt.figure()
# plt.plot(dt, apollo['RA'] - luna['RA'])
# plt.plot(dt, apollo['DEC'] - luna['DEC'])