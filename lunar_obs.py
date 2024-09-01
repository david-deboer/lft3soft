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
        #plt.figure('Pointings')
        if view_fullres:
            print("FULL")
            self.gsm.view(self.idisplay, logged=True)
        else:
            if logged:
                print("LOGGED")
                hp.visufunc.mollview(np.log10(self.map_cube[self.idisplay]))
            else:
                print("ELSE")
                hp.visufunc.mollview(self.map_cube[self.idisplay])

class Pointing:
    def __init__(self, view=30.0, obs_lon=182, obs_lat=-23.0,
                 start='2028-07-01T00:00', stop='2028-12-31T23:59', step='1d'):

        mobj = Horizons(id=301, location=500, epochs={'start': start, 'stop': stop, 'step': step})
        meph = mobj.ephemerides()
        w = np.where(meph['PDObsLon'] > 180.0)
        meph['PDObsLon'][w] = meph['PDObsLon'][w] - 360.0
        RA = meph['RA'] + meph['PDObsLon'] + (obs_lon - 180.0)
        DEC = meph['DEC'] + meph['PDObsLat'] + obs_lat

        plt.plot(meph['datetime_jd'], RA / 15.0, '.')
        plt.plot(meph['datetime_jd'], DEC, '.')
        this_fig = plt.figure('radec mw')
        ax = this_fig.add_subplot(111)
        ax.plot(RA / 15.0, DEC, '.')
        gp_l = np.arange(0.0, 359.9, 0.1) * u.deg
        gp_b = np.ones(len(gp_l)) * 0.0 * u.deg
        gal_plane = SkyCoord(frame='galactic', l=gp_l, b=gp_b)
        radec = gal_plane.transform_to('icrs')
        ax.plot(radec.ra.to_value() / 15.0, radec.dec.to_value(), '.', label='GalPlane')
        gal_cen = SkyCoord(frame='galactic', l=0.0 * u.deg, b=0.0 * u.deg)
        radecc = gal_cen.transform_to('icrs')
        ax.plot(radecc.ra.to_value() / 15.0, radecc.dec.to_value(), 'o', markersize=8)
        plt.xlabel('RA [hr]')
        plt.ylabel('Dec [deg]')

        next_fig = plt.figure('galactic')
        next_ax = next_fig.add_subplot(111, projection='mollweide')
        pts = SkyCoord(frame='icrs', ra=RA * u.deg, dec = DEC * u.deg)
        lb = pts.transform_to('galactic')
        cl = lb.l.to_value()
        w = np.where(cl > 180.0)
        cl[w] = cl[w] - 360.0
        next_ax.plot(-1.0 * cl, lb.b.to_value(), 'o')

        self.lowview = SkyCoord(ra=np.array(RA)*u.degree, dec=(np.array(DEC) - view/2.0)*u.degree, frame='icrs')
        self.midview = SkyCoord(ra=np.array(RA)*u.degree, dec=np.array(DEC)*u.degree, frame='icrs')
        self.hiview = SkyCoord(ra=np.array(RA)*u.degree, dec=(np.array(DEC) + view/2.0)*u.degree, frame='icrs')

        next_ax.plot(self.midview.galactic.l.to_value(), self.midview.galactic.b.to_value(), 'o')

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