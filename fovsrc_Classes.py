from astroquery.jplhorizons import Horizons
import numpy as np
from astropy import units as u
import pandas as pd
from astropy.coordinates import SkyCoord
import healpy as hp
import matplotlib.pyplot as plt
from pygdsm import GlobalSkyModel
from copy import copy


DATA_DIR = 'data'


class Position:
    """Keep position parameters: lat, lon, x, y, z, etc"""
    def __init__(self, **kwargs):
        self.r = 1.0  # Make unit sphere unless redefined
        for key, val in kwargs.items():
            setattr(self, key, val)
    
    def latlon2xyz(self):
        self.x = self.r * np.sin(np.deg2rad(90.0 - self.lat)) * np.cos(np.deg2rad(self.lon))
        self.y = self.r * np.sin(np.deg2rad(90.0 - self.lat)) * np.sin(np.deg2rad(self.lon))
        self.z = self.r * np.cos(np.deg2rad(90.0 - self.lat))

LFT3 = Position(lon=182.13737, lat=-23.78930)


def get_sublunar(start_time, end_time, step):
    """Return the sublunar point"""
    mobj = Horizons(id=301, location=500, epochs={'start': start_time, 'stop': end_time, 'step': step})
    meph = mobj.ephemerides()
    sublunar = Position(lon=meph['PDObsLon'], lat = meph['PDObsLat'], times=meph['datetime_jd'])
    for i in range(0, len(sublunar.lon)):
        if sublunar.lon[i] > 180:
            sublunar.lon[i] -= 360
    return sublunar

def angles(tele, origin):
    offset = Position(ra=[], dec=[])
    for i in range(len(origin.lat)):
        anti = Position(lon=origin.lon[i] + 180.0, lat=-1.0 * origin.lat[i])
        if anti.lon >= 360.0:
            anti.lon -= 360.0
        offset.ra.append(tele.lon - anti.lon)
        offset.dec.append(tele.lat - anti.lat)
    return offset

def overpole(rd, resample_type=None):
    if resample_type is not None:
        RES = 1.0  # 1deg
        RAresample = np.arange(0, 360.0, RES)
        kernel = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1])
        kernel = kernel / np.sum(kernel)
    lo_flip = np.where(rd.DEC < -90.0)
    if len(lo_flip[0]):
        lowrap = True
        rd.DEC[lo_flip] = -1.0 * (rd.DEC[lo_flip] + 180.0)
        rd.RA[lo_flip] = rd.RA[lo_flip] + 180.0
        # ra_wrap should equal lo_flip
        ra_wrap = np.where(rd.RA > 360.0)
        rd.RA[ra_wrap] = rd.RA[ra_wrap] - 360.0
        ret = Position(RA=rd.RA, DEC=rd.DEC)
    else:
        lowrap = False
        ret = rd
    # Don't have for now
    # hi_flip = np.where(rd.DEC > 90.0)
    # rd.DEC[hi_flip] = 180.0 - rd.DEC[hi_flip]
    # rd.RA[hi_flip] = rd.RA[hi_flip] + 180.0
    # ra_wrap = np.where(rd.RA > 360.0)
    # rd.RA[ra_wrap] = rd.RA[ra_wrap] - 360.0
    if resample_type is not None:
        DECresample = []
        rdchk = (ret.RA / RES).astype(int)
        for rsam in RAresample:
            key = int(rsam / RES)
            entries = np.where(rdchk == key)
            if len(entries[0]) > 1:
                DECresample.append(resample_type(ret.DEC[entries]))
            else:  # this is pretty ad hoc...
                if lowrap:
                    DECresample.append(-90.0)
                else:
                    DECresample.append(np.average(ret.DEC))
                    print("ADHOC")
        DECrepl = DECresample[-50:] + DECresample + DECresample[:50]
        DECresample = np.convolve(DECrepl, kernel, mode='same')[50:-50]
        ret = Position(RA=RAresample, DEC=DECresample)
    return ret

class Pointing:
    """Get where the telescope will be pointing"""
    def __init__(self, telescope_lon, telescope_lat, start_time, end_time, step, show_plots=True):
        self.telescope_position = Position(lon=telescope_lon, lat=telescope_lat)
        telescope = {'lon': telescope_lon * u.deg, 'lat': telescope_lat * u.deg, 'elevation': 0 * u.km, 'body': 301}
        self.telescope = telescope
        self.start_time = start_time
        self.end_time = end_time
        self.step = step
        self.get_zenith_pointing()
        # Apply field-of-view extents
        lo = self.fov_offset(self.zen, DEC=-60)
        hi = self.fov_offset(self.zen, DEC=60)
        if show_plots:
            plt.figure('RA')
            plt.plot(self.times, self.zen.RA, '.')
            if len(self.times) == len(lo.RA):
                plt.plot(self.times, lo.RA, '.')
                plt.plot(self.times, hi.RA, '.')
            plt.figure('DEC')
            plt.plot(self.times, self.zen.DEC, '.')
            if len(self.times) == len(lo.DEC):
                plt.plot(self.times, lo.DEC, '.')
                plt.plot(self.times, hi.DEC, '.')
            plt.figure("RADEC")
            plt.plot(self.zen.RA, self.zen.DEC, '.')
            plt.plot(lo.RA, lo.DEC, '.')
            plt.plot(hi.RA, hi.DEC, '.')

        self.zenith = SkyCoord(ra=np.array(self.zen.RA) * u.degree, dec=np.array(self.zen.DEC) * u.degree, frame='icrs')
        self.loview = SkyCoord(ra=np.array(lo.RA) * u.degree, dec=np.array(lo.DEC) * u.degree, frame='icrs')
        self.hiview = SkyCoord(ra=np.array(hi.RA) * u.degree, dec=np.array(hi.DEC) * u.degree, frame='icrs')

    def get_zenith_pointing(self):
        # Get sublunar point for all times
        self.sublunar = get_sublunar(self.start_time, self.end_time, self.step)
        # Get pointing to the center of the moon for all times
        mobj = Horizons(id=301, location=500, epochs={'start': self.start_time, 'stop': self.end_time, 'step': self.step})
        meph = mobj.ephemerides()
        self.times = meph['datetime_jd']
        # Get offset angles for sub-lunar point and telescope location
        offset = angles(self.telescope_position, self.sublunar)
        # Apply offsets and correct for wraps
        self.zen = Position(RA=[], DEC=[])
        for rac, decc, rad, decd in zip(meph['RA'], meph['DEC'], offset.ra, offset.dec):
            self.zen.RA.append(rac + rad)
            self.zen.DEC.append(decc + decd)
        self.zen.RA = np.array(self.zen.RA)
        rawrap = np.where(self.zen.RA >= 360.0)
        self.zen.RA[rawrap] = self.zen.RA[rawrap] - 360.0
        self.zen.DEC = np.array(self.zen.DEC)
        # Don't need for now
        #self.zen = overpole(self.zen, resample=False)

    def fov_offset(self, orig, RA=0.0, DEC=0.0):
        newpos = Position(RA=orig.RA + RA, DEC=orig.DEC + DEC)
        if DEC < 0.0:
            resample_type = min
        else:
            resample_type = max
        newpos = overpole(newpos, resample_type=resample_type)
        return newpos
    
class Sources:
    def __init__(self, start_time, end_time, step, sources=['sun', 'jupiter', 'mars'], method='compute&write', show_plots=False):
        self.beam = Pointing(LFT3.lon, LFT3.lat, start_time, end_time, step, show_plots)
        if method.startswith('compute'):
            print(f"Computing positions for {', '.join(sources)}")
            beam = {}
            for i in range(len(self.beam.times)):
                beam[int(1440.0 * self.beam.times[i])] = (self.beam.zen.RA[i], self.beam.zen.DEC[i])
            if 'sun' in sources:
                sun = Horizons(id='sun', location=500, epochs={'start': start_time, 'stop': end_time, 'step': step})
                self.sun = self.get_visible(beam, sun)
                self.sunview = SkyCoord(ra=np.array(self.sun.RA) * u.degree, dec=np.array(self.sun.DEC) * u.degree, frame='icrs')
            if 'jupiter' in sources:
                jup = Horizons(id=599, location=500, epochs={'start': start_time, 'stop': end_time, 'step': step})
                self.jupiter = self.get_visible(beam, jup)
                self.jupiterview = SkyCoord(ra=np.array(self.jupiter.RA) * u.degree, dec=np.array(self.jupiter.DEC) * u.degree, frame='icrs')
            if 'mars' in sources:
                mars = Horizons(id=299, location=500, epochs={'start': start_time, 'stop': end_time, 'step': step})
                self.mars = self.get_visible(beam, mars)
                self.marsview = SkyCoord(ra=np.array(self.mars.RA) * u.degree, dec=np.array(self.mars.DEC) * u.degree, frame='icrs')
            if 'alphacen' in sources:
                self.alphacen = self.get_visible(beam, (219.90206, -60.83399))
                self.alphacenview = SkyCoord(ra=219.90206 * u.degree, dec=-60.83399 * u.degree, frame='icrs')
            if 'southpole' in sources:
                self.southpole = Position(RA=[0.0], DEC=[-89.99], vis_times=[0.0])
                self.southpoleview = SkyCoord(ra=0.0 * u.degree, dec=-89.9 * u.degree, frame='icrs')
            if 'stars' in sources:
                #stars2chk = pd.read_csv(f"{DATA_DIR}/hundred_parsecs-result.csv", dtype=None)
                stars2chk = pd.read_fwf(f"{DATA_DIR}/closest_stars.txt", dtype=None)
                self.stars = Position(RA=[], DEC=[], vis_times=[])
                print(f"Using {len(stars2chk)} stars")
                for i in range(len(stars2chk)):
                    thisone = (stars2chk['RA'][i], stars2chk['DEC'][i])
                    if self.get_visible(beam, thisone, check_only=True):
                        self.stars.RA.append(thisone[0])
                        self.stars.DEC.append(thisone[1])
                        self.stars.vis_times.append(0.0)
                self.starsview = SkyCoord(ra=np.array(self.stars.RA) * u.degree, dec=np.array(self.stars.DEC) * u.degree, frame='icrs')

            if method.endswith('write'):
                print("Writing positions")
                for src in sources:
                    with open("{src}.dat", 'w') as fp:
                        sv = getattr(self, src)
                        print(f"#jd,ra,dec  --  {start_time}-{end_time} step {step}", file=fp)
                        for _t, _r, _d in zip(sv.vis_times, sv.RA, sv.DEC):
                            print(f"{_t},{_r},{_d}", file=fp)
        elif method.startswith('read'):
            print(f"Reading positions for {', '.join(sources)}")
            for src in sources:
                vis_times, ra, dec = [], [], []
                fn = f"{DATA_DIR}/{src}.dat"
                print(f"Reading {fn}")
                with open(fn, 'r') as fp:
                    for line in fp:
                        if line[0] == '#':
                            continue
                        data = [float(x) for x in line.strip().split(',')]
                        vis_times.append(data[0])
                        ra.append(data[1])
                        dec.append(data[2])
                    setattr(self, src, Position(vis_times=vis_times, RA=ra, DEC=dec))
                    setattr(self, f"{src}view", SkyCoord(ra=np.array(ra) * u.degree, dec=np.array(dec) * u.degree, frame='icrs'))

    def get_visible(self, beam, body, check_only=False):
        vbod = Position(RA=[], DEC=[], vis_times=[])
        if isinstance(body, (list, tuple)):
            chk_is_single = True
            chk = body  # List/tuple of one RA/Dec
        else:
            chk_is_single = False
            chk = {}
            for i, jd in enumerate(body.ephemerides()['datetime_jd']):
                chk[int(jd * 1440.0)] = (body.ephemerides()['RA'][i], body.ephemerides()['DEC'][i])
        #print("This should use lo_RADEC/hi_RADEC...")
        for key in beam:
            try:
                if chk_is_single:
                    sRA = chk[0]
                    sDEC = chk[1]
                else:
                    sRA = chk[key][0]
                    sDEC = chk[key][1]
                bRA = beam[key][0]
                bDEC = beam[key][1]
                if (bRA - 10.0) <= sRA <= (bRA + 10.0):
                    nearestRA = np.abs(self.beam.loview.ra.value - bRA).argmin()
                    #nearestRA = 1
                    loDEC = bDEC - 60.0
                    loDEC = self.beam.loview.dec[nearestRA].value
                    hiDEC = bDEC + 60.0
                    hiDEC = self.beam.hiview.dec[nearestRA].value
                    if loDEC <= sDEC <= hiDEC:
                        if check_only:
                            return True
                        vbod.vis_times.append(key / 1440.0)
                        vbod.RA.append(sRA)
                        vbod.DEC.append(sDEC)
            except KeyError:
                continue
        if check_only:
            return False
        return vbod

# lunar_obs
class Galaxy:
    def __init__(self, freqs, fwhm=None):
        self.gsm = GlobalSkyModel(include_cmb=True)  # This uses GSM08 (?)
        self.freqs = freqs
        if fwhm is None:
            self.fwhm = None
        elif type(fwhm) == float:
            self.fwhm = np.ones(len(freqs)) * fwhm
        else:
            self.fwhm = fwhm
        self.idisplay = int((freqs[-1] - freqs[0]) / (freqs[1] - freqs[0]))

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
        self.coord = Pointing(LFT3.lon, LFT3.lat, '2028-01-01T00:00:00', '2028-12-31T23:59:59', '1h')
        if self.coord.zenith is not None:
            self.locations = hp.pixelfunc.ang2pix(NSIDE, theta=self.coord.zenith.galactic.l.to_value(),
                                              phi=self.coord.zenith.galactic.b.to_value(), lonlat=True)

    def view(self, view_fullres=True, logged=True):
        if view_fullres:
            self.gsm.view(self.idisplay, logged=True)
        else:
            if logged:
                hp.visufunc.mollview(np.log10(self.map_cube[self.idisplay]))
            else:
                hp.visufunc.mollview(self.map_cube[self.idisplay])

# observer
NSIDE = 512
plt.style.use('ggplot')

color_palette = [
    (0.12156862745098039, 0.4666666666666667, 0.7058823529411765, 1.0),
    (1.0, 0.4980392156862745, 0.054901960784313725, 1.0),
    (0.17254901960784313, 0.6274509803921569, 0.17254901960784313, 1.0),
    (0.8392156862745098, 0.15294117647058825, 0.1568627450980392, 1.0),
    (0.5803921568627451, 0.403921568627451, 0.7411764705882353, 1.0),
    (0.5490196078431373, 0.33725490196078434, 0.29411764705882354, 1.0)
]

class Observe:
    def __init__(self, start=250.0, stop=750, step=100):
        self.freqs = np.arange(start, stop, step)
        self.fwhm = 10.0
        self.galaxy = None

    def get_galaxy(self):
        self.galaxy = Galaxy(freqs=self.freqs, fwhm=self.fwhm)
        self.galaxy.gen_map_cube()
        self.galaxy.gen_pointings()
        self.galaxy.view()

    def get_sky_Tsys(self, Trcvr=30.0):
        if self.galaxy is None:
            self.get_galaxy()
        self.Gal = []
        self.Tsys = []
        for i in range(len(self.galaxy.locations)):
            self.Gal.append(self.galaxy.map_cube[:, self.galaxy.locations[i]])
            self.Tsys.append(self.galaxy.map_cube[:, self.galaxy.locations[i]] + Trcvr)

    def run(self, start='2028-01-01T00:00:00', stop='2028-12-31T23:59:59', step='2h', sources=['sun', 'jupiter', 'mars'], method='read'):
        self.get_sky_Tsys()
        x = Sources(start, stop, step, sources, method, show_plots=False)
        self.source = x
        #hp.visufunc.projplot(self.galaxy.coord.zenith.galactic.l.to_value(), self.galaxy.coord.zenith.galactic.b.to_value(),'k', linestyle=':', lonlat=True)
        hp.visufunc.projplot(self.galaxy.coord.loview.galactic.l.to_value(), self.galaxy.coord.loview.galactic.b.to_value(), 'k', lonlat=True)
        hp.visufunc.projplot(self.galaxy.coord.hiview.galactic.l.to_value(), self.galaxy.coord.hiview.galactic.b.to_value(), 'k', lonlat=True)
        if 'stars' in sources:
            hp.visufunc.projplot(x.starsview.galactic.l.to_value(), x.starsview.galactic.b.to_value(), 'w*', markersize=3, lonlat=True)
        if 'sun' in sources:
            hp.visufunc.projplot(x.sunview.galactic.l.to_value(), x.sunview.galactic.b.to_value(), 'yo', lonlat=True)
        if 'jupiter' in sources:
            hp.visufunc.projplot(x.jupiterview.galactic.l.to_value(), x.jupiterview.galactic.b.to_value(), 'co', lonlat=True)
        if 'mars' in sources:
            hp.visufunc.projplot(x.marsview.galactic.l.to_value(), x.marsview.galactic.b.to_value(), 'ro', lonlat=True)
        if 'alphacen' in sources:
            hp.visufunc.projplot(x.alphacenview.galactic.l.to_value(), x.alphacenview.galactic.b.to_value(), 'bo', lonlat=True)
        if 'southpole' in sources:
            hp.visufunc.projplot(x.southpoleview.galactic.l.to_value(), x.southpoleview.galactic.b.to_value(), 'rs', lonlat=True)