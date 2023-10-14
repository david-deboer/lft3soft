from astroquery.jplhorizons import Horizons
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord


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