import matplotlib.pyplot as plt
import fovsrc_Classes

# First run get_source_visibility with compute&write to get the appropriate files.
# Move it/them to source_data
# Then run galaxy
START_AT = "2028-01-01T00:00:00"
STOP_AT = "2028-12-31T23:59:59"
STEP = '2h'

def sublunar(start_at=START_AT, stop_at=STOP_AT, step=STEP, write=False):
    sublunar = fovsrc_Classes.get_sublunar(start_at, stop_at, step)
    plt.figure('sublunar')
    plt.plot(sublunar.lon, sublunar.lat, '.', color='r')
    plt.xlabel('Selenographic longitude')
    plt.ylabel('Selenographic latitude')
    plt.title(f'Sublunar point {start_at} - {stop_at}')
    if write:
        print("Writing sublunar.dat")
        with open('sublunar.dat', 'w') as fp:
            for t, o, a in zip(sublunar.times, sublunar.lon, sublunar.lat):
                print(f"{t},{o},{a}", file=fp)
    return sublunar

def pointing(start_at=START_AT, stop_at=STOP_AT, step=STEP):
    tele = fovsrc_Classes.Pointing(fovsrc_Classes.LFT3.lon, fovsrc_Classes.LFT3.lat, start_at, stop_at, step)
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
    ax1.plot(tele.times, tele.zen.RA, '.')
    ax2.plot(tele.times, tele.zen.DEC)
    ax3.plot(tele.zen.RA, tele.zen.DEC, '.')
    ax1.set_title('Pointing RA vs time')
    ax2.set_title('Pointing Dec vs time')
    ax3.set_title('Pointing RA vs Dec')
    ax1.set_xlabel('Time (jd)')
    ax1.set_ylabel('RA (degrees)')
    ax2.set_xlabel('Time (jd)')
    ax2.set_ylabel('Dec (degrees)')
    ax3.set_xlabel('Dec (degrees)')
    ax3.set_ylabel('RA (degrees)')
    fig.suptitle('2028')

def get_source_visibility(start_at=START_AT, stop_at=STOP_AT, step=STEP, method='compute&write', sources=['stars']):
    source = fovsrc_Classes.Sources(start_at, stop_at, step, sources=sources, method=method)
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
    if 'sun' in sources:
        ax1.plot(source.sun.vis_times, source.sun.RA, '.', color='y', label='Sun')
        ax2.plot(source.sun.vis_times, source.sun.DEC, '.', color='y', label='Sun')
        ax3.plot(source.sun.RA, source.sun.DEC, '.', color='y', label='Sun')
    if 'jupiter' in sources:
        ax1.plot(source.jupiter.vis_times, source.jupiter.RA, '.', color='b', label='Jupiter')
        ax2.plot(source.jupiter.vis_times, source.jupiter.DEC, '.', color='b', label='Jupiter')
        ax3.plot(source.jupiter.RA, source.jupiter.DEC, '.', color='b', label='Jupiter')
    if 'mars' in sources:
        ax1.plot(source.mars.vis_times, source.mars.RA, '.', color='r', label='Mars')
        ax2.plot(source.mars.vis_times, source.mars.DEC, '.', color='r', label='Mars')
        ax3.plot(source.mars.RA, source.mars.DEC, '.', color='r', label='Mars')
    if 'alphacen' in sources:
        ax1.plot(source.alphacen.vis_times, source.alphacen.RA, '.', color='y', label='Sun')
        ax2.plot(source.alphacen.vis_times, source.alphacen.DEC, '.', color='y', label='Sun')
        ax3.plot(source.alphacen.RA, source.alphacen.DEC, '.', color='y', label='Sun')
    ax1.set_title('RA vs time')
    ax2.set_title('Dec vs time')
    ax3.set_title('RA vs Dec')
    ax1.set_xlabel('Time (jd)')
    ax1.set_ylabel('RA (degrees)')
    ax2.set_xlabel('Time (jd)')
    ax2.set_ylabel('Dec (degrees)')
    ax3.set_xlabel('Dec (degrees)')
    ax3.set_ylabel('RA (degrees)')
    #fig.suptitle('Sun 6/2028 through 8/2028')
    ax1.legend()

def galaxy(start_at=START_AT, stop_at=STOP_AT, step=STEP, sources=['sun', 'jupiter', 'mars', 'alphacen', 'stars'], method='read'):
    obs = fovsrc_Classes.Observe()
    obs.run(start_at, stop_at, step, sources=sources, method=method)


from astropy.coordinates import SkyCoord
import astropy.units as u

def angular_separation(ra1_deg, dec1_deg, ra2_deg, dec2_deg):
    """
    Compute the angular separation between two celestial coordinates.

    Parameters:
    - ra1_deg, dec1_deg: Right Ascension and Declination of the first object (in degrees)
    - ra2_deg, dec2_deg: Right Ascension and Declination of the second object (in degrees)

    Returns:
    - Angular separation in degrees.
    """
    coord1 = SkyCoord(ra=ra1_deg * u.deg, dec=dec1_deg * u.deg, frame='icrs')
    coord2 = SkyCoord(ra=ra2_deg * u.deg, dec=dec2_deg * u.deg, frame='icrs')
    
    separation = coord1.separation(coord2)
    return separation.degree

START_AT = "2028-12-01T00:00:00"
STOP_AT = "2028-12-31T23:59:59"
STEP = '1m'
SIDEREAL = 15.041 / 3600.0 
def sky_track(start_at=START_AT, stop_at=STOP_AT, step=STEP):
    filename = f"sky_track_{START_AT.split('-')[1]}.dat"
    ptg = fovsrc_Classes.Pointing(fovsrc_Classes.LFT3.lon, fovsrc_Classes.LFT3.lat, start_at, stop_at, step)
    rel_time = (ptg.times - ptg.times[0]) * 3600.0 * 24.0
    rate = [0.0]
    frac = [0.0]
    for i in range(len(ptg.zen.RA)):
        try:
            dist = angular_separation(ptg.zen.RA[i+1], ptg.zen.DEC[i+1], ptg.zen.RA[i], ptg.zen.DEC[i])
            this_rate = dist / (rel_time[i+1] - rel_time[i])  # in degrees per second
            rate.append(this_rate)
            this_frac = SIDEREAL / this_rate if this_rate != 0 else 0
            frac.append(this_frac)
        except Exception as e:
            continue
    rate[0] = rate[1]  # Set the first point to the first distance
    frac[0] = frac[1]  # Set the first point to the first fraction
    with open(filename, 'w') as fp:
        for t, r, f in zip(rel_time, rate, frac):
            print(f"{t},{r},{f}", file=fp)
    print(f"Sky track data written to {filename}")


    plt.figure('Sky Track')
    plt.plot(rel_time, rate, '.')
    plt.xlabel('Relative Time (seconds)')
    plt.ylabel('Deg/sec')
    plt.title(f'Sky Track from {start_at} to {stop_at}')

    plt.plot([0, rel_time[-1]], [SIDEREAL, SIDEREAL], 'r--', label='Sidereal Rate')
    plt.grid()

    plt.figure('Sky Track Fraction')
    plt.plot(rel_time, frac, '.')
    plt.xlabel('Relative Time (seconds)')
    plt.ylabel('1 / Fraction of Sidereal Rate')
    plt.title(f'Sky Track Fraction from {start_at} to {stop_at}')
    return ptg
