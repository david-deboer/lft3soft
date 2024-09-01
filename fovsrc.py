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