This repo contains various definitions and plots for the LFT3 system.

To make WP \label{fig:freqbands}, you can just type (or import) `plot_sens.py`.  It uses pre-computed data generated from:

```
In [1]: import observer
In [2]: lft3 = observer.Observe(band='UL')  # one of HF, VL, VH, UL, UH
In [3]: lft3.plot_band()
```

It assumes that you have the pygdsm installed.
Note that using `lft3.run() instead of lft3.plot_band() does sort of the same thing, but gives some different plots and doesn't write the files.

I'll remind myself how to make WP \label{fig:fieldofview} and add that below shortly.

fovsrc.py has the general source code.

For tracking speed see fovsrc.sky_track

Most of the base classes are in `fovsrc_Classes.py`:

 - `Position`: very simple class to hold the position of LFT3
 - `Pointing`: determine the RA/Dec of where the telescope is pointing (uses Position)
 - `Sources`: determine the position of various sources (e.g. Sun, Jupiter)
 - `Galaxy`: makes the global sky model and sky temperature values
 - `Observe`: class to help make the Galaxy view

To get the pointing information, you can use this:

```
In [1]: import lft3_sys
In [1]: lft3 = lft3_sys.System([40.0, 560.0, 1420.0])  # List of observing freqs
In [1]: import fovsrc_Classes
In [2]: gal = fovsrc_Classes.Galaxy(lft3.obs_freqs, lft3.obs_fwhm)
In [3]: gal.gen_map_cube()
In [4]: gal.view()
In [5]: gal.gen_pointings()
```