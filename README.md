This repo contains various definitions and plots for the LFT3 system.

To make WP \label{fig:freqbands}, you can just type (or import) `plot_sens.py`.  It uses pre-computed data generated from:

> import observer
> lft3 = observer.Observe(band='UL')  # one of HF, VL, VH, UL, UH
> lft3.plot_band()

It assumes that you have the pygdsm installed.
Note that using `lft3.run() instead of lft3.plot_band() does sort of the same thing, but gives some different plots and doesn't write the files.

I'll remind myself how to make WP \label{fig:fieldofview} and add that below shortly.