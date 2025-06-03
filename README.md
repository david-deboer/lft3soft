This repo contains various definitions and plots for the LFT3 system.

`lft3_sys.py` contains the system description and some computed estimated parameters.

`observer.py` contains scripts to read in the lft3_sys and `lunar_obs.py`, which has the pointing direction and galaxy etc and produces various figures.

To initiate:

```
import observer
lft3 = observer.Observe(band='UL')
```
where band is in `[HF, VL, VH, UL, UH]`.

To generate the sensitivity information use:
`lft3.plot_bands()`

lft3.run()

To make the frequency band/sensitivity plot using the data saved from the `plot_bands` module, use `plot_sens.py`