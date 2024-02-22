# odna-inheritance
Different models for oDNA inheritance

Currently, `inherit.c` scans through different environments and leakage / DUI combinations, reporting a limited set of fitness time series and complete set of fitness summary statistics for each. `baseline.c` just considers two environments, no leakage or DUI, and reports the same after a finer-grained scan through oDNA population size and mutation rate (to characterise baseline mutation-selection balance); `baseline-single.c` does this for a single oDNA per cell, seeking to match simpler population genetic theory. `inherit-zoom.c` zooms in to a region of parameter space where re-entrant behaviour seems to exist.

`run-inherit.sh` is a wrapper script compiling and running these simulations with various parameter ranges corresponding to the investigation. `plotter.R` is R script visualising the output and doing some statistics for empirical model fits.
