# odna-inheritance
Different models for oDNA inheritance

Currently, `inherit.c` scans through different environments and leakage / DUI combinations, reporting a limited set of fitness time series and complete set of fitness summary statistics for each. `baseline.c` just considers two environments, no leakage or DUI, and reports the same after a finer-grained scan through oDNA population size and mutation rate (to characterise baseline mutation-selection balance).
