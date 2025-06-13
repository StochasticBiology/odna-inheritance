# odna-inheritance
Simulation investigation of (multifaceted organelle inheritance strategies) x (changing environments) x (mutational hazard).

![image](https://github.com/user-attachments/assets/5bb06499-700c-4b79-97a3-f4f84fc5db4f)

Requirements
----
For the simulation code, only GCC is needed to compile C code. The wrapper script is written for Bash but you can easily see the different calls it makes to run the code.

For visualisation, R with `reshape2`, `dplyr`, `ggplot2`, `metR`, `viridis`, and `ggpubr` libraries is required. 

Outline
----

The simulation code is split into several submodules for tractability. Each submodule can be invoked by passing its name as part of a command-line argument to `run-template.sh`. (To run `run-template.sh` from the command line, you will probably need to mark it as executable, with e.g. `chmod +x run-template.sh`). The command-line argument should be a single string, concatenated with commas or other non-whitespace symbols, determining which aspects of the pipeline will be run. For example,

`./run-template.sh default,timeseries` 

would run the first two submodules listed below.

The options are:

* `default`         -- default model structure
* `timeseries`      -- example time series
* `altmodels`       -- alternative model structures (determinism, ICs)
* `templaterepair`  -- nonzero undirect templated repair
* `altrepair`       -- nonzero directed templated repair
* `popnsize`        -- different population sizes
* `hetpenalty`      -- different heteroplasmy penalties
* `bigger`          -- heteroplasmy penalty with large population
* `fitnessdiffs`    -- different fitness differences between alleles
* `competition`     -- evolutionary competition between strategies
* `multilevel`      -- within-organism selective differences
* `mlrep`           -- mutant advantage + template repair
* `cluster`         -- different inherited cluster sizes

`timeseries` takes a few seconds; each of the others will probably take several dozen core-days (each runs one or more sets of about ten parallel simulations, most of which take several hours or so).

Following the simulation code, `plot-template.R` plots most of the figures and performs some statistical analysis; `plot-comp.R` plots summary outputs of the evolutionary competition experiments.

Details
----

`inherit-old.c` takes coarse-grained steps through param space for different model variants; `inherit-template.c` takes finer steps, and `inherit-template-ts.c` records specific time series behaviour. `inherit-comp.c` and `inherit-comp-single.c` run the competition experiments where different strategies evolve together, respectively for a repeatedly-switching and once-switching environmental change regime. 

The bulk of the work is done by `inherit-template.c`. This has command-line structure (if the executable is called `inherit-template.ce`)

`./inherit-template.ce [Npop] [initially homoplasmic organisms?] [environmental change period] [fitness scale] [heteroplasmy penalty] [deterministic reamp?] [deterministic leakage?] [templating rate] [rel int fitness of B allele] [rel int fitness of C allele] [cluster size]`

So that, for example,

`./inherit-template.ce 100 1 32 0.5 0 0 0 0.001 1 1 1`

Runs a simulation with 100 initially homoplasmic individuals, environment change after 32 generations, non-optimal allele contributing 0.5 to fitness, zero heteroplasmy penalty, stochastic reamplification and leakage, undirected repair rate 0.001, equal fitnesses of A, B, and C alleles, and individual oDNAs (cluster size 1) inherited.

Awkwardly, to distinguish *directed* repair, we use the templating rate argument but with a negative value:

`./inherit-template.ce 100 1 32 0.5 0 0 0 -0.001 1 1 1`

This is a legacy code issue and will be refactored in future. Some other legacy code is included due to the history of the project: `Inheritance_algorithm.py` is Python code for the same type of simulation, but works rather more slowly. `Inherit_comparison.py` is this code looking at a subset of parameter space, for comparison with the C code. 
