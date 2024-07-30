# odna-inheritance
Different models for oDNA inheritance

`inherit-old.c` takes coarse-grained steps through param space for different model variants; `inherit-template.c` takes finer steps, and `inherit-template-ts.c' records specific time series behaviour. 

`plot-template.R` plots everything; `run-template.sh` is a Bash script running the C simulations.

Some observations from cross-model comparison
 * deterministic leakage (as opposed to stochastic) has a substantial disadvantage to higher leakage values (heterozygosity enforced, rather than a spread of possible values)
 * heteroplasmic initial conditions (as opposed to homoplasmic) preserve heterozygosity much longer and mean the environmental change period supported is longer (or equivalently, performance is better for a given environmental change period)
 * deterministic reamplification (as opposed to stochastic) seems to have a smaller effect
 * DUI has a smaller but detectable positive effect for varying environments

`Inheritance_algorithm.py` is Python code for the same type of simulation, but works rather more slowly. `Inherit_comparison.py` is this code looking at a subset of parameter space, for comparison with the C code. 
