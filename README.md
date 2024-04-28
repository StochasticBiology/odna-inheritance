# odna-inheritance
Different models for oDNA inheritance

Currently (Apr 26th) the repo structure is a bit complicated. `inherit-old.c` takes coarse-grained steps through param space; `inherit-align.c` takes finer steps. The two should otherwise be identical. This is to help rapid investigation of different model structures (old) and in-depth investigation of a chosen protocol (align). 

Some observations from cross-model comparison
 * deterministic leakage (as opposed to stochastic) has a substantial disadvantage to higher leakage values (heterozygosity enforced, rather than a spread of possible values)
 * heteroplasmic initial conditions (as opposed to homoplasmic) preserve heterozygosity much longer and mean the environmental change period supported is longer (or equivalently, performance is better for a given environmental change period)
 * deterministic reamplification (as opposed to stochastic) seems to have a smaller effect
 * DUI has a smaller but detectable positive effect for varying environments

