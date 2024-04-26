# odna-inheritance
Different models for oDNA inheritance

Currently (Apr 26th) the repo structure is a bit complicated. `inherit-old.c` takes coarse-grained steps through param space; `inherit-align.c` takes finer steps. The two should otherwise be identical. This is to resolve an apparent inconsistency between older and newer results about the influence of high leakage. The old version shows that high leakage is generally bad; the new version appears not to. Several degrees of freedom could exist. Bug; change from deterministic to stochastic reamplification; specifics about the param values considered. Aim is to reproduce the discrepancy then gradually shift the old to the new protocol to see where the effect vanishes.

Broadly:
 * deterministic leakage (as opposed to stochastic) has a substantial disadvantage to higher leakage values (not sure why?)
 * heteroplasmic initial conditions (as opposed to homoplasmic) preserve heterozygosity much longer and mean the environmental change period is longer
 * deterministic reamplification (as opposed to stochastic) seems to have a smaller effect
