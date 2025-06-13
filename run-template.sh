# master script for organelle DNA inheritance simulation
# requirements: GCC for C compilation

# takes a command-line argument -- a single string, concatenated with commas or other non-whitespace symbols, determining which aspects of the pipeline will be run
# the options are:

# default         -- default model structure
# timeseries      -- example time series
# altmodels       -- alternative model structures (determinism, ICs)
# templaterepair  -- nonzero undirected templated repair
# altrepair       -- nonzero directed templated repair
# popnsize        -- different population sizes
# hetpenalty      -- different heteroplasmy penalties
# bigger          -- heteroplasmy penalty with large population
# fitnessdiffs    -- different fitness differences between alleles
# competition     -- evolutionary competition between strategies
# multilevel      -- within-organism selective differences
# mlrep           -- mutant advantage + template repair
# cluster         -- different inherited cluster sizes

# process command-line arguments
if [ $# -eq 0 ]
then
    echo "No modules selected! Not running anything. See script preamble for options."
    exit 1
fi

commandstr=$1

gcc -o3 inherit-old.c -lm -o inherit-old.ce
gcc -o3 inherit-template.c -lm -o inherit-template.ce
gcc -o3 inherit-template-ts.c -lm -o inherit-template-ts.ce

if [[ $commandstr == *timeseries* ]]; then
  ./inherit-template-ts.ce 10 1 0 0.5 0 0 0 0 > tmp &
  ./inherit-template-ts.ce 10 1 8 0.5 0 0 0 0 > tmp &  
fi

# inherit-old takes coarser steps through parameter space; inherit-template finer
# old arguments: [organismal population size] [homoplasmic ICs] [env change period] [fitness of bad allele] [heteroplasmy penalty] [deterministic reamplification] [deterministic leakage] [for inherit-template, template repair rate]
# output files are labelled by these arguments (and "old" in the case of inherit-old)

if [[ $commandstr == *altmodels* ]]; then
  # compare the various degrees of freedom
  ./inherit-old.ce 100 0 0 0.5 0 0 0 >tmp &
  ./inherit-old.ce 100 0 0 0.5 0 0 1 >tmp &
  ./inherit-old.ce 100 0 0 0.5 0 1 0 >tmp &
  ./inherit-old.ce 100 0 0 0.5 0 1 1 >tmp &
  ./inherit-old.ce 100 1 0 0.5 0 0 0 >tmp &
  ./inherit-old.ce 100 1 0 0.5 0 0 1 >tmp &
  ./inherit-old.ce 100 1 0 0.5 0 1 0 >tmp &
  ./inherit-old.ce 100 1 0 0.5 0 1 1 >tmp &
  
  ./inherit-old.ce 100 0 10 0.5 0 0 0 >tmp &
  ./inherit-old.ce 100 0 10 0.5 0 0 1 >tmp &
  ./inherit-old.ce 100 0 10 0.5 0 1 0 >tmp &
  ./inherit-old.ce 100 0 10 0.5 0 1 1 >tmp &
  ./inherit-old.ce 100 1 10 0.5 0 0 0 >tmp &
  ./inherit-old.ce 100 1 10 0.5 0 0 1 >tmp &
  ./inherit-old.ce 100 1 10 0.5 0 1 0 >tmp &
  ./inherit-old.ce 100 1 10 0.5 0 1 1 >tmp &
fi

# inherit-template arguments:
# [Npop] [initially homoplasmic organisms?] [environmental change period] [fitness scale] [heteroplasmy penalty] [deterministic reamp?] [deterministic leakage?] [templating rate] [rel int fitness of B allele] [rel int fitness of C allele] [cluster size]
if [[ $commandstr == *default* ]]; then
  # default experiment (homoplasmic ICs, stochastic leakage, stochastic reamplification)
  ./inherit-template.ce 100 1 0 0.5 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 2 0.5 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 4 0.5 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 8 0.5 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 16 0.5 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 32 0.5 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 64 0.5 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 128 0.5 0 0 0 0 1 1 1 > tmp &
fi

if [[ $commandstr == *templaterepair* ]]; then
  # default experiment with limited template repair
  ./inherit-template.ce 100 1 0 0.5 0 0 0 0.001 1 1 1 > tmp &
  ./inherit-template.ce 100 1 2 0.5 0 0 0 0.001 1 1 1 > tmp &
  ./inherit-template.ce 100 1 4 0.5 0 0 0 0.001 1 1 1 > tmp &
  ./inherit-template.ce 100 1 8 0.5 0 0 0 0.001 1 1 1 > tmp &
  ./inherit-template.ce 100 1 16 0.5 0 0 0 0.001 1 1 1 > tmp &
  ./inherit-template.ce 100 1 32 0.5 0 0 0 0.001 1 1 1 > tmp &
  ./inherit-template.ce 100 1 64 0.5 0 0 0 0.001 1 1 1 > tmp &
  ./inherit-template.ce 100 1 128 0.5 0 0 0 0.001 1 1 1 > tmp &
  
  # default experiment with stronger template repair
  ./inherit-template.ce 100 1 0 0.5 0 0 0 0.1 1 1 1 > tmp &
  ./inherit-template.ce 100 1 2 0.5 0 0 0 0.1 1 1 1 > tmp &
  ./inherit-template.ce 100 1 4 0.5 0 0 0 0.1 1 1 1 > tmp &
  ./inherit-template.ce 100 1 8 0.5 0 0 0 0.1 1 1 1 > tmp &
  ./inherit-template.ce 100 1 16 0.5 0 0 0 0.1 1 1 1 > tmp &
  ./inherit-template.ce 100 1 32 0.5 0 0 0 0.1 1 1 1 > tmp &
  ./inherit-template.ce 100 1 64 0.5 0 0 0 0.1 1 1 1 > tmp &
  ./inherit-template.ce 100 1 128 0.5 0 0 0 0.1 1 1 1 > tmp &
fi

if [[ $commandstr == *altrepair* ]]; then
  ./inherit-template.ce 100 1 0 0.5 0 0 0 -0.001 1 1 1 > tmp &
  ./inherit-template.ce 100 1 2 0.5 0 0 0 -0.001 1 1 1 > tmp &
  ./inherit-template.ce 100 1 4 0.5 0 0 0 -0.001 1 1 1 > tmp &
  ./inherit-template.ce 100 1 8 0.5 0 0 0 -0.001 1 1 1 > tmp &
  ./inherit-template.ce 100 1 16 0.5 0 0 0 -0.001 1 1 1 > tmp &
  ./inherit-template.ce 100 1 32 0.5 0 0 0 -0.001 1 1 1 > tmp &
  ./inherit-template.ce 100 1 64 0.5 0 0 0 -0.001 1 1 1 > tmp &
  ./inherit-template.ce 100 1 128 0.5 0 0 0 -0.001 1 1 1 > tmp &
fi

if [[ $commandstr == *popnsize* ]]; then
  # smaller population size
  ./inherit-template.ce 50 1 0 0.5 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 50 1 2 0.5 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 50 1 4 0.5 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 50 1 8 0.5 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 50 1 16 0.5 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 50 1 32 0.5 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 50 1 64 0.5 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 50 1 128 0.5 0 0 0 0 1 1 1 > tmp &
  
  # larger population size
  ./inherit-template.ce 200 1 0 0.5 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 200 1 2 0.5 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 200 1 4 0.5 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 200 1 8 0.5 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 200 1 16 0.5 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 200 1 32 0.5 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 200 1 64 0.5 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 200 1 128 0.5 0 0 0 0 1 1 1 > tmp &
fi

if [[ $commandstr == *hetpenalty* ]]; then
  # heteroplasmy penalty
  ./inherit-template.ce 100 1 0 0.5 0.05 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 2 0.5 0.05 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 4 0.5 0.05 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 8 0.5 0.05 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 16 0.5 0.05 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 32 0.5 0.05 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 64 0.5 0.05 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 128 0.5 0.05 0 0 0 1 1 1 > tmp &
  
  # bigger heteroplasmy penalty 
  ./inherit-template.ce 100 1 0 0.5 0.25 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 2 0.5 0.25 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 4 0.5 0.25 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 8 0.5 0.25 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 16 0.5 0.25 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 32 0.5 0.25 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 64 0.5 0.25 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 128 0.5 0.25 0 0 0 1 1 1 > tmp &
fi

if [[ $commandstr == *bigger* ]]; then
  # bigger heteroplasmy penalty, bigger population 
  ./inherit-template.ce 500 1 0 0.5 0.25 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 500 1 2 0.5 0.25 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 500 1 4 0.5 0.25 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 500 1 8 0.5 0.25 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 500 1 16 0.5 0.25 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 500 1 32 0.5 0.25 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 500 1 64 0.5 0.25 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 500 1 128 0.5 0.25 0 0 0 1 1 1 > tmp &
fi

if [[ $commandstr == *fitnessdiffs* ]]; then
  # more comparable mitotype fitness
  ./inherit-template.ce 100 1 0 0.9 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 2 0.9 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 4 0.9 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 8 0.9 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 16 0.9 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 32 0.9 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 64 0.9 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 128 0.9 0 0 0 0 1 1 1 > tmp &
  
  # more distinct mitotype fitness
  ./inherit-template.ce 100 1 0 0.1 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 2 0.1 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 4 0.1 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 8 0.1 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 16 0.1 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 32 0.1 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 64 0.1 0 0 0 0 1 1 1 > tmp &
  ./inherit-template.ce 100 1 128 0.1 0 0 0 0 1 1 1 > tmp &
fi

### post review
# evolutionary experiment -- organisms can have different strategies, which ones are evolutionarily stable when competed under different conditions?
if [[ $commandstr == *competition* ]]; then
  gcc -o3 inherit-comp.c -lm -o inherit-comp.ce
#  ./inherit-comp.ce 1000 0. 0 0 1 1 0 > tmp &
  #  ./inherit-comp.ce 1000 0.5 0 0 1 1 0 > tmp &
  ./inherit-comp.ce 100 0.9 0 0 1 1 0 > tmp &
  ./inherit-comp.ce 1000 0.9 0 0 1 1 0 > tmp &
  gcc -o3 inherit-comp-single.c -lm -o inherit-comp-single.ce
#  ./inherit-comp.ce 1000 0. 0 0 1 1 0 > tmp &
  #  ./inherit-comp.ce 1000 0.5 0 0 1 1 0 > tmp &
  ./inherit-comp-single.ce 100 0.9 0 0 1 1 0 > tmp &
  ./inherit-comp-single.ce 1000 0.9 0 0 1 1 0 > tmp &
fi

# organelle alleles have different amplification rates (within-organism selection)
if [[ $commandstr == *multilevel* ]]; then
gcc -o3 inherit-template.c -lm -o inherit-template.ce
  ./inherit-template.ce 100 1 0 0.5 0 0 0 0 0.91 1 1 > tmp &
  ./inherit-template.ce 100 1 2 0.5 0 0 0 0 0.91 1 1 > tmp &
  ./inherit-template.ce 100 1 4 0.5 0 0 0 0 0.91 1 1 > tmp &
  ./inherit-template.ce 100 1 8 0.5 0 0 0 0 0.91 1 1 > tmp &
  ./inherit-template.ce 100 1 16 0.5 0 0 0 0 0.91 1 1 > tmp &
  ./inherit-template.ce 100 1 32 0.5 0 0 0 0 0.91 1 1 > tmp &
  ./inherit-template.ce 100 1 64 0.5 0 0 0 0 0.91 1 1 > tmp &
  ./inherit-template.ce 100 1 128 0.5 0 0 0 0 0.91 1 1 > tmp &

  ./inherit-template.ce 100 1 0 0.5 0 0 0 0 1.1 1 1 > tmp &
  ./inherit-template.ce 100 1 2 0.5 0 0 0 0 1.1 1 1 > tmp &
  ./inherit-template.ce 100 1 4 0.5 0 0 0 0 1.1 1 1 > tmp &
  ./inherit-template.ce 100 1 8 0.5 0 0 0 0 1.1 1 1 > tmp &
  ./inherit-template.ce 100 1 16 0.5 0 0 0 0 1.1 1 1 > tmp &
  ./inherit-template.ce 100 1 32 0.5 0 0 0 0 1.1 1 1 > tmp &
  ./inherit-template.ce 100 1 64 0.5 0 0 0 0 1.1 1 1 > tmp &
  ./inherit-template.ce 100 1 128 0.5 0 0 0 0 1.1 1 1 > tmp &

  ./inherit-template.ce 100 1 0 0.5 0 0 0 0 1 0.91 1 > tmp &
  ./inherit-template.ce 100 1 2 0.5 0 0 0 0 1 0.91 1 > tmp &
  ./inherit-template.ce 100 1 4 0.5 0 0 0 0 1 0.91 1 > tmp &
  ./inherit-template.ce 100 1 8 0.5 0 0 0 0 1 0.91 1 > tmp &
  ./inherit-template.ce 100 1 16 0.5 0 0 0 0 1 0.91 1 > tmp &
  ./inherit-template.ce 100 1 32 0.5 0 0 0 0 1 0.91 1 > tmp &
  ./inherit-template.ce 100 1 64 0.5 0 0 0 0 1 0.91 1 > tmp &
  ./inherit-template.ce 100 1 128 0.5 0 0 0 0 1 0.91 1 > tmp &

  ./inherit-template.ce 100 1 0 0.5 0 0 0 0 1 1.1 1 > tmp &
  ./inherit-template.ce 100 1 2 0.5 0 0 0 0 1 1.1 1 > tmp &
  ./inherit-template.ce 100 1 4 0.5 0 0 0 0 1 1.1 1 > tmp &
  ./inherit-template.ce 100 1 8 0.5 0 0 0 0 1 1.1 1 > tmp &
  ./inherit-template.ce 100 1 16 0.5 0 0 0 0 1 1.1 1 > tmp &
  ./inherit-template.ce 100 1 32 0.5 0 0 0 0 1 1.1 1 > tmp &
  ./inherit-template.ce 100 1 64 0.5 0 0 0 0 1 1.1 1 > tmp &
  ./inherit-template.ce 100 1 128 0.5 0 0 0 0 1 1.1 1 > tmp &
fi

# organelle alleles have different amplification rates (within-organism selection)
if [[ $commandstr == *mlrep* ]]; then
  ./inherit-template.ce 100 1 0 0.5 0 0 0 0.001 1 1.1 1 > tmp &
  ./inherit-template.ce 100 1 2 0.5 0 0 0 0.001 1 1.1 1 > tmp &
  ./inherit-template.ce 100 1 4 0.5 0 0 0 0.001 1 1.1 1 > tmp &
  ./inherit-template.ce 100 1 8 0.5 0 0 0 0.001 1 1.1 1 > tmp &
  ./inherit-template.ce 100 1 16 0.5 0 0 0 0.001 1 1.1 1 > tmp &
  ./inherit-template.ce 100 1 32 0.5 0 0 0 0.001 1 1.1 1 > tmp &
  ./inherit-template.ce 100 1 64 0.5 0 0 0 0.001 1 1.1 1 > tmp &
  ./inherit-template.ce 100 1 128 0.5 0 0 0 0.001 1 1.1 1 > tmp &
fi

if [[ $commandstr == *cluster* ]]; then
  # inherit pairs
  ./inherit-template.ce 100 1 0 0.5 0 0 0 0 1 1 2 > tmp &
  ./inherit-template.ce 100 1 2 0.5 0 0 0 0 1 1 2 > tmp &
  ./inherit-template.ce 100 1 4 0.5 0 0 0 0 1 1 2 > tmp &
  ./inherit-template.ce 100 1 8 0.5 0 0 0 0 1 1 2 > tmp &
  ./inherit-template.ce 100 1 16 0.5 0 0 0 0 1 1 2 > tmp &
  ./inherit-template.ce 100 1 32 0.5 0 0 0 0 1 1 2 > tmp &
  ./inherit-template.ce 100 1 64 0.5 0 0 0 0 1 1 2 > tmp &
  ./inherit-template.ce 100 1 128 0.5 0 0 0 0 1 1 2 > tmp &
  
  # inherit clusters
  ./inherit-template.ce 100 1 0 0.5 0 0 0 0 1 1 10 > tmp &
  ./inherit-template.ce 100 1 2 0.5 0 0 0 0 1 1 10 > tmp &
  ./inherit-template.ce 100 1 4 0.5 0 0 0 0 1 1 10 > tmp &
  ./inherit-template.ce 100 1 8 0.5 0 0 0 0 1 1 10 > tmp &
  ./inherit-template.ce 100 1 16 0.5 0 0 0 0 1 1 10 > tmp &
  ./inherit-template.ce 100 1 32 0.5 0 0 0 0 1 1 10 > tmp &
  ./inherit-template.ce 100 1 64 0.5 0 0 0 0 1 1 10 > tmp &
  ./inherit-template.ce 100 1 128 0.5 0 0 0 0 1 1 10 > tmp &
fi
