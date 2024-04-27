gcc -o3 inherit-old.c -lm -o inherit-old.ce
gcc -o3 inherit-align.c -lm -o inherit-align.ce

# inherit-old takes coarser steps through parameter space; inherit-align finer
# arguments: [organismal population size] [homoplasmic ICs] [env change period] [fitness of bad allele] [heteroplasmy penalty] [deterministic reamplification] [deterministic leakage]
# output files are labelled by these arguments (and "old" in the case of inherit-old)

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

# align with BGP version (heteroplasmic ICs, deterministic leakage, deterministic reamplification)
./inherit-old.ce 100 0 0 0.5 0 1 1 >tmp &
./inherit-old.ce 100 0 5 0.5 0 1 1 >tmp &
./inherit-old.ce 100 0 10 0.5 0 1 1 >tmp &
./inherit-old.ce 100 0 20 0.5 0 1 1 >tmp &
./inherit-old.ce 100 0 30 0.5 0 1 1 >tmp &
./inherit-old.ce 100 0 40 0.5 0 1 1 >tmp &
./inherit-old.ce 100 0 50 0.5 0 1 1 >tmp &
./inherit-old.ce 100 0 60 0.5 0 1 1 >tmp &
./inherit-old.ce 100 0 70 0.5 0 1 1 >tmp &
./inherit-old.ce 100 0 80 0.5 0 1 1 >tmp &
./inherit-old.ce 100 0 90 0.5 0 1 1 >tmp &
./inherit-old.ce 100 0 100 0.5 0 1 1 >tmp &
./inherit-old.ce 100 0 150 0.5 0 1 1 >tmp &
./inherit-old.ce 100 0 200 0.5 0 1 1 >tmp &

# align with stepping-stone version (homoplasmic ICs, deterministic leakage, deterministic reamplification)
./inherit-old.ce 100 1 0 0.5 0 1 1 > tmp &
./inherit-old.ce 100 1 1 0.5 0 1 1 > tmp &
./inherit-old.ce 100 1 2 0.5 0 1 1 > tmp &
./inherit-old.ce 100 1 3 0.5 0 1 1 > tmp &
./inherit-old.ce 100 1 4 0.5 0 1 1 > tmp &
./inherit-old.ce 100 1 5 0.5 0 1 1 > tmp &
./inherit-old.ce 100 1 6 0.5 0 1 1 > tmp &
./inherit-old.ce 100 1 7 0.5 0 1 1 > tmp &
./inherit-old.ce 100 1 8 0.5 0 1 1 > tmp &
./inherit-old.ce 100 1 9 0.5 0 1 1 > tmp &
./inherit-old.ce 100 1 10 0.5 0 1 1 > tmp &
./inherit-old.ce 100 1 15 0.5 0 1 1 > tmp &
./inherit-old.ce 100 1 20 0.5 0 1 1 > tmp &

# default experiment (homoplasmic ICs, stochastic leakage, stochastic reamplification)
./inherit-align.ce 100 1 0 0.5 0 0 0 > tmp &
./inherit-align.ce 100 1 2 0.5 0 0 0 > tmp &
./inherit-align.ce 100 1 4 0.5 0 0 0 > tmp &
./inherit-align.ce 100 1 8 0.5 0 0 0 > tmp &
./inherit-align.ce 100 1 16 0.5 0 0 0 > tmp &
./inherit-align.ce 100 1 32 0.5 0 0 0 > tmp &
./inherit-align.ce 100 1 64 0.5 0 0 0 > tmp &
./inherit-align.ce 100 1 128 0.5 0 0 0 > tmp &

# smaller population size
./inherit-align.ce 50 1 0 0.5 0 0 0 > tmp &
./inherit-align.ce 50 1 2 0.5 0 0 0 > tmp &
./inherit-align.ce 50 1 4 0.5 0 0 0 > tmp &
./inherit-align.ce 50 1 8 0.5 0 0 0 > tmp &
./inherit-align.ce 50 1 16 0.5 0 0 0 > tmp &
./inherit-align.ce 50 1 32 0.5 0 0 0 > tmp &
./inherit-align.ce 50 1 64 0.5 0 0 0 > tmp &
./inherit-align.ce 50 1 128 0.5 0 0 0 > tmp &

# larger population size
./inherit-align.ce 200 1 0 0.5 0 0 0 > tmp &
./inherit-align.ce 200 1 2 0.5 0 0 0 > tmp &
./inherit-align.ce 200 1 4 0.5 0 0 0 > tmp &
./inherit-align.ce 200 1 8 0.5 0 0 0 > tmp &
./inherit-align.ce 200 1 16 0.5 0 0 0 > tmp &
./inherit-align.ce 200 1 32 0.5 0 0 0 > tmp &
./inherit-align.ce 200 1 64 0.5 0 0 0 > tmp &
./inherit-align.ce 200 1 128 0.5 0 0 0 > tmp &

# heteroplasmy penalty
./inherit-align.ce 100 1 0 0.5 0.05 0 0 > tmp &
./inherit-align.ce 100 1 2 0.5 0.05 0 0 > tmp &
./inherit-align.ce 100 1 4 0.5 0.05 0 0 > tmp &
./inherit-align.ce 100 1 8 0.5 0.05 0 0 > tmp &
./inherit-align.ce 100 1 16 0.5 0.05 0 0 > tmp &
./inherit-align.ce 100 1 32 0.5 0.05 0 0 > tmp &
./inherit-align.ce 100 1 64 0.5 0.05 0 0 > tmp &
./inherit-align.ce 100 1 128 0.5 0.05 0 0 > tmp &

# more comparable mitotype fitness
./inherit-align.ce 100 1 0 0.9 0 0 0 > tmp &
./inherit-align.ce 100 1 2 0.9 0 0 0 > tmp &
./inherit-align.ce 100 1 4 0.9 0 0 0 > tmp &
./inherit-align.ce 100 1 8 0.9 0 0 0 > tmp &
./inherit-align.ce 100 1 16 0.9 0 0 0 > tmp &
./inherit-align.ce 100 1 32 0.9 0 0 0 > tmp &
./inherit-align.ce 100 1 64 0.9 0 0 0 > tmp &
./inherit-align.ce 100 1 128 0.9 0 0 0 > tmp &

# more distinct mitotype fitness
./inherit-align.ce 100 1 0 0.1 0 0 0 > tmp &
./inherit-align.ce 100 1 2 0.1 0 0 0 > tmp &
./inherit-align.ce 100 1 4 0.1 0 0 0 > tmp &
./inherit-align.ce 100 1 8 0.1 0 0 0 > tmp &
./inherit-align.ce 100 1 16 0.1 0 0 0 > tmp &
./inherit-align.ce 100 1 32 0.1 0 0 0 > tmp &
./inherit-align.ce 100 1 64 0.1 0 0 0 > tmp &
./inherit-align.ce 100 1 128 0.1 0 0 0 > tmp &
