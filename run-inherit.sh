# arguments to the basic simulation code are
# population size, environment, fitness of less good allele, heteroplasmy penalty

gcc -o3 inherit.c -lm -o inherit.ce

# basic simulation, sweeping params for different environments
./inherit.ce 100 0 0.5 0 > tmp &
./inherit.ce 100 1 0.5 0 > tmp &
./inherit.ce 100 2 0.5 0 > tmp &
./inherit.ce 100 3 0.5 0 > tmp &
./inherit.ce 100 4 0.5 0 > tmp &
./inherit.ce 100 5 0.5 0 > tmp &
./inherit.ce 100 6 0.5 0 > tmp &

# smaller organismal population size
./inherit.ce 50 0 0.5 0 > tmp &
./inherit.ce 50 1 0.5 0 > tmp &
./inherit.ce 50 2 0.5 0 > tmp &
./inherit.ce 50 3 0.5 0 > tmp &
./inherit.ce 50 4 0.5 0 > tmp &
./inherit.ce 50 5 0.5 0 > tmp &
./inherit.ce 50 6 0.5 0 > tmp &

# larger organismal population size
./inherit.ce 200 0 0.5 0 > tmp &
./inherit.ce 200 1 0.5 0 > tmp &
./inherit.ce 200 2 0.5 0 > tmp &
./inherit.ce 200 3 0.5 0 > tmp &
./inherit.ce 200 4 0.5 0 > tmp &
./inherit.ce 200 5 0.5 0 > tmp &
./inherit.ce 200 6 0.5 0 > tmp &

# no fitness contribution from less good allele
./inherit.ce 100 0 0 0 > tmp &
./inherit.ce 100 1 0 0 > tmp &
./inherit.ce 100 2 0 0 > tmp &
./inherit.ce 100 3 0 0 > tmp &
./inherit.ce 100 4 0 0 > tmp &
./inherit.ce 100 5 0 0 > tmp &
./inherit.ce 100 6 0 0 > tmp &

# heteroplasmy penalty
./inherit.ce 100 0 0.5 1 > tmp &
./inherit.ce 100 1 0.5 1 > tmp &
./inherit.ce 100 2 0.5 1 > tmp &
./inherit.ce 100 3 0.5 1 > tmp &
./inherit.ce 100 4 0.5 1 > tmp &
./inherit.ce 100 5 0.5 1 > tmp &
./inherit.ce 100 6 0.5 1 > tmp &

# arguments to the baseline simulation code are
# population size, DUI, fitness of less good allele, heteroplasmy penalty
# "baseline" means no leakage or DUI, constant environment

gcc -o3 baseline.c -lm -o baseline.ce

# basic baseline
./baseline.ce 100 0 0 0 > tmp &
./baseline.ce 100 1 0 0 > tmp &
./baseline.ce 100 0 0.5 0 > tmp &
./baseline.ce 100 1 0.5 0 > tmp &
./baseline.ce 100 0 0 1 > tmp &
./baseline.ce 100 1 0 1 > tmp &
./baseline.ce 100 0 0.5 1 > tmp &
./baseline.ce 100 1 0.5 1 > tmp &

# smaller organismal population
./baseline.ce 50 0 0 0 > tmp &
./baseline.ce 50 1 0 0 > tmp &
./baseline.ce 50 0 0.5 0 > tmp &
./baseline.ce 50 1 0.5 0 > tmp &
./baseline.ce 50 0 0 1 > tmp &
./baseline.ce 50 1 0 1 > tmp &
./baseline.ce 50 0 0.5 1 > tmp &
./baseline.ce 50 1 0.5 1 > tmp &

# larger organismal population
./baseline.ce 200 0 0 0 > tmp &
./baseline.ce 200 1 0 0 > tmp &
./baseline.ce 200 0 0.5 0 > tmp &
./baseline.ce 200 1 0.5 0 > tmp &
./baseline.ce 200 0 0 1 > tmp &
./baseline.ce 200 1 0 1 > tmp &
./baseline.ce 200 0 0.5 1 > tmp &
./baseline.ce 200 1 0.5 1 > tmp &

gcc -o3 baseline-single.c -lm -o baseline-single.ce

./baseline-single.ce 100 0 0.5 0 > tmp &

# arguments to the zoomed-in simulation code are
# environment, fitness of less good allele, heteroplasmy penalty

gcc -o3 inherit-zoom.c -lm -o inherit-zoom.ce

./inherit-zoom.ce 0 0.5 0 > tmp &
./inherit-zoom.ce 1 0.5 0 > tmp &
./inherit-zoom.ce 2 0.5 0 > tmp &
./inherit-zoom.ce 3 0.5 0 > tmp &
./inherit-zoom.ce 4 0.5 0 > tmp &
./inherit-zoom.ce 5 0.5 0 > tmp &
./inherit-zoom.ce 6 0.5 0 > tmp &
./inherit-zoom.ce 7 0.5 0 > tmp &
./inherit-zoom.ce 8 0.5 0 > tmp &
./inherit-zoom.ce 9 0.5 0 > tmp &
./inherit-zoom.ce 10 0.5 0 > tmp &
./inherit-zoom.ce 20 0.5 0 > tmp &
./inherit-zoom.ce 50 0.5 0 > tmp &
./inherit-zoom.ce 100 0.5 0 > tmp &
./inherit-zoom.ce 200 0.5 0 > tmp &



