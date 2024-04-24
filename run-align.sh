gcc -o3 inherit-align.c -lm -o inherit-align.ce

# align with BGP version
./inherit-align.ce 100 0 0 0.5 0 > tmp &
./inherit-align.ce 100 0 5 0.5 0 > tmp &
./inherit-align.ce 100 0 10 0.5 0 > tmp &
./inherit-align.ce 100 0 20 0.5 0 > tmp &
./inherit-align.ce 100 0 30 0.5 0 > tmp &
./inherit-align.ce 100 0 40 0.5 0 > tmp &
./inherit-align.ce 100 0 50 0.5 0 > tmp &
./inherit-align.ce 100 0 60 0.5 0 > tmp &
./inherit-align.ce 100 0 70 0.5 0 > tmp &
./inherit-align.ce 100 0 80 0.5 0 > tmp &
./inherit-align.ce 100 0 90 0.5 0 > tmp &
./inherit-align.ce 100 0 100 0.5 0 > tmp &
./inherit-align.ce 100 0 150 0.5 0 > tmp &
./inherit-align.ce 100 0 200 0.5 0 > tmp &

# default experiment
./inherit-align.ce 100 1 0 0.5 0 > tmp &
./inherit-align.ce 100 1 1 0.5 0 > tmp &
./inherit-align.ce 100 1 2 0.5 0 > tmp &
./inherit-align.ce 100 1 3 0.5 0 > tmp &
./inherit-align.ce 100 1 4 0.5 0 > tmp &
./inherit-align.ce 100 1 5 0.5 0 > tmp &
./inherit-align.ce 100 1 6 0.5 0 > tmp &
./inherit-align.ce 100 1 7 0.5 0 > tmp &
./inherit-align.ce 100 1 8 0.5 0 > tmp &
./inherit-align.ce 100 1 9 0.5 0 > tmp &
./inherit-align.ce 100 1 10 0.5 0 > tmp &
./inherit-align.ce 100 1 15 0.5 0 > tmp &
./inherit-align.ce 100 1 20 0.5 0 > tmp &

# smaller population size
./inherit-align.ce 50 0 0.5 0 > tmp &
./inherit-align.ce 50 5 0.5 0 > tmp &
./inherit-align.ce 50 10 0.5 0 > tmp &
./inherit-align.ce 50 20 0.5 0 > tmp &
./inherit-align.ce 50 30 0.5 0 > tmp &
./inherit-align.ce 50 40 0.5 0 > tmp &
./inherit-align.ce 50 50 0.5 0 > tmp &
./inherit-align.ce 50 60 0.5 0 > tmp &
./inherit-align.ce 50 70 0.5 0 > tmp &
./inherit-align.ce 50 80 0.5 0 > tmp &
./inherit-align.ce 50 90 0.5 0 > tmp &
./inherit-align.ce 50 100 0.5 0 > tmp &
./inherit-align.ce 50 150 0.5 0 > tmp &
./inherit-align.ce 50 200 0.5 0 > tmp &

# larger population size
./inherit-align.ce 200 0 0.5 0 > tmp &
./inherit-align.ce 200 5 0.5 0 > tmp &
./inherit-align.ce 200 10 0.5 0 > tmp &
./inherit-align.ce 200 20 0.5 0 > tmp &
./inherit-align.ce 200 30 0.5 0 > tmp &
./inherit-align.ce 200 40 0.5 0 > tmp &
./inherit-align.ce 200 50 0.5 0 > tmp &
./inherit-align.ce 200 60 0.5 0 > tmp &
./inherit-align.ce 200 70 0.5 0 > tmp &
./inherit-align.ce 200 80 0.5 0 > tmp &
./inherit-align.ce 200 90 0.5 0 > tmp &
./inherit-align.ce 200 100 0.5 0 > tmp &
./inherit-align.ce 200 150 0.5 0 > tmp &
./inherit-align.ce 200 200 0.5 0 > tmp &

# heteroplasmy penalty
./inherit-align.ce 100 0 0.5 0.05 > tmp &
./inherit-align.ce 100 5 0.5 0.05 > tmp &
./inherit-align.ce 100 10 0.5 0.05 > tmp &
./inherit-align.ce 100 20 0.5 0.05 > tmp &
./inherit-align.ce 100 30 0.5 0.05 > tmp &
./inherit-align.ce 100 40 0.5 0.05 > tmp &
./inherit-align.ce 100 50 0.5 0.05 > tmp &
./inherit-align.ce 100 60 0.5 0.05 > tmp &
./inherit-align.ce 100 70 0.5 0.05 > tmp &
./inherit-align.ce 100 80 0.5 0.05 > tmp &
./inherit-align.ce 100 90 0.5 0.05 > tmp &
./inherit-align.ce 100 100 0.5 0.05 > tmp &
./inherit-align.ce 100 150 0.5 0.05 > tmp &
./inherit-align.ce 100 200 0.5 0.05 > tmp &

# more comparable mitotype fitness
./inherit-align.ce 100 0 0.9 0 > tmp &
./inherit-align.ce 100 5 0.9 0 > tmp &
./inherit-align.ce 100 10 0.9 0 > tmp &
./inherit-align.ce 100 20 0.9 0 > tmp &
./inherit-align.ce 100 30 0.9 0 > tmp &
./inherit-align.ce 100 40 0.9 0 > tmp &
./inherit-align.ce 100 50 0.9 0 > tmp &
./inherit-align.ce 100 60 0.9 0 > tmp &
./inherit-align.ce 100 70 0.9 0 > tmp &
./inherit-align.ce 100 80 0.9 0 > tmp &
./inherit-align.ce 100 90 0.9 0 > tmp &
./inherit-align.ce 100 100 0.9 0 > tmp &
./inherit-align.ce 100 150 0.9 0 > tmp &
./inherit-align.ce 100 200 0.9 0 > tmp &

# more distinct mitotype fitness
./inherit-align.ce 100 0 0.1 0 > tmp &
./inherit-align.ce 100 5 0.1 0 > tmp &
./inherit-align.ce 100 10 0.1 0 > tmp &
./inherit-align.ce 100 20 0.1 0 > tmp &
./inherit-align.ce 100 30 0.1 0 > tmp &
./inherit-align.ce 100 40 0.1 0 > tmp &
./inherit-align.ce 100 50 0.1 0 > tmp &
./inherit-align.ce 100 60 0.1 0 > tmp &
./inherit-align.ce 100 70 0.1 0 > tmp &
./inherit-align.ce 100 80 0.1 0 > tmp &
./inherit-align.ce 100 90 0.1 0 > tmp &
./inherit-align.ce 100 100 0.1 0 > tmp &
./inherit-align.ce 100 150 0.1 0 > tmp &
./inherit-align.ce 100 200 0.1 0 > tmp &

