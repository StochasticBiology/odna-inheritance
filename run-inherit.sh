gcc -o3 inherit.c -lm -o inherit.ce

./inherit.ce 0 0 0 > tmp &
./inherit.ce 1 0 0 > tmp &
./inherit.ce 2 0 0 > tmp &
./inherit.ce 3 0 0 > tmp &
./inherit.ce 4 0 0 > tmp &
./inherit.ce 5 0 0 > tmp &
./inherit.ce 6 0 0 > tmp &

./inherit.ce 0 0.5 0 > tmp &
./inherit.ce 1 0.5 0 > tmp &
./inherit.ce 2 0.5 0 > tmp &
./inherit.ce 3 0.5 0 > tmp &
./inherit.ce 4 0.5 0 > tmp &
./inherit.ce 5 0.5 0 > tmp &
./inherit.ce 6 0.5 0 > tmp &

./inherit.ce 0 0 1 > tmp &
./inherit.ce 1 0 1 > tmp &
./inherit.ce 2 0 1 > tmp &
./inherit.ce 3 0 1 > tmp &
./inherit.ce 4 0 1 > tmp &
./inherit.ce 5 0 1 > tmp &
./inherit.ce 6 0 1 > tmp &

./inherit.ce 0 0.5 1 > tmp &
./inherit.ce 1 0.5 1 > tmp &
./inherit.ce 2 0.5 1 > tmp &
./inherit.ce 3 0.5 1 > tmp &
./inherit.ce 4 0.5 1 > tmp &
./inherit.ce 5 0.5 1 > tmp &
./inherit.ce 6 0.5 1 > tmp &

gcc -o3 baseline.c -lm -o baseline.ce

./baseline.ce 0 0 0 > tmp &
./baseline.ce 1 0 0 > tmp &
./baseline.ce 0 0.5 0 > tmp &
./baseline.ce 1 0.5 0 > tmp &
./baseline.ce 0 0 1 > tmp &
./baseline.ce 1 0 1 > tmp &
./baseline.ce 0 0.5 1 > tmp &
./baseline.ce 1 0.5 1 > tmp &

gcc -o3 inherit-zoom.c -lm -o inherit-zoom.ce
./inherit-zoom.ce 3 0.5 0 > tmp &
./inherit-zoom.ce 4 0.5 0 > tmp &
./inherit-zoom.ce 5 0.5 0 > tmp &
./inherit-zoom.ce 7 0.5 0 > tmp &
./inherit-zoom.ce 8 0.5 0 > tmp &
./inherit-zoom.ce 9 0.5 0 > tmp &
./inherit-zoom.ce 10 0.5 0 > tmp &
