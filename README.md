# trajectory_calc
Ballistic trajectory emulator

Do:

g++ trajectory_calc.cpp -o trajectory_calc

./trajectory_calc > trajectory.dat

Then do:

gnuplot

> splot './trajectory.dat' using 1:2:3 with lines
