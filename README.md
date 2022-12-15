# trajectory_calc
Ballistic trajectory emulator

Do:

g++ trajectory_calc.cpp -o trajectory_calc

Then, if you want to see a single trajectory, run:

./trajectory_calc > trajectory.dat
python3 ./calcTrackAdv.py

If you want to create a table of simulations with 0.5 degrees launching
angle step run the program with some dummy parameter, like

./trajectory_calc qqq > table.txt



You may need:
pip3 install --user --upgrade matplotlib
pip3 install --user --upgrade numpy