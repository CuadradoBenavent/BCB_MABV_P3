To execute this program:

1) Ensure that the desired input and output files are specified in the MAKEFILE file. Modify section:

    input_file="tests/inp.txt"
    output_trajectory_file="tests/trajectory.xyz"
    output_energy_file="tests/energy.dat"


2) Run:    ./MAKEFILE

This file compiles the C code, executes the generated executable. 
The molecular dynamics is performed including the writing of a trajectory file with geometries and energies every 10 steps. 
The energies are extracted onto a .dat file to be used for plotting.
