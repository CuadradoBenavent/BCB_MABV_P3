This program computes the molecular dynamics of N-Argon atoms.

Molecular Dynamics (MD) simulates the movement of atoms based on their initial positions and velocities.
In this project, a MD program to illustrate key concepts in MD simulations is developed.
The program is structured in different functions. First,functions for allocating and freeing two-dimensional arrays are defined.
Then, the input file is read and the number of atoms and the coordinates are extracted.
The next step involves the calculation of the Lennard-Jones potential, the kinetic energy, the potential energy and the total energy.
Next, the acceleration is calculated solving Newton's equations.
Finally, the MD is implemented using the Verlet algorithm. A final movie can be visualized in software like VMD or Molden.

#####################################################################

Within this folder there are the following documents:

- AUTHORS file containing the names of all contributors.

- bin directory contains the compiled executable file.

- dynamics.pdf file containing the instructions followed to peform this homework.

- INSTALL.md file containing instructions on how to compile and run the program.

- LICENSE file specifying the license for the code, in this case a GNU GENERAL PUBLIC LICENSE. 

- MAKEFILE file used to compile and execute the molecular dynamics code, extract the trajectories and energies every 10 steps.

- README.md - this file containing information.

- src directory contains the main source code file in C necessary to develop the project.

- tests directory containing different inputs for different number of Argon atoms and their initial trajectories, also including the given example "inp.txt".
  Included also their trajectory and energy outputs.


#####################################################################
