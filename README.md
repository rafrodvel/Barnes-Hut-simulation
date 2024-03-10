# Barnes-Hut-simulation

This c program implements a quad-tree to solve Barnes Hut simulation with Verlet Velocity integration method and uses local parallelization with openmp. It reads input data from input_data folder and writes out result.gal with the positions and velocities at the end of the simulation.

Usage: ./galsim N_bodies input_file Time_steps Delta_t theta