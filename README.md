# Barnes-Hut-simulation

This c program implements a quad-tree to solve Barnes Hut simulation with Verlet Velocity integration method and uses local parallelization with openmp. It reads input data from input_data folder and writes out result.gal with the positions and velocities at the end of the simulation.

Usage: ./galsim N_bodies filename nsteps delta_t theta_0 processes

• N_bodies: The number of bodies in the simulation.
• filename: The name of the file containing the initial positions, velocities, masses, and
brightness of the bodies.
• nsteps: The number of time steps to run the simulation.
• delta_t: The time step of the simulation.
• theta_0: The parameter θ of the Barnes-Hut algorithm.
• processes: The number of threads to use in the parallelization.

Read the report for more information
