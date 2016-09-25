# CS420
This repository contains coding projects written for the CS420: Biologically-inspired Computation course offered at the University of Tennessee, Knoxville. Data spreadsheets, code files, and report documents for projects 2-8 are included. Each project involved creating simulation software for a particular topic covered in the course and collecting meaningful data from these simulation experiments. Conclusions were then drawn about how particular simulation parameters contributed to the overall behavior of the system being modeled. 

### Project 2
For Project 2, I wrote a simulation program to model various 2D activator/inhibitor cellular automata (AICA). Measurements were taken to investigate the formation of spatial structures within these systems. Several parameters were considered for predicting which AICA would produce more "interesting" structural patterns.

### Project 3
Project 3 involved the creation of a basic Hopfield neural network. Experiments were conducted to evaluate the associative memory capacity of a general Hopfield net.

### Project 4
In Project 4, I designed a Back-Propagation Software System that allowed for the specification of the following parameters: number of layers, number of neurons in each layer, learning rate, number of training epochs, and the data files to use for training, validation, and testing. The parameters were adjusted to find the optimal values for solving two algebraic functions, one that took 2 input values and one that took 3 input values.

### Project 5
For Project 5, I created a system that implemented a basic genetic algorithm. It started by generating a population of individual "organisms", each of which was defined by a particular genetic string (comprised of some number of binary bits). New generations were spawned by influencing each subsequent population with some probability for individual genetic mutation or crossover between two genetic strings. A simple fitness function was incorporated such that more "fit" individuals had a greater chance of passing on their genetic string to offspring in the next generation.

### Project 6
Project 6 dealt with the implementation of the Particle Swarm Optimization (PSO) algorithm. This algorithm models a system in which candidate solutions("particles") to a problem move around a search-space according to position and velocity functions. Each particle is influenced by its own local best known position in the space, but it is also influenced by the current best known position found by any other particle in the space. Over time, the majority of particles will usually migrate toward an optimal solution. The simulation software I wrote took the following arguments: number of epochs, number of particles, inertia, cognitive influence weight, and social influence weight. These values were manipulated across many experiments, and the resulting changes in swarm behavior were observed.  
