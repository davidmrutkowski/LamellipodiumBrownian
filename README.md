Brownian dynamics code that simulates a lamellipodium at the discrete filament level undergoing retrograde flow in the presence of a nascent focal adhesion. This code was developed for the simulations in Rutkowski and Vavylonis, PLOS Comp. Biol. 17: e1009506 (2021), https://doi.org/10.1371/journal.pcbi.1009506

Compilation using g++:
g++ -O2 -fopenmp -o main main.cpp Filament.cpp Grid.cpp Particle.cpp ParticleInfo.cpp Simulation.cpp SimulationBox.cpp SpringBondInfo.cpp TaggedVector.cpp

Input parameter file (input.txt):
This file lists various parameters used by the simulation. 
In particular the numThreads parameter of 10 should be changed to the number of threads the code is desired to be run on and should not exceed the number of threads available on a given computer.
The parameters in the given input.txt are for simulating a kappa_FA value of 100 (etaInNascentAdhesion / eta) in the pull uniform mode (uniformPullingForceOn true and backPullinForceOn false).

Command to run code:
./main input.txt

Restart from previously run simulation:
Run GetOutputFrameFromXYZ.py specifying the xyz output file and the last frame number
Run GetOutputBondsFromBND.py specifying the bnd output file and the last frame number
Add the following lines to input.txt
positions kFA-100_1-frame###.xyz
bonds kFA-100_1-frame###.bnd
Change the outputName parameter to a new name to keep from overwriting previous output.
