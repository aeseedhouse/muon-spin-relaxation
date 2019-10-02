# muon-spin-relaxation
Modelling the magnetic field of a material felt by nearby spins

## Getting Started

For Windows the cygwin1.dll application needs to be installed to run C programs.


## File Explanations
The C files used to create the executables:
gauss_muon.c: models an x by y, 2D, "sheet" of spins at a distance z from a magnet
gauss_muon_3d.c:  models an x by y by z, 3D, map of spins near a magnet, z is the distance from the magnet
and gauss_random_3D.c :  models an x by y by z, 3D, map of spins near a magnet, z is the distance from the magnet

The config.csv can be edited to choose the input parameters for the desired simulation.

An output file will be created if the executables are ran, called "c_grid_mag_field_1.csv". This gives the x, y, and z locations of the spin particles (muons) and the Bx, By, and Bz calculated values of the magnetic field felt by that particular spin particle location.
