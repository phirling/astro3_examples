# Astrophysics III (EPFL) - Course Examples
This repo contains the scripts and instructions to reproduce the N-Body simulations used as examples during the lectures.
In the following, "collisionless" or "dark matter" means that the particles interact only gravitationally, as opposed to
"hydrodynamic", when particles interact with gas-dynamic effects as well (only simulation 3).
The current examples are:
1. Razor-thin (2D), infinite (periodic) homogeneous slab of particles interacting through their self-gravity.
2. homogeneous 3D box of particles interacting through their self-gravity.
3. homogeneous 3D box of particles interacting through their self-gravity AND hydrodynamics.
4. Plummer model with only circular orbits, using *external gravity*.
5. Plummer model with only circular orbits, using *internal* or *self-gravity*.
6. plummer model with intermediate anisotropy in the orbits, using self-gravity.
7. Miyamoto-Nagai disk with only circular orbits and external gravity.
8. Miyamoto-Nagai disk with only circular orbits and self-gravity.
9. Miyamoto-Nagai disk with circular orbits rotating in both directions (50:50) and external-gravity.
10. Miyamoto-Nagai disk with circular orbits rotating in both directions (50:50) and self-gravity.
11. Thin slab of particles of nonzero thickness, interacting through self-gravity.

## Usage

### Requirements
The scripts to generate initial conditions rely on
* `Numpy`
* `Scipy`
* `pNbody`: https://gitlab.com/revaz/pNbody

To run the simulations, we use the [SWIFT](https://swift.dur.ac.uk/) N-Body code, whose build instructions can be found [here](https://swift.dur.ac.uk/docs/GettingStarted/compiling_code.html).

### Running the examples
Each example requires SWIFT to be compiled with slightly different configurations, which can be found in the Readme files of the
respective example subdirectories. The general compilation process looks like this (in the SWIFT source directory):
```
# Generate Config
./autogen.sh
# Configure SWIFT
./configure --disable-compiler-warnings --disable-doxygen-doc --disable-hand-vec --with-hdf5=[path/to/hdf5/lib] --with-metis=[path/to/METIS/lib] [Additional, example dependent flags]
# Compile using 4 threads
make -j 4
```
### Output Visualization
TODO