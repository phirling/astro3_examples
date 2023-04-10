# Astrophysics III (EPFL) - Course Examples
This repo contains the scripts and instructions to reproduce the N-Body simulations used as examples during the lectures.
In the following, "collisionless" or "dark matter" means that the particles interact only gravitationally, as opposed to
"hydrodynamic", when particles interact with gas-dynamic effects as well (only example #3).
The current examples are:
1. Razor-thin (2D), infinite (periodic) homogeneous slab of particles interacting through their self-gravity.
2. Homogeneous 3D box of particles interacting through their self-gravity.
3. Homogeneous 3D box of particles interacting through their self-gravity AND hydrodynamics.
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
* `pNbody` (https://gitlab.com/revaz/pNbody)

To run the simulations, we use the [SWIFT](https://swift.dur.ac.uk/) N-Body code, whose build instructions can be found [here](https://swift.dur.ac.uk/docs/GettingStarted/compiling_code.html).

### Build Instructions
To build SWIFT for the examples provided here, [clone](https://gitlab.cosma.dur.ac.uk/swift/swiftsim.git) the repo, enter the `swiftsim` directory and run:
```
./autogen.sh

./configure --disable-compiler-warnings --disable-doxygen-doc --disable-hand-vec --with-ext-potential=[External Potential] --with-hdf5=[path/to/hdf5/lib] --with-metis=[path/to/METIS/lib]

# Compile using 4 threads (more if you wish)
make -j 4
```
Here, *[External Potential]* is be the name of an external potential used by SWIFT, for the examples that require it. It will be called *--with-ext-potential=point-mass-softened* for Ex. 4 and *--with-ext-potential=nfw-mn* for Ex. 7-10. For the other examples, it can be either one since the external potential will be disabled at runtime. The location of the *hdf5* and *metis* libraries are system-dependent.

### Running the examples
The examples contain run scripts that assume the `swiftsim` directory (containing the compiled SWIFT executable) to be at the root of the repo. If you installed SWIFT elsewhere, please change the location at the start of the scripts.
Furthermore, each example directory contains:
* A SWIFT parameter file `params.yml` with appropriate options for the example in question.
* A python script `makeIC.py`to generate the initial conditions.
* A README file with example-specific usage instructions

To try an example, assuming SWIFT is installed at the default location, simply run:
```
mkdir snap                      # Default directory where SWIFT stores its output files
python makeIC.py [options]      # Generate initial conditions
./run.sh                        # Run SWIFT
```

### Output Visualization
TODO