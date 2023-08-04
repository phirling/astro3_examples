# Astrophysics III (EPFL) - Course Examples
This repo contains the scripts and instructions to reproduce the N-Body simulations used as examples during the lecture
"Astrophysics III: Stellar and Galactic Dynamics" at EPFL. In the following, "collisionless" or "dark matter" means that
the particles interact only gravitationally, as opposed to "hydrodynamic", when particles interact with gas-dynamic
effects as well (only example 3).
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
The scripts to generate initial conditions and to visualize results depend on
* `Numpy`
* `Scipy`
* `h5py`
* `tqdm`
* `pNbody` (https://gitlab.com/revaz/pNbody)

To run the simulations, we use the [SWIFT](https://swift.dur.ac.uk/) N-Body code, whose build instructions and dependencies can be found [here](https://swift.dur.ac.uk/docs/GettingStarted/compiling_code.html).

### Build Instructions
To build SWIFT for the examples provided here, [clone](https://gitlab.cosma.dur.ac.uk/swift/swiftsim.git) the repo, enter the `swiftsim` directory and run:
```
./autogen.sh

./configure --with-ext-potential=[External Potential] --with-hdf5=[path/to/hdf5] --with-metis=[path/to/METIS]

# Compile using 4 threads (more if you wish)
make -j 4
```
Here, *[External Potential]* is be the name of an external potential used by SWIFT, for the examples that require it. It will be called *--with-ext-potential=point-mass-softened* for example 4 and *--with-ext-potential=nfw-mn* for examples 7-10. For the other examples, it can be either one since the external potential will be disabled at runtime. The location of the *hdf5* and *metis* libraries are
system-dependent, but on a personal linux install are typically in /usr/. On macOS, they are typically in /usr/local/. If the compilation fails, especially on macOS, you can try appending the following flags to the configure call: `--disable-compiler-warnings --disable-doxygen-doc --disable-hand-vec`

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
The provided python script `visualize_swift.py` can be used to easily produce images and films of the simulation results. Using the default config, SWIFT will output a set of
`snapshot_*.hdf5` files that can be read in by the script. When a single file name is passed, the script will produce an image, when multiple files are provided, it will produce an animation.
The details of the script arguments can be viewed by typing `python visualize_swift.py --help`. As an example, to read the output of example 1 and encode it into a mp4 file, one would run
```
python visualize_swift.py 1_razor_thin_slab_collisionless/snap/snapshot_*.hdf5 -shift 0.5 -nbins 400 -interp kaiser --nounits --savefilm
```
Here,
* The `-shift` flag shifts the positions of the particles in all 3 dimensions, since SWIFT works only on positive $x,y,z$ values while the visualization script is centered on the
origin.
* `-nbins` sets the number of bins used in each dimension, and so effectively the resolution of the image. The more particles used, the higher the resolution can be set while
preserving the relevant information, but for the examples here (with $1\cdot 10^5 - 3\cdot 10^5$ particles), using 400 bins (the default value) is usually fine.
* `-interp` sets the interpolation used (none by default) to produce the image. *kaiser* produces a nice detailed image, but e.g. *gaussian* or *bicubic* can be used to
achieve smoother results.
* The `--nounits` flag disables units in the axes and time indicator, since here we are dealing with non-physical units.
* The `--savefilm` flag is used to encode the output to a mp4 file. Omitting it will result in the animation being shown directly, although probably with a low framerate. This can be useful to
quickly visualize results.


### Notes
All examples use a unit system in which $G=1$, and a unity box size. Other parameters relevant to the model are scaled accordingly, but keep in
mind that this means that both the units and the scales displayed here are unphysical in the sense that they're not connected in any way to the
actual scales of the physical systems that one might think of (globular clusters, galactic disks, etc).
