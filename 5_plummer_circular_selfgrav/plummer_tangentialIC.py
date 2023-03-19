################################################################################
# Copyright (c) 2022 Patrick Hirling (patrick.hirling@epfl.ch)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################


##################################################
# IC for a Plummer model with only circular orbits
##################################################

import numpy as np
import argparse
from swiftsimio import Writer
import unyt

parser = argparse.ArgumentParser(
    description="Generate discrete realization of Anisotropic Plummer Model"
)

# Parse Main Arguments
parser.add_argument("-a", type=float, default=0.05, help="Softening Length (in kpc)")
parser.add_argument("-M", type=float, default=1.0e-5, help="Total Mass of Model (in 10^10 solar masses)")
parser.add_argument("-N", type=int, default=1e5, help="Number of particles in the simulation")
parser.add_argument("-bound", type=float, default=2.0, help="Upper bound of radius sampling")
parser.add_argument("-boxsize", type=float, default=4.0, help="Box Size of simulation")
parser.add_argument("-o",type=str,default="plummer.hdf5",help="Output file name")

args = parser.parse_args()

#### Parameters
# Plummer Model
G = 4.299581e04  # Gravitational constant [kpc / 10^10 M_s * (kms)^2]
a = args.a  # Plummer softening length [kpc]
M = args.M  # Total Mass [10^10 M_s]

# IC File
N = int(args.N)  # Number of Particles
fname = args.o  # Name of the ic file


### Print parameters
print("\n################################################################################# \n")
print("Generating Plummer ICs with the following parameters::\n")
print("Softening length:                                  a = " + "{:.3e} kpc".format(a))
print("Total Mass:                                        M = "+ "{:.3e} Solar Masses".format(M * 1e10))
print("Anisotropy parameter:                              q = 0")
print("Output file:                                       " + fname)
print("Number of particles:                               N = " + "{:.1e}".format(N))
print("\n################################################################################# \n")

# Estimate good softening length (density at origin, sphere s.t. contains on avg 1 particle)
epsilon0 = a * N ** (-1.0 / 3.0)
print(
    "Recommended Softening length (times conventional factor): "
    + "{:.4e}".format(epsilon0)
    + " kpc"
)

### Generate Positions
# Inverse transform sampling from analytical inverse of M(<r)
def r_of_m(m, a, M):
    return a * ((M / m) ** (2.0 / 3.0) - 1) ** (-1.0 / 2.0)


m_rand = M * np.random.uniform(0.0, 1.0, N)
r_rand = r_of_m(m_rand, a, M)
phi_rand = np.random.uniform(0.0, 2 * np.pi, N)
theta_rand = np.arccos(np.random.uniform(-1.0, 1.0, N))

x = r_rand * np.sin(theta_rand) * np.cos(phi_rand)
y = r_rand * np.sin(theta_rand) * np.sin(phi_rand)
z = r_rand * np.cos(theta_rand)

X = np.array([x, y, z]).transpose()

### Generate Velocities for purely circular orbits
vv = np.sqrt(G * M * r_rand**2 * (r_rand**2 + a**2)**(-3. / 2.))

# Convert to Cartesian
# First: project vt on e_theta, e_phi with random orientation
alph = np.random.uniform(0, 2 * np.pi, N)
vphi = vv * np.cos(alph)
vtheta = vv * np.sin(alph)

# Convert Spherical to cartesian coordinates
v_x = (
    + np.cos(theta_rand) * np.cos(phi_rand) * vtheta
    - np.sin(phi_rand) * vphi
)
v_y = (
    + np.cos(theta_rand) * np.sin(phi_rand) * vtheta
    + np.cos(phi_rand) * vphi
)
v_z = - np.sin(theta_rand) * vtheta

# Create velocity array
V = np.array([v_x, v_y, v_z]).transpose()

### Generate masses

m = M / N * np.ones(N)

### Exclude extreme outliers
idx = r_rand < args.bound
X = X[idx]
V = V[idx]
m = m[idx]
new_N = len(m)

### Write to hdf5
galactic_units = unyt.UnitSystem(
    "galactic",
    unyt.kpc,
    unyt.unyt_quantity(1e10, units=unyt.Solar_Mass),
    unyt.unyt_quantity(1.0, units=unyt.s * unyt.Mpc / unyt.km).to(unyt.Gyr),
)
wrt = Writer(galactic_units, args.boxsize * unyt.kpc)
wrt.dark_matter.coordinates = X * unyt.kpc
wrt.dark_matter.velocities = V * (unyt.km / unyt.s)
wrt.dark_matter.masses = m * 1e10 * unyt.msun
wrt.dark_matter.particle_ids = np.arange(new_N)
wrt.write(fname)
print("Writing IC file...")