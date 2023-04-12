##################################################
# IC for a Plummer model with only circular orbits
##################################################

import numpy as np
import argparse
from pNbody import Nbody

parser = argparse.ArgumentParser()

# Parse Main Arguments
parser.add_argument("-a", type=float, default=0.1, help="Softening Length)")
parser.add_argument("-M", type=float, default=1.0, help="Total Mass of Model")
parser.add_argument("-N", type=int, default=250000, help="Number of particles in the simulation")
parser.add_argument("-bound", type=float, default=0.5, help="Upper bound of radius sampling")
parser.add_argument("-o",type=str,default="plummer.hdf5",help="Output file name")
parser.add_argument("-f",type=str,default='swift',help="Format of the IC file (swift,gadget")

args = parser.parse_args()

# Model Parameters
G = 1               # Gravitational constant
a = args.a          # Plummer softening length [kpc]
M = args.M          # Total Mass [10^10 M_s]
N = int(args.N)     # Number of Particles

# Generate positions: Inverse transform sampling from analytical inverse of M(<r)
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

# Generate Velocities for purely circular orbits
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

V = np.array([v_x, v_y, v_z]).transpose()

# Generate masses
m = M / N * np.ones(N)

# Exclude extreme outliers
idx = r_rand < args.bound
X = X[idx]
V = V[idx]
m = m[idx]
new_N = len(m)

# Create Nbody object, assign G=1 units and write to hdf5
nb = Nbody(pos=X,vel=V,ftype=args.f,p_name=args.o,verbose=1,status='new')
nb.set_tpe(1)
nb.UnitLength_in_cm = 3.086e+21
nb.UnitVelocity_in_cm_per_s = 9.78469e+07
nb.UnitMass_in_g = 4.4356e+44
nb.boxsize = 1.0
nb.write()