import numpy as np
import argparse
from swiftsimio import Writer
import unyt

parser = argparse.ArgumentParser(
    description=""
)

parser.add_argument("-N", type=int, default=250000, help="Number of particles in the simulation")
parser.add_argument("-sigmaxy",type=float,default=1)
parser.add_argument("-sigmaz",type=float,default=1)
parser.add_argument("-rho0",type=float,default=1)
parser.add_argument("-o",type=str,default="slab_thick.hdf5",help="Output Filename")

args = parser.parse_args()

a = 1.0          # Box Size
N = int(args.N)  # Number of particles
G = 1.0          # Gravitational constant
M = 1.0          # Total Mass
rho0 = 1.0
max_z = 0.5

# Parameters
z0 = M/(4.*rho0)
sigma_z = M*np.sqrt(np.pi*G/(2*rho0)) #M*np.sqrt(np.pi*G/2.)

# Sample x & y
x = np.random.uniform(0.,a,N)
y = np.random.uniform(0.,a,N)

# Sample z
def z_of_M(M):
    return z0 * np.log((2*args.rho0*z0 + M)/(2*args.rho0*z0 - M))

print(z0)
Mgen = np.random.uniform(0.,2*rho0*z0,N)
z = np.random.choice([-1,1],N)*z_of_M(Mgen)

X = np.array([x, y, z]).transpose()

# Sample velocities
v_x = np.random.normal(0.,args.sigmaxy,N)
v_y = np.random.normal(0.,args.sigmaxy,N)
v_z = np.random.normal(0.,args.sigmaz,N)

V = np.array([v_x, v_y, v_z]).transpose()

m = M / N * np.ones(N)

# Estimate good softening length radius of circle that contains on avg one particle
epsilon0 = a / np.sqrt(N*np.pi)
print(
    "Recommended Softening length (times conventional factor): "
    + "{:.4e}".format(epsilon0)
    + " kpc"
)

# Convention: exclude particles with z > +/- 0.5
idx = np.abs(X[:,2]) < max_z
X = X[idx]
V = V[idx]
print("Excluded {:n} outliers".format(N - len(idx)))

### Write to hdf5
galactic_units = unyt.UnitSystem(
    "G_is_one",
    unyt.kpc,
    unyt.unyt_quantity(1e10, units=unyt.Solar_Mass),
    unyt.unyt_quantity(1.0, units=unyt.s * unyt.Mpc / unyt.km).to(unyt.Gyr),
)
boxsize = unyt.kpc * np.array([1,1,1]) # args.a * unyt.kpc
wrt = Writer(galactic_units, boxsize)
wrt.dark_matter.coordinates = X * unyt.kpc
wrt.dark_matter.velocities = V * (unyt.km / unyt.s)
wrt.dark_matter.masses = m * 1e10 * unyt.msun
wrt.dark_matter.particle_ids = np.arange(N)
wrt.write(args.o)
print("Writing IC file...")

# For Debugging
import matplotlib.pyplot as plt
plt.hist(X[:,2],100)
plt.show()