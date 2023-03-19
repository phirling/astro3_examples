import numpy as np
from pNbody import Nbody
import argparse


"""
================== EXAMPLE 11: Slab of Finite, Nonzero Thickness ==================

The model is described in Ex. 4.21 of Binney & Tremaine.
In x and y, the particles are uniformly distributed and follow a gaussian
isotropic velocity distribution, while in the z direction they follow a
sech^2 profile.

This model is supposed to run with periodic boundary conditions (for x and y).
Since SWIFT can only run periodic BC in either all or no dimensions, this means
that in z, the model has to be thin enough for the contribution from periodicity
to be negligible. 
===================================================================================
"""

parser = argparse.ArgumentParser(description="Create IC of an infinite slab of finite nonzero thickness")
parser.add_argument("-z0",type=float,default=0.03,help="Scaling height of slab")
parser.add_argument("-xydisp",type=float,default=1.5,help="Velocity Dispersion (sigma^2) in the xy direction (isotropic)")
parser.add_argument("-zmax",type=float,default=0.3,help="Maximum allowed actual z (in model, this is potentially infinite)")
parser.add_argument("-N",type=int,default=100000,help="Number of Particles")
parser.add_argument("-o",type=str,default="slab.hdf5",help="Output Filename")
parser.add_argument("-f",type=str,default='swift',help="Format of the IC file (swift,gadget")

args = parser.parse_args()


# Model Parameters
L = 1.                  # Box size
zmax = args.zmax        # Maximum allowed value of z         
n = int(args.N)         # Number of Particles
G = 1.                  # Gravitational Constant
M = 1.                  # Total Mass
z0 = args.z0            # Scaling height
L2 = L*L                # Area
rho0 = M/(4*L2*z0)      # Scaling Density

"""
=== Positions ===
Since the problem is symmetric in z, we consider only the upper
plane for sampling and then add a random sign.
"""

# Integrated mass on one side of the plane
def M_of_z(z): 
    return 2*L2*rho0*z0*np.tanh(0.5*z/z0)

# Inverse function
def z_of_M(M):
    return 2*z0*np.arctanh(M/(2*L2*rho0*z0))

# Maximum (half-)mass to enforce the zmax parameter
Mmax = M_of_z(zmax)

# Inverse transform sampling for z, with random sign
Mgen = np.random.uniform(0.,Mmax,n)
z = np.random.choice([-1,1],n) * z_of_M(Mgen)

# x and y are just uniformly distributed on the plane
x = np.random.uniform(0.,L,n)
y = np.random.uniform(0.,L,n)

# Position array
pos = np.vstack((x,y,z)).T

"""
=== Velocities ===
In the xy directions, the velocities are isotropic and follow a
gaussian distribution. In the z direction, also gaussian but with a
dispersion that is determined by other parameters through the DF
"""

sigz = z0*np.sqrt(8*np.pi*G*rho0)
sigxy = np.sqrt(args.xydisp/2)

v_x = np.random.normal(0.0,sigxy,n)
v_y = np.random.normal(0.0,sigxy,n)
v_z = np.random.normal(0.0,sigz,n)
vel = np.vstack((v_x,v_y,v_z)).T

# Create Nbody object, assign G=1 units and write to hdf5
nb = Nbody(pos=pos,vel=vel,ftype=args.f,p_name=args.o,verbose=1,status='new')
nb.set_tpe(1)
nb.UnitLength_in_cm = 3.086e+21
nb.UnitVelocity_in_cm_per_s = 9.78469e+07
nb.UnitMass_in_g = 4.4356e+44
nb.boxsize = 1.0
nb.write()

# ====================================
# Testing, use only with few particles
# ====================================

if 0:
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(pos[:,0],pos[:,1],pos[:,2],s=1)
    ax.set_aspect('equal')

    plt.show()