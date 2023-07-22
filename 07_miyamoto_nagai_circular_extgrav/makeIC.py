import numpy as np
from pNbody import ic
from pNbody.mass_models import miyamoto
import argparse

"""
================== EXAMPLE 7: Miyamoto-Nagai Disk, circular orbits ==================

[Add Description]

=====================================================================================
"""

parser = argparse.ArgumentParser(description="Create IC of a Miyamoto-Nagai Disk")
parser.add_argument("-a",type=float,default=1.0,help="A parameter of M-N Potential")
parser.add_argument("-b",type=float,default=1.0,help="B parameter of M-N Potential")
parser.add_argument("-Rmax",type=float,default=10.0,help="")
parser.add_argument("-Zmax",type=float,default=10.0,help="")
parser.add_argument("-N",type=int,default=25000,help="Number of Particles")
parser.add_argument("-o",type=str,default="miyamoto.hdf5",help="Output Filename")
parser.add_argument("-f",type=str,default='swift',help="Format of the IC file (swift,gadget")
args = parser.parse_args()


# Miyamoto-Nagai Parameters
G = 1.                  # Gravitational Constant
M = 1.                  # Total Mass
n = int(args.N)         # Number of particles
a = args.a              # Scaling 1
b = args.b              # Scaling 2
boxsize = 20.0          # Simulation Box size

# Create Nbody object from IC module and set unit system
nb = ic.miyamoto_nagai(n,a,b,args.Rmax,args.Zmax,ftype=args.f,name=args.o)
nb.set_tpe(1)
nb.UnitLength_in_cm = 3.086e+21
nb.UnitVelocity_in_cm_per_s = 9.78469e+07
nb.UnitMass_in_g = 4.4356e+44
nb.boxsize = boxsize

# Convert to Cylindrical coords
r2 = nb.pos[:,0]**2 + nb.pos[:,1]**2
r = np.sqrt(nb.pos[:,0]**2 + nb.pos[:,1]**2)
r = nb.rxy()
z = nb.pos[:,2]
phi = nb.phi_xy()

# Circular Velocities
F1 = np.sqrt(G*M)
def Vc(r,z):
    return F1 * r / (r**2 + (np.sqrt(z**2 + b**2) + a)**2)**(3./4)

def Vc2(r2,z):
    return G*M*r2 / (r2 + (np.sqrt(z**2 + b**2) + a)**2)**(3./2)

vc = miyamoto.Vcirc(M,a,b,r)
#vc = np.sqrt(Vc2(r2,z))
v_x = -np.sin(phi) * vc
v_y = np.cos(phi) * vc
v_z = np.zeros(n)

# Set Velocities
nb.vel[:,0] = v_x
nb.vel[:,1] = v_y
nb.vel[:,2] = v_z

# Set Mass
nb.mass = np.ones(n) * M/n

nb.write()

