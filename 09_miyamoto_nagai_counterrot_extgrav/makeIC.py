import numpy as np
from pNbody import ic
from pNbody.mass_models import miyamoto
import argparse

"""
Script to generate IC of a Miyamoto-Nagai Disk with only circular orbits,
such that half of the particles are rotating in one direction and the other
half in the other. Can be used either in self-gravity mode (unstable) or
external gravity mode.

This script assumes a specific choice of internal units (see below), that
are chosen so that G=1.
"""

parser = argparse.ArgumentParser(description="Create IC of a Miyamoto-Nagai Disk with counter-rotation")
parser.add_argument("-a",type=float,default=1.0,help="Disk scale length of M-N Potential")
parser.add_argument("-b",type=float,default=1.0,help="Disk scale height of M-N Potential")
parser.add_argument("-Rmax",type=float,default=10.0,help="Maximal radius to sample")
parser.add_argument("-Zmax",type=float,default=10.0,help="Maximal height (z) to sample")
parser.add_argument("-N",type=int,default=25000,help="Number of Particles")
parser.add_argument("-o",type=str,default="miyamoto.hdf5",help="Output Filename")
parser.add_argument("-f",type=str,default='swift',help="Format of the IC file (swift,gadget")
args = parser.parse_args()

# Physical Parameters for the Miyamoto-Nagai model
G = 1.
M = 1.
n = int(args.N)
a = args.a
b = args.b

# Create NBody object with MN distribution, set particle type and box size
nb = ic.miyamoto_nagai(n,a,b,args.Rmax,args.Zmax,ftype=args.f,name=args.o)
nb.set_tpe(1)
nb.boxsize = 2*args.Rmax

# Fix internal units s.t. G=1
nb.UnitLength_in_cm = 3.086e+21
nb.UnitVelocity_in_cm_per_s = 9.78469e+07
nb.UnitMass_in_g = 4.4356e+44

# Get cylindrical coordinates
r = nb.rxy()
z = nb.pos[:,2]
phi = nb.phi_xy()

# Compute Circular Velocity
vc = miyamoto.Vcirc(M,a,b,r)

# Project on cartesian coordinates
v_x = -np.sin(phi) * vc
v_y = np.cos(phi) * vc
v_z = np.zeros(n)

# Invert rotation direction for half of the particles
mid = int(n / 2)
v_x[mid:] *= -1
v_y[mid:] *= -1

# Set Velocities
nb.vel[:,0] = v_x
nb.vel[:,1] = v_y
nb.vel[:,2] = v_z

# Set Mass
nb.mass = np.ones(n) * M/n

# Write to IC file
nb.write()


# ===========
# Unused Code
# ===========

# Circular Velocities
# F1 = np.sqrt(G*M)
# def Vc(r,z):
#     return F1 * r / (r**2 + (np.sqrt(z**2 + b**2) + a)**2)**(3./4)
# 
# def Vc2(r2,z):
#     return G*M*r2 / (r2 + (np.sqrt(z**2 + b**2) + a)**2)**(3./2)