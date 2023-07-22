import numpy as np
from pNbody import ic
import argparse

parser = argparse.ArgumentParser(description="Create IC of a box with optional density perturbation")

parser.add_argument("-l",type=float,default=0.25,help="wavelength of the perturbation")
parser.add_argument("-lJ",type=float,default=0.25,help="Jeans Wavelength")
parser.add_argument("-o",type=str,default="box.hdf5",help="Output Filename")
parser.add_argument("--hydro",action='store_true',help="Use Hydro Mode")
parser.add_argument("-N",type=int,default=250000,help="Number of Particles")
parser.add_argument("-f",type=str,default='swift',help="Format of the IC file (swift,gadget")

args = parser.parse_args()

# Box Parameters
l = args.l           # Wavelength of the perturbation  if < lJ => stable
lJ = args.lJ         # Jeans wavelength
n = int(args.N)      # Number of particles
L = 1.               # Simulation Box Size
m = 5                # Perturbation Mode
G = 1.               # Gravitational Constant
rho = 1.             # Mean mass density of box

# Find dispersion from Jeans length
if lJ != 0.0:
   kJ = 2*np.pi/lJ
   sigma = np.sqrt(4*np.pi*G*rho)/ kJ
else:
   sigma = 0.0

print(f"Sigma = {sigma:.3f}")

k = 2*np.pi/l # wave number


# Create Nbody object from IC module and set unit system
nb = ic.box(n, L/2.,L/2.,L/2., irand=1, name=args.o,ftype=args.f)
nb.set_tpe(1)
nb.verbose = 1
nb.UnitLength_in_cm = 3.086e+21
nb.UnitVelocity_in_cm_per_s = 9.78469e+07
nb.UnitMass_in_g = 4.4356e+44
nb.boxsize = L

# Shift Positions to fit in box
nb.pos += L/2.0

# Perturbations
#nb.pos[:,0] = nb.x() + e*np.cos(k*nb.x()) / k

# Set Velocities
nb.vel[:,0] = np.random.normal(scale=sigma/np.sqrt(3),size=n)
nb.vel[:,1] = np.random.normal(scale=sigma/np.sqrt(3),size=n)
nb.vel[:,2] = np.random.normal(scale=sigma/np.sqrt(3),size=n)
print("vmax = {}".format(nb.vx().max()))

# Set Mass
nb.mass = np.ones(n)/n * rho * L**3

nb.write()