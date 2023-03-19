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

l = args.l # wavelength of the perturbation  if < lJ => stable
lJ = args.lJ # Jeans wavelength
n = int(args.N)     # number of particles
L = 1.        # box size

#N = 50        # number of points ????
m = 5         # mode
G = 1.
rho = 1.

if lJ != 0.0:
   kJ = 2*np.pi/lJ
   sigma = np.sqrt(4*np.pi*G*rho)/ kJ
else:
   sigma = 0.0

print("sigma = {}".format(sigma))

k = 2*np.pi/l # wave number

nb = ic.box(n, L/2.,L/2.,L/2., irand=1, name=args.o,ftype=args.f)
nb.mass = np.ones(n)/n * rho * L**3

nb.pos[:,0] = nb.pos[:,0] + L/2.
nb.pos[:,1] = nb.pos[:,1] + L/2.
#nb.pos[:,2] = nb.pos[:,2] * 0

#nb.pos[:,0] = nb.x() + e*np.cos(k*nb.x()) / k

nb.verbose = 1
if args.hydro:
   nb.set_tpe(0)
   nb.u = np.ones(n)*sigma**2
   hsml = nb.get_rsp_approximation()
   nb.rsp = hsml
else:
   nb.set_tpe(1)
   # nb.vel[:,0] = sigma* np.sqrt(3) *np.random.standard_normal([nb.nbody])
   # nb.vel[:,1] = sigma* np.sqrt(3) *np.random.standard_normal([nb.nbody])
   # nb.vel[:,2] = sigma* np.sqrt(3) *np.random.standard_normal([nb.nbody])
np.random.normal(scale=args.sigma,size=N)
   nb.vel[:,0] = np.random.normal(scale=sigma/np.sqrt(3),size=N)
   nb.vel[:,1] = np.random.normal(scale=sigma/np.sqrt(3),size=N)
   nb.vel[:,2] = np.random.normal(scale=sigma/np.sqrt(3),size=N)
   print("vmax = {}".format(nb.vx().max()))

nb.boxsize = L
nb.rebox()
nb.write()
