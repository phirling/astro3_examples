###############################################################################################################################
# IC for a Plummer model with constant nonzero beta
# See https://obswww.unige.ch/~revaz/pNbody/rst/GeneratingVelocities.html#example-set-a-plummer-sphere-to-the-jeans-equilibrium
###############################################################################################################################

import numpy as np
import argparse
from pNbody import Nbody,ic

parser = argparse.ArgumentParser()

# Parse Main Arguments
parser.add_argument("-a", type=float, default=1.0, help="Softening Length)")
parser.add_argument("-N", type=int, default=250000, help="Number of particles in the simulation")
parser.add_argument("-beta", type=float, default=0.0, help="Binney Anisotropy Parameter")
parser.add_argument("-o",type=str,default="plummer.hdf5",help="Output file name")
parser.add_argument("-f",type=str,default='swift',help="Format of the IC file (swift,gadget")

args = parser.parse_args()

# Model Parameters
M = 100.0
G = 1.0
a = args.a          # Plummer softening length [kpc]
N = int(args.N)     # Number of Particles
rmax = 10.0

# Estimate Dynamical time:
v_a = np.sqrt(G * M * a**2 * (a**2 + a**2)**(-3. / 2.))
tdyn = 2*np.pi*a / v_a
print("Dynamical time estimate (internal units): ",tdyn)

def r_of_m(m, a, M):
    return a * ((M / m) ** (2.0 / 3.0) - 1) ** (-1.0 / 2.0)
def m_of_r(r,a,M):
    return M*(r**3)/(r**2 + a**2)**(3./2)

Mmax = m_of_r(rmax,a,M)
m_rand = np.random.uniform(0.,Mmax,N)
r_rand = r_of_m(m_rand,a,M)

phi_rand = np.random.uniform(0.0, 2 * np.pi, N)
theta_rand = np.arccos(np.random.uniform(-1.0, 1.0, N))

x = r_rand * np.sin(theta_rand) * np.cos(phi_rand)
y = r_rand * np.sin(theta_rand) * np.sin(phi_rand)
z = r_rand * np.cos(theta_rand)

X = np.array([x, y, z]).transpose()


V = np.zeros_like(X)

mass = M/N * np.ones(N)

nb = Nbody(status='new',
        p_name=args.o,
        pos=X,
        vel=V,
        mass=mass,
        ftype=args.f)

# Generate mass distribution
# nb = ic.plummer(N,1,1,1,eps=a,rmax=rmax,name=args.o,ftype=args.f)

# Generate velocities
grmin = 0         # grid minimal radius
grmax = rmax*1.05 # grid maximal radius
nr = 128          # number of radial bins
rc = a
g  = lambda r:np.log(r/rc+1)       # the function
gm = lambda r:rc*(np.exp(r)-1)     # and its inverse

eps = 0.01 # gravitational softening
nb,phi,stats = nb.Get_Velocities_From_Spherical_Grid(eps=eps,nr=nr,rmax=grmax,phi=None,g=g,gm=gm,UseTree=True,ErrTolTheta=0.5,beta=0)

nb.vel *= 1
print(nb.vel)

nb.set_tpe(1)
#nb.UnitLength_in_cm = 3.085e+21
#nb.UnitMass_in_g = 1.9890E+40
#nb.UnitVelocity_in_cm_per_s = 1e5

nb.UnitLength_in_cm = 3.085e+21
nb.UnitMass_in_g = 4.4356e+44
nb.UnitVelocity_in_cm_per_s = 9.78469e+07

nb.boxsize = 2*rmax # Need to adapt shift in params.yml if this is changed
nb.write()