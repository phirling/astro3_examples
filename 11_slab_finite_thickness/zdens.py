import numpy as np
from matplotlib import pyplot as plt
import argparse
import h5py

"""
Script to plot the mass density against z, for the slab of finite thickness.
The analytical form of this density is given by rho(z) = rho0 * sech^2 (z/2z0),
(Eq. 4.302b in Binney & Tremaine). The script is used to test the stability of the
numerical model.
"""

# Parse user input
parser = argparse.ArgumentParser(description="Plot multiple density profiles against theoretical prediction")
parser.add_argument("files", nargs="+", help="snapshot files to be imaged")
parser.add_argument("-shift",type=float,default=0.5)
parser.add_argument("-z0",type=float,default=0.03)

args = parser.parse_args()
fnames = args.files

# Limit number of snapshots to plot
if len(fnames) > 20:
    raise ValueError(
        "Too many ({:d}) files provided (cannot plot more than 20).".format(len(fnames))
    )

# Figure Parameters
figsize = 6

# Model Parameters
G = 1.0
L = 1.                   # Box size
M = 1.                   # Total Mass
z0 = args.z0
rho0 = M/(4*L*L*z0)

# Z-bins to compute the mass density
nbins = 150
zsp = np.linspace(0,0.5,nbins)
dz = np.diff(zsp)[0]

# Analytical Density
def dens_analytical(z):
    return rho0 * np.cosh(0.5*z/z0)**(-2)

# Plot the analytical density
fig, ax = plt.subplots(figsize=(figsize*1.3, figsize))
ax.plot(zsp,dens_analytical(zsp),color='black')
ax.plot(-zsp,dens_analytical(-zsp),color='black')

# full z, negative to positive, in correct order, to plot density (exclude zeros)
fz = np.hstack((np.flip(-zsp[1:]),zsp[1:]))

# Plot numerical densities for all snapshots
for fname in fnames:
    print(fname)
    f = h5py.File(fname, "r")
    pos = np.array(f["DMParticles"]["Coordinates"]) - args.shift
    mass = np.array(f["DMParticles"]["Masses"])[0]
    z = pos[:, 2]
    sgns = np.sign(z)
    abv = sgns == 1
    blw = sgns == -1

    m_ins_abv = np.empty(nbins)
    m_ins_blw = np.empty(nbins)

    for k,zz in enumerate(zsp):
        mm1 = ((z[abv] < zz) * mass).sum()
        mm2 = ((-z[blw] < zz) * mass).sum()
        m_ins_abv[k] = mm1
        m_ins_blw[k] = mm2

    dm_abv = np.diff(m_ins_abv)
    dm_blw = np.diff(m_ins_blw)

    dens_abv = dm_abv / dz
    dens_blw = dm_blw / dz

    dens = np.hstack((np.flip(dens_blw),dens_abv))
    ax.plot(fz, dens,lw=0.8)

ax.set_xlabel("$z$ [Arbitrary Units]",fontsize=13)
ax.set_ylabel(r"$\rho(z)$ [Arbitrary Units]",fontsize=13)

plt.show()
#fig.savefig('plummerdens.eps',bbox_inches='tight')