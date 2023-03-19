import numpy as np
from matplotlib import pyplot as plt
import argparse
import h5py

# Parse user input
parser = argparse.ArgumentParser(description="Plot multiple density profiles against theoretical prediction")
parser.add_argument("files", nargs="+", help="snapshot files to be imaged")
parser.add_argument("-shift",type=float,default=0.0)

args = parser.parse_args()
fnames = args.files

# Limit number of snapshots to plot
if len(fnames) > 20:
    raise ValueError(
        "Too many ({:d}) files provided (cannot plot more than 20).".format(len(fnames))
    )

# Set parameters
figsize = 8

# Model Parameters (Mestel surface density)
G = 1.0
nbins = 100
zsp = np.linspace(0,0.5,nbins)
dz = np.diff(zsp)[0]

z0 = 0.03
rho0 = 8.333333333333334
def dens_analytical(z):
    return rho0 * np.cosh(0.5*z/z0)**(-2)

# Plot densities
fig, ax = plt.subplots(figsize=(figsize, figsize))
ax.plot(zsp,dens_analytical(zsp),color='black')
for fname in fnames:
    print(fname)
    f = h5py.File(fname, "r")
    pos = np.array(f["DMParticles"]["Coordinates"]) - args.shift
    time = (
        f["Header"].attrs["Time"][0]
        * f["Units"].attrs["Unit time in cgs (U_t)"][0]
        / 31557600.0e9
    )
    mass = np.array(f["DMParticles"]["Masses"])
    z = pos[:, 2]
    print(z)
    m_ins = np.empty(nbins)
    rho_z = np.empty(nbins)
    for k,zz in enumerate(zsp):
        mm = ((np.abs(z) < zz) * mass).sum()
        m_ins[k] = mm/2
    dm = np.diff(m_ins)
    dens = dm / dz
    ax.plot(zsp[1:], dens,lw=0.8)

ax.set_xlabel("$z$",fontsize=25)

#ax.set_ylabel(r"$\rho(r)$ [$M_{\odot}$ kpc$^{-3}$]",fontsize=25)
plt.show()
#fig.savefig('plummerdens.eps',bbox_inches='tight')