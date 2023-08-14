import numpy as np
from matplotlib import pyplot as plt
import argparse
import h5py

fontsz = 14
plt.rcParams.update({"text.usetex": True, "font.size": fontsz, "font.family": "serif"})
figsize = (7,5)

# Parse user input
parser = argparse.ArgumentParser(
    description="Plot multiple density profiles against theoretical prediction"
)
parser.add_argument("files", nargs="+", help="snapshot files to be imaged")
parser.add_argument("-Rmin", type=float, default=0.01, help="Min Radius")
parser.add_argument("-Rmax", type=float, default=10.0, help="Max Radius")
parser.add_argument("-nbins", type=int, default=64, help="Number of radii to sample (bins)")
parser.add_argument("-shift", type=float, default=10.0, help="Shift applied to particles in params.yml")
parser.add_argument("--saveplot", action='store_true')

args = parser.parse_args()
fnames = args.files

# Limit number of snapshots to plot
if len(fnames) > 20:
    raise ValueError(
        "Too many ({:d}) files provided (cannot plot more than 20).".format(len(fnames))
    )

rbins = np.logspace(np.log10(args.Rmin), np.log10(args.Rmax), args.nbins)

# Theoretical Density profile (change values of a and M if you used different ones to generate the IC)
M = 100
a = 1
def plummer_analytical(r):
    return (3.0*M/ (4.0 * np.pi * a ** 3)* (1.0 + r ** 2 / a ** 2) ** (-2.5))

t = np.empty(len(fnames))
mdens =   np.empty((len(fnames),args.nbins))

# Calculate & Plot Density profile
fig, ax = plt.subplots(figsize=figsize)
ax.set_xlabel("$r$ [kpc]",fontsize=fontsz+4)
ax.set_ylabel(r"$\rho(r)$",fontsize=fontsz+4)
ax.loglog()

ax.plot(rbins, plummer_analytical(rbins), c="magenta", label="Plummer Theoretical")

for j,fname in enumerate(fnames):
    f = h5py.File(fname, "r")
    pos = np.array(f["DMParticles"]["Coordinates"]) - args.shift
    vel = np.array(f["DMParticles"]["Velocities"])
    time = (f["Header"].attrs["Time"][0])
    t[j] = time
    mass = np.array(f["DMParticles"]["Masses"])
    r = np.sqrt(np.sum(pos ** 2, 1))

    # Methods to compute density profile
    def mass_ins(R):
        return ((r < R) * mass).sum()

    mass_ins_vect = np.vectorize(mass_ins)

    def density(R):
        return np.diff(mass_ins_vect(R)) / np.diff(R) / (4.0 * np.pi * R[1:] ** 2)

    dens = density(rbins)
    rs = rbins[1:]

    # remove empty bins
    c = dens > 0
    dens = np.compress(c, dens)
    rs = np.compress(c, rs)

    # Plot
    # ax.loglog(rsp[1:], density(rsp), "o", ms=1.7, label=r"$t=$ {:.3f} Gyr".format(time))
    if j == 0:
        ax.plot(rs, dens, label=r"$t=$ {:.1f}".format(t[j]),color='black',lw=2)
    else:
        ax.plot(rs, dens, label=r"$t=$ {:.1f}".format(t[j]),lw=1)


#ax.set_ylabel(r"$\rho(r)$ [$M_{\odot}$ kpc$^{-3}$]",fontsize=fontsz+4)

ax.legend()
if args.saveplot:
    fig.savefig("density.png",dpi=300,bbox_inches='tight')
else:
    plt.show()