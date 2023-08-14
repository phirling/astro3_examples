import numpy as np
from matplotlib import pyplot as plt
import argparse
import h5py
import cycler

fontsz = 14
plt.rcParams.update({"text.usetex": True, "font.size": fontsz, "font.family": "serif"})
figsize = (7,5)

# Parse user input
parser = argparse.ArgumentParser(
    description="Plot multiple density profiles against theoretical prediction"
)
parser.add_argument("files", nargs="+", help="snapshot files to be imaged")
parser.add_argument("-Rmin", type=float, default=0.04, help="Min Radius")
parser.add_argument("-Rmax", type=float, default=10.0, help="Max Radius")
parser.add_argument("-nbins", type=int, default=32, help="Number of radii to sample (bins)")
parser.add_argument("-shift", type=float, default=10.0, help="Shift applied to particles in params.yml")
parser.add_argument("--saveplot", action='store_true')

args = parser.parse_args()
fnames = args.files

# Limit number of snapshots to plot
if len(fnames) > 20:
    raise ValueError(
        "Too many ({:d}) files provided (cannot plot more than 20).".format(len(fnames))
    )

# Set color cycle
ncolor = len(fnames)
color_cycle = plt.cm.plasma(np.linspace(0, 1,ncolor))
plt.rcParams['axes.prop_cycle'] = cycler.cycler('color', color_cycle)

# Estimate Characteristic time:
G = 1.
a = 1.0
b = 1.0
M = 1.0
rchar = a
v_a = np.sqrt(G*M*rchar**2 / ((a+b)**2 + rchar**2)**(3./2))
tdyn = 2*np.pi*rchar / v_a

rbins = np.logspace(np.log10(args.Rmin), np.log10(args.Rmax), args.nbins)

t = np.empty(len(fnames))
mdens =   np.empty((len(fnames),args.nbins))

# Calculate & Plot Density profile
fig, ax = plt.subplots(figsize=figsize)
ax.set_xlabel("Radius [kpc]",fontsize=fontsz+4)
ax.set_ylabel(r"Projected Surface Density",fontsize=fontsz+4)
ax.loglog()

for j,fname in enumerate(fnames):
    f = h5py.File(fname, "r")
    pos = np.array(f["DMParticles"]["Coordinates"]) - args.shift
    time = (f["Header"].attrs["Time"][0])
    t[j] = time / tdyn
    mass = np.array(f["DMParticles"]["Masses"])
    x = pos[:, 0]
    y = pos[:, 1]
    r = np.sqrt(pos[:,0]**2 + pos[:,1]**2)

    # Methods to compute density profile
    def mass_ins(R):
        return ((r < R) * mass).sum()

    mass_ins_vect = np.vectorize(mass_ins)

    def density(R):
        return np.diff(mass_ins_vect(R)) / np.diff(R) / (2.0 * np.pi * R[1:])

    dens = density(rbins)
    rs = rbins[1:]

    # remove empty bins
    c = dens > 0
    dens = np.compress(c, dens)
    rs = np.compress(c, rs)

    # Plot
    if j == 0:
        ax.plot(rs, dens, label=r"$t={:.1f} \ t_{{c}}$".format(t[j]),color='black',lw=2)
    else:
        ax.plot(rs, dens, label=r"$t={:.1f} \ t_{{c}}$".format(t[j]))

#ax.plot(rbins, plummer_analytical(rbins), c="darkturquoise",ls='--' ,lw=2.3,label="Plummer Theoretical")

ax.legend()
if args.saveplot:
    fig.savefig("density.png",dpi=300,bbox_inches='tight')
else:
    plt.show()