import numpy as np
from matplotlib import pyplot as plt
import argparse
import h5py
import matplotlib.pyplot as plt

# Figure Parameters
plt.rcParams.update({"text.usetex": True, "font.size": 22, "font.family": "serif"})
figsize = (7,7)

# ===============================================================
# Plot orbit of an individual particle over time
# ===============================================================

# Parse user input
parser = argparse.ArgumentParser(description="Plot orbit of an individual particle over time")
parser.add_argument("files", nargs="+", help="snapshot files")
parser.add_argument("-lim",type=float,default=0.5,help="Limit of the figure in internal units")
parser.add_argument("-shift", type=float, default=2.0, help="Shift applied to particles in params.yml (to center the potential)")
parser.add_argument("-np",default=0,help="Index of the particle to follow (must be < than number of particles in simulation)")
parser.add_argument("-axis",type=int,default=2,help="Axis from which to view orbit (default: 2 == z axis == view xy plane)")

args = parser.parse_args()
fnames = args.files

particle_number = int(args.np)
orbit = np.empty((3,len(fnames)))

# Extract orbit of particle
for j,fname in enumerate(fnames):
    f = h5py.File(fname, "r")
    pos = np.array(f["DMParticles"]["Coordinates"]) - args.shift
    orbit[:,j] = pos[particle_number,:]

# ============
# Plot Results
# ============
fig, ax = plt.subplots(figsize=figsize)
axnames = ['$x$','$y$','$z$']

if args.axis == 0:
    i1 = 1
    i2 = 2
elif args.axis == 1:
    i1 = 0
    i2 = 2
else:
    i1 = 0
    i2 = 1
ax.plot(orbit[i1],orbit[i2])
ax.set_xlim(-args.lim,args.lim)
ax.set_ylim(-args.lim,args.lim)
ax.set_xlabel(axnames[i1])
ax.set_ylabel(axnames[i2])

plt.show()