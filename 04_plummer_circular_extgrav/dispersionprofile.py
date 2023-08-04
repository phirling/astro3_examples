import numpy as np
from matplotlib import pyplot as plt
import argparse
from scipy.optimize import curve_fit
import h5py
import matplotlib.pyplot as plt

# Figure Parameters
#plt.rcParams.update({"text.usetex": True, "font.size": 22, "font.family": "serif"})
figsize = (7,7)

# ===============================================================
# Plot the velocity dispersion (radial, tangential, total)
# of a spherically symmetric system.
# ===============================================================

# Use G = 1 units
G = 1

# Parse user input
parser = argparse.ArgumentParser(
    description="Plot multiple density profiles against theoretical prediction"
)
parser.add_argument("files", nargs="+", help="snapshot files to be imaged")
parser.add_argument("-Rmin", type=float, default=0.0, help="Min Radius")
parser.add_argument("-Rmax", type=float, default=0.5, help="Max Radius")
parser.add_argument("-nbins", type=int, default=64, help="Number of radial bins")
parser.add_argument("-shift", type=float, default=2.0, help="Shift applied to particles in params.yml (to center the potential)")
parser.add_argument("--beta",action='store_true',help="Plot beta parameter rather than dispersions")
parser.add_argument("--log",action='store_true',help="Make loglog plot")
parser.add_argument("--v",action='store_true',help="Plot total dispersion")
parser.add_argument("--vr",action='store_true',help="Plot radial component of dispersion")
parser.add_argument("--vt",action='store_true',help="Plot tangential component of dispersion")

args = parser.parse_args()
fnames = args.files

# Exit if no flags haven been passed
if not (args.v or args.vr or args.vt):
    print("No components have been passed, exiting... (use --v, --vr, --vt flags)")
    quit()

# Define radial bins
rsp = np.linspace(args.Rmin,args.Rmax,args.nbins)

# Spherical angles of cartesian position
def sph_angles(x, y, z):
    hxy = np.hypot(x, y)
    theta = np.arctan2(hxy,z)
    phi = np.arctan2(y, x)
    return phi, theta

# Calculate Velocity Dispersion
t = np.empty(len(fnames))
sigma =   np.empty((len(fnames),args.nbins))
sigma_t = np.empty((len(fnames),args.nbins))
sigma_r = np.empty((len(fnames),args.nbins))
for j,fname in enumerate(fnames):
    f = h5py.File(fname, "r")
    pos = np.array(f["DMParticles"]["Coordinates"]) - args.shift
    vel = np.array(f["DMParticles"]["Velocities"])
    time = (f["Header"].attrs["Time"][0])
    t[j] = time

    # Compute radii and place particles in correct radial bins
    r = np.sqrt(np.sum((pos**2),axis=1))
    iid = np.digitize(r,rsp,right=False)

    # Loop over radial bins and compute dispersions in each bin
    for i in range(args.nbins):
        v = vel[iid == i]
        pp = pos[iid == i]
        phi, theta = sph_angles(pp[:,0],pp[:,1],pp[:,2])

        # Velocity components in spherical coordinates
        v_r = v[:,0]*np.sin(theta)*np.cos(phi) + v[:,1]*np.sin(theta)*np.sin(phi) + v[:,2]*np.cos(theta)
        v_theta = v[:,0]*np.cos(theta)*np.cos(phi) + v[:,1]*np.cos(theta)*np.sin(phi) - v[:,2]*np.sin(theta)
        v_phi = -v[:,0]*np.sin(phi) + v[:,1]*np.cos(phi)

        # If the bin is empty, set invalid value, if not, set dispersions
        if len(v) == 0:
            sigma_t[j,i] = -1
            sigma_r[j,i] = -1
            sigma[j,i] = -1
        else:
            sigma_t[j,i] = np.var(v_theta) + np.var(v_phi)
            sigma_r[j,i] = np.var(v_r)
            sigma[j,i] =   sigma_t[j,i] + sigma_r[j,i]

# ============
# Plot Results
# ============
fig, ax = plt.subplots(figsize=figsize)

for j in range(len(fnames)):
    idx = sigma[j,:] != -1
    if args.vr: ax.plot(rsp[idx],sigma_r[j,idx],'--',label='Radial')
    if args.vt: ax.plot(rsp[idx],sigma_t[j,idx],'-.',label='Tangential')
    if args.v:  ax.plot(rsp[idx],sigma[j,idx],'-',label='Total')
    ax.set_xlabel('$r$')
    ax.set_ylabel('$\sigma^2$')

if args.log:
        ax.set_xscale('log')
        ax.set_yscale('log')

ax.legend()
plt.show()