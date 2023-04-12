import numpy as np
from matplotlib import pyplot as plt
import argparse
from scipy.optimize import curve_fit
import h5py

# ===============================================================
# Plot the velocity dispersion (radial, tangential, total)
# of a Plummer model.
# When a single file is given, can plot all 3 components OR beta
# When multiple files are given, can only plot 1 component
# OR beta per file. Use -comp [both/r/t].
# ===============================================================

# Parse user input
parser = argparse.ArgumentParser(
    description="Plot multiple density profiles against theoretical prediction"
)
parser.add_argument("files", nargs="+", help="snapshot files to be imaged")
parser.add_argument("--notex", action="store_true", help="Flag to not use LaTeX markup")
parser.add_argument("-Rmin", type=float, default=0.0, help="Min Radius")
parser.add_argument("-Rmax", type=float, default=1.0, help="Max Radius")
parser.add_argument("-nbins", type=int, default=64, help="Number of radii to sample (bins)")
parser.add_argument("-shift", type=float, default=2.0, help="Shift applied to particles in params.yml")
parser.add_argument("--beta",action='store_true',help="Plot beta parameter rather than dispersions")
parser.add_argument("--log",action='store_true',help="Make loglog plot")
parser.add_argument("-comp",type=str,default='both',help="Component to Plot")

args = parser.parse_args()

fnames = args.files

# Set parameters
tex = not args.notex
if tex:
    plt.rcParams.update({"text.usetex": tex, "font.size": 22, "font.family": "serif"})
else:
    plt.rcParams.update({"font.size": 12})
figsize = 8

# Model Parameters (Mestel surface density)
G = 4.299581e04
#rsp = np.logspace(np.log10(args.Rmin), np.log10(args.Rmax), args.nbins)
rsp = np.linspace(args.Rmin,args.Rmax,args.nbins)

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
    time = (
        f["Header"].attrs["Time"][0]
        * f["Units"].attrs["Unit time in cgs (U_t)"][0]
        / 31557600.0e9
    )
    t[j] = time
    r = np.sqrt(np.sum((pos**2),axis=1))
    iid = np.digitize(r,rsp,right=False)
    for i in range(args.nbins):
        v = vel[iid == i]
        #vnorm = np.sqrt(np.sum(v**2,axis=1))
        pp = pos[iid == i]
        phi, theta = sph_angles(pp[:,0],pp[:,1],pp[:,2])

        v_r = v[:,0]*np.sin(theta)*np.cos(phi) + v[:,1]*np.sin(theta)*np.sin(phi) + v[:,2]*np.cos(theta)
        v_theta = v[:,0]*np.cos(theta)*np.cos(phi) + v[:,1]*np.cos(theta)*np.sin(phi) - v[:,2]*np.sin(theta)
        v_phi = -v[:,0]*np.sin(phi) + v[:,1]*np.cos(phi)
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

import prettyfancyplot as pfp

if not args.beta:
    fig, ax = pfp.makefig((7,7),xlabel='$r$',ylabel=r'$\sigma^2$')
    for j in range(len(fnames)):
        idx = sigma[j,:] != -1
        ax.plot(rsp[idx],sigma_r[j,idx],label='Radial')
        ax.plot(rsp[idx],sigma_t[j,idx],label='Tangential')
        ax.plot(rsp[idx],sigma[j,idx],label='Total')
        #ax.plot(rsp[idx],sigma_r[j,idx] + sigma_t[j,idx],label='Total')
    if args.log:
        ax.set_xscale('log')
        ax.set_yscale('log')
    ax.legend()
else:
    fig, ax = pfp.makefig((7,7),xlabel='$r$',ylabel=r'$\beta$')
    for j in range(len(fnames)):
        idx = sigma[j,:] != -1
        beta = 1. - sigma_t[j,idx]/(2*sigma_r[j,idx])
        ax.plot(rsp[idx],beta)
    if args.log:
        ax.set_xscale('log')
        ax.set_yscale('log')
pfp.show()