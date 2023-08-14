import numpy as np
from matplotlib import pyplot as plt
import argparse
import h5py
import matplotlib.pyplot as plt

# Figure Parameters
fontsz = 14
plt.rcParams.update({"text.usetex": True, "font.size": fontsz, "font.family": "serif"})
figsize = (7,5)

# ===============================================================
# Plot the velocity dispersion (radial, tangential, total)
# of a spherically symmetric system.
# ===============================================================

# Parse user input
parser = argparse.ArgumentParser(
    description="Plot multiple density profiles against theoretical prediction"
)
parser.add_argument("files", nargs="+", help="snapshot files to be imaged")
parser.add_argument("-Rmin", type=float, default=0.0, help="Min Radius")
parser.add_argument("-Rmax", type=float, default=10.0, help="Max Radius")
parser.add_argument("-nbins", type=int, default=64, help="Number of radial bins")
parser.add_argument("-shift", type=float, default=10.0, help="Shift applied to particles in params.yml (to center the potential)")
parser.add_argument("--saveplot", action='store_true')

args = parser.parse_args()
fnames = args.files

# Define radial bins
rbins = np.linspace(args.Rmin,args.Rmax,args.nbins)

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
    iid = np.digitize(r,rbins,right=False)

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
fig1, ax1 = plt.subplots(figsize=figsize)
fig2, ax2 = plt.subplots(figsize=figsize)

for j in range(len(fnames)):
    idx = sigma[j,:] != -1
    if j == 0:
        ax1.plot(rbins[idx],sigma_r[j,idx],label=f"$t={t[j]:.1f}$",color='black',lw=2)
        ax2.plot(rbins[idx],sigma_t[j,idx],label=f"$t={t[j]:.1f}$",color='black',lw=2)
    else:
        ax1.plot(rbins[idx],sigma_r[j,idx],label=f"$t={t[j]:.1f}$")
        ax2.plot(rbins[idx],sigma_t[j,idx],label=f"$t={t[j]:.1f}$")
    #if args.v:  ax.plot(rbins[idx],sigma[j,idx],'-',label='Total')
    ax1.set_xlabel('$r$',fontsize=fontsz+4)
    ax2.set_xlabel('$r$',fontsize=fontsz+4)
    ax1.set_ylabel('$\sigma_r^2$',fontsize=fontsz+4)
    ax2.set_ylabel('$\sigma_t^2$',fontsize=fontsz+4)

ax1.set_title("Radial Velocity Dispersion")
ax2.set_title("Tangential Velocity Dispersion")

ax1.legend()
ax2.legend()

if args.saveplot:
    fig1.savefig("radial_velocity_dispersion.png",dpi=300,bbox_inches='tight')
    fig2.savefig("tangential_velocity_dispersion.png",dpi=300,bbox_inches='tight')
else:
    plt.show()