#!/usr/bin/env python

###########################################################################
# Copyright (c) 2023 Yves Revaz & Patrick Hirling
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###########################################################################


import numpy as np
import argparse
import h5py
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize, LogNorm, PowerNorm
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import pickle as pkl
from tqdm import tqdm
from pNbody.libutil import set_ranges
from pNbody import Nbody

# ================
# Parse User input
# ================

parser = argparse.ArgumentParser(description="Extract data from SWIFT snapshot(s) and create film/image")

# Snapshot Parameters
parser.add_argument("files",nargs='+',help="snapshot files to be imaged")
parser.add_argument("-shift",nargs='*', type=float, default=[0.0], help="Shift applied to particles in params.yml [Either a single value or triple (x y z)]")
parser.add_argument("--openextract",action='store_true',help="Open a previous pickled extract file (in this case the first argument is interpreted as the file)")
parser.add_argument("--saveextract",action='store_true',help="Save the extracted data as a pickle file")

# Histogram Parameters
parser.add_argument("-nbins",type=int,default=400,help="Number of bins in each dimension (x,y)")
parser.add_argument("-lim",type=float,default=0.5,help="Limits of histogram (in kpc)")
parser.add_argument("-view",type=str,default='xy',help="Which plane to image")

# Image Parameters (Color, normalization,etc)
parser.add_argument("-cmap",type=str,default='YlGnBu_r',help="Colormap (try YlGnBu_r, Magma !)")
parser.add_argument("-cmin",type=float,default=None,help=
                    """Minimum Physical value in the Histogram.
                    This effectively sets the contrast of the images. Default: data minimum""")
parser.add_argument("-cmax",type=float,default=None,help="Maximum physical value in the Histogram. Default: data maximum")
parser.add_argument("-interp",type=str,default='none',help="Interpolation used ('none','kaiser','gaussian',...). Default: none")

# Figure Parameterrs
parser.add_argument("--notex",action='store_true')
parser.add_argument("-figheight",type=float,default=8,help="Height of figure in inches")
parser.add_argument("-figwidth",type=float,default=8,help="Width of figure in inches")
parser.add_argument("--savefilm",action='store_true',help="Save film/image")
parser.add_argument("--saveframes",action='store_true',help="Save as frames rather than film")
parser.add_argument("--nolabels",action='store_true',help="Do not display axes, ticks, labels")
parser.add_argument("--nounits",action='store_true',help="Do not display axis and time units")

# Encoding Parameters
parser.add_argument("-fps",type=float,default=24,help="Frames per second")
parser.add_argument("--headless",action='store_true',help="Flag to run on headless server, uses non-GUI Agg backend")
parser.add_argument("-bitrate",type=float,default=-1,help="Birate to use for ffmpeg")
parser.add_argument("-dpi",type=float,default=300,help="DPI if an image is saved")

# Corotation Parameters
parser.add_argument("-rcr",type=float,default=0.0,
    help="Radius of the corotating origin, set to 0 for no corotation (default)")
parser.add_argument("-acr",type=float,default=0,
    help="Initial angle of the corotating origin")
parser.add_argument("-vcr",type=float,default=200,
    help="Circular velocity of the corotation")
parser.add_argument("--set_origin",action='store_true',
    help="Set the origin of the figure to the corotating origin (default: no)")

args = parser.parse_args()


# =====================================
# Define Functions and useful variables
# =====================================
# Convenience
fnames = args.files
lim = args.lim
sh = args.shift

# Corotation
corot = args.rcr > 0
if corot:
    omega = -args.vcr / args.rcr
    x0, y0 = 0.0,0.0
    if args.set_origin:
        x0 = args.rcr*np.cos(args.acr)
        y0 = args.rcr*np.sin(args.acr)
    def corotate_pos(X,Y,t):
        x = np.cos(omega*t)*X - np.sin(omega*t)*Y - x0
        y = np.sin(omega*t)*X + np.cos(omega*t)*Y - y0
        return x,y
    
# Shift
if len(sh) == 1:
    shift = np.array([sh[0],sh[0],sh[0]])
elif len(sh) == 3:
    shift = np.array([sh[0],sh[1],sh[2]])
else:
    raise ValueError("Shift can be scalar or 3-array")

# TODO: adapt
scale   = "log"
mn      = None
mx      = None
cd      = None
frsp    = 0 #2.5
params = {}
params['size']  	= (lim,lim)
params['shape'] 	= (args.nbins,args.nbins)
params['rendering']	= 'map'
params['obs'] 	= None
params['xp']  	= None
params['view']	= args.view
params['mode']	= 'm'
params['exec']	= None
params['macro']	= None
params['frsp']	= frsp
params['filter_name'] = None 

# ============================================================
# Extract Data from the snapshot(s) or load a previous extract
# ============================================================
if not args.openextract:
    nframes = len(fnames)
    histdata = np.zeros((nframes,args.nbins,args.nbins),dtype=np.uint32)
    times = np.empty(nframes)
    print("Extracting Data...")
    for i,fn in enumerate(tqdm(fnames)):
        nb = Nbody(fn,ftype='swift')
        nb.translate(-shift)
        hh = 1e6*nb.CombiMap(params) # Multiply by 1M to ensure that close to 0 values are saved correctly
        histdata[i] = hh
        times[i] = nb.time

    print("done.")

    # If requested, save extracted histogram data to pickle file
    if args.saveextract:
            print("Saving Extracted Data...")
            infotext = "Time in Myrs, lim in kpc, 'hist' obj is a (nb_frames * nbins_x * nbins_y) array of int32"
            output = {
                'distunit': 'kpc',
                'timeunit': 'Myr',
                'lim' : lim,
                'info' : infotext,
                't' : times,
                'hist' : histdata,
                'rcr' : args.rcr,
                'acr' : args.acr,
                'vcr' : args.vcr,
                'corot_origin' : args.set_origin
            }
            with open('extract.pkl','wb') as f:
                pkl.dump(output,f)

else:
    if len(fnames) == 1:
        print("Reading "+fnames[0]+"...")
        with open(fnames[0],'rb') as f:
            extract = pkl.load(f)
        histdata = extract['hist']
        nframes = histdata.shape[0]
        lim = extract['lim']
        times = extract['t']
    else:
        print(fnames)
        raise RuntimeError("Can only open one extract")
    
# ====================
# Create Image or Film
# ====================

# Configure colour & normalization: by default, normalize to min,max value across ALL snapshots
if args.cmin is not None:
    vmin = args.cmin
else:
    vmin = np.amin(histdata)
if args.cmax is not None:
    vmax = args.cmax
else:
    vmax = np.amax(histdata)

if vmin >= vmax:
    raise ValueError("vmin ({:n}) is >= vmax ({:n})".format(vmin,vmax))

# Configure Matplotlib
if args.headless:
    import matplotlib
    matplotlib.use('Agg')
if not args.notex: plt.rcParams.update({"text.usetex": True,'font.size':2*args.figheight,'font.family': 'serif'})
else: plt.rcParams.update({'font.size':15})


# Construct figure
fig, ax = plt.subplots(figsize=(args.figwidth,args.figheight))

# Units
if args.nounits:
    time_unitstr = ""
    space_unitstr = ""
else:
    time_unitstr = " Myr"
    space_unitstr = " [kpc]"

# Define imshow extent
ext = (-lim,lim,-lim,lim)


# Set axes limits and add labels if desired
ax.set_xlim(-lim,lim)
ax.set_ylim(-lim,lim)
ax.set_aspect('equal')
if args.nolabels:
    ax.set_xticks([])
    ax.set_yticks([])
else:
    ax.set_xlabel(args.view[0] + space_unitstr)
    ax.set_ylabel(args.view[1] + space_unitstr)


# Initialize Image
im = ax.imshow(np.zeros((args.nbins,args.nbins)),
                interpolation = args.interp,
                vmin=0,vmax=255,
                extent = ext, cmap = args.cmap,origin='lower')

# Initialize Title
ttl = ax.text(0.01, 0.99,"$t={:.2f}$".format(0)+ time_unitstr,
    horizontalalignment='left',
    verticalalignment='top',
    color = 'white',
    transform = ax.transAxes)

# If a single file is passed, fill image and show
if nframes == 1:
    hdat,mn_opt,mx_opt,cd_opt = set_ranges(histdata[0],scale=scale,cd=cd,mn=vmin,mx=vmax)
    im.set_data(hdat)

# If multiple files are passed and a dynamical output is wished, launch animation
elif not args.saveframes:
    # Closure function to use global image & title object with blitting
    def prepare_anim(im,title):
        def update(i):    
            # Title
            title.set_text("$t={:.2f}$".format(times[i])+ time_unitstr)
        
            # Update image
            hdat,mn_opt,mx_opt,cd_opt = set_ranges(histdata[i],scale=scale,cd=cd,mn=vmin,mx=vmax)
            im.set_data(hdat)
            return im, title

        return update

    # Show encoding progress if --savefilm
    if args.savefilm:
        print("Encoding...")
        itrtr = tqdm(range(nframes))
    else:
        itrtr = range(nframes)
    # Animation
    ani = FuncAnimation(fig, prepare_anim(im,ttl),frames = itrtr,
                    blit=True)#cache_frame_data=False)

# If multiple files are passed but a per-frame output is wished, fill image and save png for each frame
else:
    for i,h in enumerate(histdata):
        hdat,mn_opt,mx_opt,cd_opt = set_ranges(h,scale=scale,cd=cd,mn=vmin,mx=vmax)
        im.set_data(hdat)
        ttl.set_text("$t={:.2f}$".format(times[i])+ time_unitstr)
        outfn = "frame_{:04n}.png".format(i)
        fig.savefig(outfn,dpi=args.dpi,bbox_inches='tight')

# Show or save. TODO: distribute above for nicer structure
if not args.savefilm and not args.saveframes:
    plt.show()
else:
    if nframes == 1:
        fig.savefig("image.png",dpi=args.dpi,bbox_inches='tight')
    else:
        ani.save('film.mp4',writer='ffmpeg',fps=args.fps,bitrate=args.bitrate)