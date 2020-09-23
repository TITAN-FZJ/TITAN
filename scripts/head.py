################################################################################
# With this file we look to unify the styles of the plots on TITAN. If the
# plot requires specific options those will go in their respective files
################################################################################

import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib as mpl                      # Plotting library
from matplotlib import rc                     # Improve math fonts
import matplotlib.colors as colors            # Color selection and manipulation
from scipy import integrate
import argparse
import pandas as pd                           # Python Data Analysis Library
import matplotlib.gridspec as gridspec

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
# rc('font',**{'family':'serif','serif':['Palatino']})
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}',r'\usepackage{siunitx}']
# Default fonts
mpl.rcParams['font.size']        = 12
mpl.rcParams['font.family']      = 'serif'
mpl.rcParams['figure.titlesize'] = 'large'
mpl.rcParams['axes.titlepad'] = 10
mpl.rcParams['lines.linewidth']  = 2
#Legends:
mpl.rcParams['legend.fontsize']  = 'medium'
mpl.rcParams['legend.fancybox'] = False
# rcParams['legend.loc'] = 'upper left'`
mpl.rcParams['legend.framealpha'] = None
mpl.rcParams['legend.edgecolor'] = 'inherit'
mpl.rcParams['legend.handlelength'] = 2
mpl.rcParams["font.weight"] = "bold"
mpl.rcParams["axes.labelweight"] = "bold"

mpl.rcParams["xtick.major.width"] = 1.5
mpl.rcParams["xtick.minor.width"] = 1.0
mpl.rcParams["ytick.major.width"] = 1.5
mpl.rcParams["ytick.minor.width"] = 1.0
linewidth = 1.5
colors = [(0.07, 0.19, 0.60),(0.57, 0.05, 0.07),(0.10, 0.46, 0.13)]

parser = argparse.ArgumentParser(description="Parse bool")
parser.add_argument("files", nargs="+", help="File(s) to plot")
parser.add_argument('--output', help='Output file', default="output.png")
parser.add_argument('--title' , help='Title of the graph', default="")
parser.add_argument('--xlim' , help='Limit of x-axis', default="")
parser.add_argument('--ylim' , help='Limit of y-axis', default="")
parser.add_argument('--domain' , help='String for the x-range given as \"star,end\" ', default="")
parser.add_argument("--mev", default=False, action="store_true" , help="use meV")
parser.add_argument("--ev", default=False, action="store_true" , help="use eV")
parser.add_argument('--orbitals', help='Array of integers', default="[]")
parser.add_argument('--legends' , help='Array of legends', default="[]")
parser.add_argument('--total' , default=False, action="store_true" , help='Plot the total DOS')
parser.add_argument('--integrated' , default=False, action="store_true" , help='Plot the integrated DOS')
parser.add_argument("--superconductivity", default=False, action="store_true" , help="Flag to plot the bands as for superconductors")
parser.add_argument("--onlyS", default=False, action="store_true" , help="Plot only the s orbital in the LDOS")
parser.add_argument("--zoom", default=0.0, type=float)
parser.add_argument("--together", default=False, action="store_true" , help="Plot superconductor bands in the same plot")
parser.add_argument("--noef", default=True, action="store_false" , help="Do not plot the line at Ef")
parser.add_argument("--separate", default=False, action="store_true" , help="Plot electron and hole bands in different plots")
parser.add_argument("--guide", default=False, action="store_true" , help="Guide for band structures")
parser.add_argument('--project', help='Array of integers', default="[[1,2,3,4,5,6,7,8,9],[10,11,12,13,14,15,16,17,18]]")
parser.add_argument('--show' , default=False, action="store_true" , help='Show the plot at the end')
# In case we need to add arguments in the specific files
# args, options = parser.parse_known_args()
args = parser.parse_args()

ry2ev = 1.0

if args.mev:
    ry2ev = 13.6057*1000 # Conversion of energy units

if args.ev:
    ry2ev = 13.6057 # Conversion of energy units
