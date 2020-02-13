################################################################################
# This script with all its options is to be run like this
# ipython <script_location> <band_file> <ldosu_file> <ldosd_file> -- --superconductivity --mev --title="desired title"
# the arguments starting with "--" are optional and can be ommited
################################################################################

import numpy as np
import matplotlib.pyplot as plt
#from matplotlib2tikz import save as tikz_save
import matplotlib.gridspec as gridspec
import sys
from matplotlib import rc
import argparse
rc('mathtext', default='regular')

# Fonts
font = {'family': 'arial',
         'color': 'black',
         'weight':'normal',
         'size': 9,
        }

nSites=1

parser = argparse.ArgumentParser(description="Parse bool")
parser.add_argument("fileband", help="File to plot")
parser.add_argument("fileu", help="File to plot")
parser.add_argument("filed", help="File to plot")
parser.add_argument("--superconductivity", default=False, action="store_true" , help="Flag to plot the bands as for superconductors")
parser.add_argument("--mev", default=False, action="store_true" , help="Plot superconductor bands in the same plot")
parser.add_argument("--title", action="store", dest="title", default="")
parser.add_argument("--onlyS", default=False, action="store_true" , help="Plot only the s orbital in the LDOS")
parser.add_argument("--zoom", default=0.0, type=float)
args = parser.parse_args()

ry2ev = 13.6057

if args.mev:
    ry2ev = 13.6057*1000 # Conversion of energy units

if(ry2ev != 1.0):
  label = r'$E-E_F$ [eV]'
else:
  label = r'$E-E_F$ [Ry]'

if args.superconductivity:
    # colors = brewer2mpl.get_map('Set1', 'qualitative', 5).mpl_colors
    colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#a65628', '#b98600', '#000000']
    # colors = np.array(['k', 'g', 'r', 'b']) # The colors you wish to cycle through
    legends = np.array(['Total', r'$u_s$', r'$u_p$', r'$u_d$',r'$v_s$', r'$v_p$', r'$v_d$'])
else:
    # colors = brewer2mpl.get_map('Set1', 'qualitative', 5).mpl_colors
    colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#a65628']
    # colors = np.array(['k', 'g', 'r', 'b']) # The colors you wish to cycle through
    legends = np.array(['Total', 's', 'p', 'd'])

bsstruct = args.fileband

ldosu = args.fileu
ldosd = args.filed
filename = 'BS_LDOS.pdf'
fig, ax = plt.subplots(1,3, sharey=True, gridspec_kw = {'width_ratios':[1,4,1]})
fig.subplots_adjust(wspace=0.15)

if args.title != "":
    fig.suptitle(args.title,fontsize=24)

ax[0].tick_params(axis='y', direction='in', left=True, right=True)
ax[0].set_ylabel(label, fontsize=14)
ax[2].tick_params(axis='y', direction='in', left=True, right=True)


fermi = 0.0

def read_header(file):
  with open(file, "r") as f:
    count = [s for s in f.readline().split()]
    npoints = int(count[1])
    # print npoints
    name = np.empty(npoints, dtype=str)
    point = np.empty(npoints, dtype=float)
    for i in range(npoints):
      name[i], point[i] = f.readline().split()[1:]
      # print name[i], point[i]

    Ef_line = f.readline().split()
    fermi = None
    if "Ef" in Ef_line[1]:
      fermi = float(Ef_line[2])
      # print fermi
  return npoints, name, point, fermi

def read_data(filename):
  # open a file using with statement
  data = []
  with open(filename,'r') as file:
    for line in file:
      # check if the current line
      # starts with "#"
      if not (line.startswith("#") or line.startswith(" #")):
        data.append([float(x) for x in line.split()])
  ndata = np.array(data)
  return ndata

def read_header_ldos(file):
   with open(file, "r") as f:
     Ef_line = f.readline().split()
     fermi = None
     if "Ef" in Ef_line[1]:
       fermi = float(Ef_line[2])
       # print fermi
   return fermi

def read_data_ldos(filename):
  # open a file using with statement
  data = []
  with open(filename,'r') as file:
    for line in file:
      # check if the current line
      # starts with "#"
      if not (line.startswith("#") or line.startswith(" #")):
        data.append([float(x) for x in line.split()])
  ndata = np.array(data)
  return ndata

npoints, name, point, fermi = read_header(args.fileband)
table = read_data(args.fileband)
fermi_ldos = read_header_ldos(args.fileu)

ndatau=read_data_ldos(args.fileu)
ndatad=read_data_ldos(args.filed)

ndatau = ndatau[ndatau[:,0].argsort()]
ndatad = ndatad[ndatad[:,0].argsort()]

ax[1].set_xticks(point)
ax[1].set_xticklabels(name)
for i in point:
    ax[1].axvline(x=i, color='k', linewidth=0.75)
ax[1].axhline(y=0.0, xmin=point[0], xmax=point[npoints-1], color='k', linestyle='--', linewidth=0.75)

if (fermi == None): # susceptibility
  ax[1].plot(table[:,1],-1.0/table[:,2])
else: # band structure
  if args.superconductivity:
      ax[1].plot(table[:,0],table[:,(table.shape[1]-1)/2+1:]*ry2ev-fermi_ldos*ry2ev, color='r', linewidth=1.0, linestyle='-')
      ax[1].plot(table[:,0],table[:,1:(table.shape[1]-1)/2+1]*ry2ev-fermi_ldos*ry2ev, color='k', linewidth=1.0, linestyle='-')
  else:
      ax[1].plot(table[:,0],table[:,1:]*ry2ev-fermi_ldos*ry2ev, color='k', linewidth=1.0, linestyle='-')

ax[1].tick_params(axis='y', direction='in', left=True, right=True)

ax[1].set_xlim([point[0],point[npoints-1]])

x = ndatau[:,0]

if args.onlyS:
    ax[0].plot(-ndatau[:,2]/ry2ev,(x-fermi_ldos)*ry2ev,linestyle='-', color=colors[2-1], label=legends[2-1],marker=1,markersize=1)
    if args.superconductivity:
        ax[0].plot(-ndatau[:,5]/ry2ev,(x-fermi_ldos)*ry2ev,linestyle='-', color=colors[5-1], label=legends[5-1],marker=1,markersize=1)
        ax[0].plot(-ndatau[:,5]/ry2ev - ndatau[:,2]/ry2ev,(x-fermi_ldos)*ry2ev,linestyle='-', color=colors[1-1], label="sum",marker=1,markersize=1)
else:
    for i in range(1,len(ndatau[0,:])):
        # ax[0].plot(-ndatau[:,i]/ry2ev,(x-fermi_ldos)*ry2ev,linestyle='-', color=colors[i-1], label=legends[i-1],marker=1,markersize=1)
        ax[0].plot(-ndatau[:,i]/ry2ev,(x-fermi_ldos)*ry2ev,linestyle='-', color=colors[i-1], label=legends[i-1],marker=1,markersize=1)



x = ndatad[:,0]
if args.onlyS:
    ax[2].plot(ndatad[:,2]/ry2ev,(x-fermi_ldos)*ry2ev,linestyle='-', color=colors[2-1], label=legends[2-1],marker=1,markersize=1)
    if args.superconductivity:
        ax[2].plot(ndatad[:,5]/ry2ev,(x-fermi_ldos)*ry2ev,linestyle='-', color=colors[5-1], label=legends[5-1],marker=1,markersize=1)
        ax[2].plot(ndatad[:,5]/ry2ev + ndatad[:,2]/ry2ev,(x-fermi_ldos)*ry2ev,linestyle='-', color=colors[1-1], label="sum",marker=1,markersize=1)
else:
    for i in range(1,len(ndatad[0,:])):
        # ax[2].plot(ndatad[:,i]/ry2ev,(x-fermi_ldos)*ry2ev,linestyle='-', color=colors[i-1],marker=1,markersize=1)
        ax[2].plot(ndatad[:,i]/ry2ev,(x-fermi_ldos)*ry2ev,linestyle='-', color=colors[i-1],marker=1,markersize=1)

#
a1 = max(ndatau[:,1]/ry2ev)
a2 = max(ndatad[:,1]/ry2ev)

max_ldos = max(a1,a2)
xlim = 1.1*abs(max_ldos)
if args.superconductivity:
    # ax[2].set_xlim([0.0,xlim])
    # ax[0].set_xlim([-xlim,0.0])
    if args.zoom != 0.0 :
        # ax[2].set_xlim([0.0,xlim])
        # ax[0].set_xlim([-1.5,0.0])
        ax[1].set_ylim([-args.zoom,args.zoom])
    else:
        ax[1].set_ylim([(x-fermi_ldos)[0]*ry2ev,(x-fermi_ldos)[-1]*ry2ev])
else:
    ax[2].set_xlim([0.0,xlim])
    ax[0].set_xlim([-xlim,0.0])
    ax[1].set_ylim([(x-fermi_ldos)[0]*ry2ev,(x-fermi_ldos)[-1]*ry2ev])
#
# ax[2].axhline(y=1000.0, xmin=point[0], xmax=point[npoints-1], color='b', linestyle='--', linewidth=0.75)
# ax[2].axhline(y=-1000.0, xmin=point[0], xmax=point[npoints-1], color='b', linestyle='--', linewidth=0.75)
#
# ax[0].axhline(y=1000.0, xmin=point[0], xmax=point[npoints-1], color='b', linestyle='--', linewidth=0.75)
# ax[0].axhline(y=-1000.0, xmin=point[0], xmax=point[npoints-1], color='b', linestyle='--', linewidth=0.75)

ax[2].axhline(y=0.0, xmin=point[0], xmax=point[npoints-1], color='k', linestyle='--', linewidth=0.75)
ax[0].axhline(y=0.0, xmin=point[0], xmax=point[npoints-1], color='k', linestyle='--', linewidth=0.75)
ax[0].legend(loc=2, prop={'size': 12})


ax[0].tick_params(axis='both', which='major', labelsize=14)
ax[1].tick_params(axis='both', which='major', labelsize=14)
ax[2].tick_params(axis='both', which='major', labelsize=14)


plt.savefig(filename)

plt.show()
