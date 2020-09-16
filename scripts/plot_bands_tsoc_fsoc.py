# import numpy as np
# import matplotlib.pyplot as plt
# import pandas as pd                           # Python Data Analysis Library
# import sys
# import matplotlib as mpl                      # Plotting library
# from matplotlib import rc                     # Improve math fonts
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# ## for Palatino and other serif fonts use:
# # rc('font',**{'family':'serif','serif':['Palatino']})
# mpl.rcParams['text.usetex'] = True
# mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}',r'\usepackage{siunitx}']
# # matplotlib.rcParams['text.latex.preamble'] = [r'\renewcommand{\seriesdefault}{\bfdefault}',r'\boldmath']
# #matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
# # rc('mathtext', default='regular')
# # Default fonts
# mpl.rcParams['font.size']        = 12
# mpl.rcParams['font.family']      = 'Arial'
# mpl.rcParams['figure.titlesize'] = 'large'
# mpl.rcParams['axes.titlepad'] = 10
# mpl.rcParams['lines.linewidth']  = 2
# #Legends:
# mpl.rcParams['legend.fontsize']  = 'medium'
# mpl.rcParams['legend.fancybox'] = False
# # rcParams['legend.loc'] = 'upper left'`
# mpl.rcParams['legend.framealpha'] = None
# mpl.rcParams['legend.edgecolor'] = 'inherit'
# mpl.rcParams['legend.handlelength'] = 2
# mpl.rcParams["font.weight"] = "bold"
# mpl.rcParams["axes.labelweight"] = "bold"

################################################################################
# Import the arguments and style
################################################################################
from head import *


################################################################################
# Get the header from the file
################################################################################
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
    fermi = 0.0
    if "Ef" in Ef_line[1]:
      fermi = float(Ef_line[2])
      # print fermi
  return npoints, name, point, fermi


################################################################################
# Get the data from the file and save it into a matrix
################################################################################
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



################################################################################
# Main program
################################################################################
if __name__ == "__main__":
  numplots = 2
  titles = ["Band Structure w/o SOI", "Band Structure with SOI"]

  fig, axs = plt.subplots(1, numplots, sharey=True, squeeze=False, figsize=(5*numplots, 5))
  axs[0,0].set_ylabel("Energy [Ry]")

  for i in range(numplots):
    npoints, name, point, fermi = read_header(args.files[i])
    table = read_data(args.files[i])
    axs[0,i].set_title(titles[i])
    axs[0,i].set_xlim([point[0],point[npoints-1]])
    # axs[0,i].set_ylim(table[:,1:].min(),table[:,1:].max())
    axs[0,i].set_ylim(table[:,1:].min(),fermi+0.2)

    axs[0,i].set_xticks(point)
    axs[0,i].set_xticklabels(name)
    for j in point:
      axs[0,i].axvline(x=j, color='k', linewidth=0.5)

    # Ploting the Fermi level or a line at y=0.0
    if not (fermi == None):
      axs[0,i].axhline(y=fermi, xmin=point[0], xmax=point[npoints-1], color='k', linestyle='--')
    else:
      axs[0,i].axhline(y=0.0, xmin=point[0], xmax=point[npoints-1], color='k', linestyle='-', linewidth=0.5)

    for j in range(1,10):
      axs[0,i].scatter(table[:,0],table[:,j], marker='.', c='r', s=0.1)
    for j in range(10,19):
      axs[0,i].scatter(table[:,0],table[:,j], marker='.', c='b', s=0.1)


  plt.tight_layout()
  plt.show()
