import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib as mpl                      # Plotting library
from matplotlib import rc                     # Improve math fonts
import argparse
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
# rc('font',**{'family':'serif','serif':['Palatino']})
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}',r'\usepackage{siunitx}']
# matplotlib.rcParams['text.latex.preamble'] = [r'\renewcommand{\seriesdefault}{\bfdefault}',r'\boldmath']
#matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
# rc('mathtext', default='regular')
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
    fermi = None
    if "Ef" in Ef_line[1]:
      fermi = float(Ef_line[2])
      # print(fermi)
    f.readline().split()
    dim_line = f.readline().split()
    dimbs = None
    if "dimbs" in dim_line[1]:
      dimbs = int(dim_line[2])
      # print(dimbs)
  return npoints, name, point, fermi, dimbs

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
#
# IMPORTANT: If you want to plot the superconductor bands, then the way to
# execute this script is
#
# $ ipython < route to script.py> <route to datafile> -- --superconductivity
#
# If the user wants to have the electron and hole bands in the same plot
# then, the script has to be run as
# $ ipython < route to script.py> <route to datafile> -- --superconductivity -- together
#

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Parse bool")
    parser.add_argument("file", help="File to plot")
    parser.add_argument("--superconductivity", default=False, action="store_true" , help="Flag to plot the bands as for superconductors")
    parser.add_argument("--together", default=False, action="store_true" , help="Plot superconductor bands in the same plot")
    parser.add_argument("--zoom", default=0.0, type=float)
    parser.add_argument("--title", action="store", dest="title", default="")
    args = parser.parse_args()

    if args.superconductivity:
         if args.together:
             numplots = 1
             titles = [r"Band Structure"] #$\Delta = 0.001$ Ry"]
             # titles = [r"$\#_{k}=100$M, $\eta=5\times10^{-4}$, 2nn", r"$\#_{k}=100$M, $\eta=5\times10^{-4}$, 3nn"]
             # titles = [r"$\#_{k}=10$M", r"$\#_{k}=100$k"])
             # titles = [r"$\eta=5\times10^{-3}$", r"$\eta=5\times10^{-4}$"])

             fig, axs = plt.subplots(1, numplots, sharey=True, squeeze=False, figsize=(6*numplots, 5))
             axs[0,0].set_ylabel("Energy [Ry]")

             for i in range(numplots):
               npoints, name, point, fermi, dimbs = read_header(args.file)
               table = read_data(args.file)
               if args.title != "":
                   axs[0,i].set_title(args.title)
               else:
                   axs[0,i].set_title(titles[i])
               axs[0,i].set_xlim([point[0],point[npoints-1]])
               if args.zoom != 0.0:
                   axs[0,i].set_ylim(-args.zoom,args.zoom)
               # axs[0,i].set_ylim(table[:,1:].min(),table[:,1:].max())
               # axs[0,i].set_ylim(1.0/(table[:,1:].min()),1.0/(table[:,1:].max()))

               axs[0,i].set_xticks(point)
               axs[0,i].set_xticklabels(name)
               for j in point:
                 axs[0,i].axvline(x=j, color='k', linewidth=0.5)

               # Ploting the Fermi level or a line at y=0.0

               if (fermi == None): # susceptibility
                 axs[0,i].axhline(y=0.0, xmin=point[0], xmax=point[npoints-1], color='k', linestyle='-', linewidth=0.5)
               else: # band structure
                 axs[0,i].axhline(y=0.0, color='b', linestyle='--')

               # for j in enumerate(table[:,1]):
               #   if j[1] == point[4]:
               #     print 'Value at ',name[4],': ',-1.0/table[j[0],2]

               # Plotting the results
               if (fermi == None): # susceptibility
                 axs[0,i].plot(table[:,1],-1.0/table[:,2])
               else: # band structure
                     # axs[0,i].plot(table[:,0],table[:,(table.shape[1]-1)/2+1:]-fermi, color='r', linewidth=1.0, linestyle='-')
                     # axs[0,i].plot(table[:,0],table[:,1:(table.shape[1]-1)/2] + fermi, color='k', linewidth=1.0, linestyle='-')
                     axs[0,i].plot(table[:,0],table[:,(table.shape[1]-1)/2+1:], color='r', linewidth=1.0, linestyle='-')
                     axs[0,i].plot(table[:,0],table[:,1:(table.shape[1]-1)/2+1], color='k', linewidth=1.0, linestyle='-')

             plt.tight_layout()
             # plt.savefig("gaps2.png",dpi=900)
             plt.show()
         else:
             numplots = len(sys.argv)-1
             titles = [r"Electron Bands", r"Hole Bands"]
             # titles = [r"$\#_{k}=100$M, $\eta=5\times10^{-4}$, 2nn", r"$\#_{k}=100$M, $\eta=5\times10^{-4}$, 3nn"]
             # titles = [r"$\#_{k}=10$M", r"$\#_{k}=100$k"])
             # titles = [r"$\eta=5\times10^{-3}$", r"$\eta=5\times10^{-4}$"])

             fig, axs = plt.subplots(1, numplots, sharey=True, squeeze=False, figsize=(6*numplots, 5))
             axs[0,0].set_ylabel("Energy [Ry]")

             for i in range(numplots):
               npoints, name, point, fermi, dimbs = read_header(args.file)
               table = read_data(args.file)
               axs[0,i].set_title(titles[i])
               axs[0,i].set_xlim([point[0],point[npoints-1]])
               # axs[0,i].set_ylim(table[:,1:].min(),table[:,1:].max())
               # axs[0,i].set_ylim(1.0/(table[:,1:].min()),1.0/(table[:,1:].max()))

               axs[0,i].set_xticks(point)
               axs[0,i].set_xticklabels(name)
               for j in point:
                 axs[0,i].axvline(x=j, color='k', linewidth=0.5)

               # Ploting the Fermi level or a line at y=0.0

               if (fermi == None): # susceptibility
                 axs[0,i].axhline(y=0.0, xmin=point[0], xmax=point[npoints-1], color='k', linestyle='-', linewidth=0.5)
               else: # band structure
                 axs[0,i].axhline(y=0.0, xmin=point[0], xmax=point[npoints-1], color='k', linestyle='--')
                 if(i==0):
                     axs[0,i].axhline(y=fermi, xmin=point[0], xmax=point[npoints-1], color='b', linestyle='--')
                 else:
                     axs[0,i].axhline(y=-fermi, xmin=point[0], xmax=point[npoints-1], color='b', linestyle='--')

               # for j in enumerate(table[:,1]):
               #   if j[1] == point[4]:
               #     print 'Value at ',name[4],': ',-1.0/table[j[0],2]

               # Plotting the results
               if (fermi == None): # susceptibility
                 axs[0,i].plot(table[:,1],-1.0/table[:,2])
               else: # band structure
                 if(i==0):
                     axs[0,i].plot(table[:,0],table[:,(table.shape[1]-1)/2+1:], color='r', linewidth=1.0, linestyle='-')
                 else:
                     axs[0,i].plot(table[:,0],table[:,1:(table.shape[1]-1)/2], color='k', linewidth=1.0, linestyle='-')

             plt.tight_layout()
             plt.show()
    else:
         numplots = 1
         titles = [r"Band structure"]
         # titles = [r"$\#_{k}=100$M, $\eta=5\times10^{-4}$, 2nn", r"$\#_{k}=100$M, $\eta=5\times10^{-4}$, 3nn"]
         # titles = [r"$\#_{k}=10$M", r"$\#_{k}=100$k"])
         # titles = [r"$\eta=5\times10^{-3}$", r"$\eta=5\times10^{-4}$"])

         fig, axs = plt.subplots(1, numplots, sharey=True, squeeze=False, figsize=(6*numplots, 5))
         axs[0,0].set_ylabel("Energy [Ry]")

         for i in range(numplots):
           npoints, name, point, fermi, dimbs = read_header(args.file)
           table = read_data(args.file)
           axs[0,i].set_title(titles[i])
           axs[0,i].set_xlim([point[0],point[npoints-1]])
           # axs[0,i].set_ylim(table[:,1:].min(),table[:,1:].max())
           # axs[0,i].set_ylim(1.0/(table[:,1:].min()),1.0/(table[:,1:].max()))

           if args.title != "":
               axs[0,i].set_title(args.title)
           else:
               axs[0,i].set_title(titles[i])

           axs[0,i].set_xticks(point)
           axs[0,i].set_xticklabels(name)
           for j in point:
             axs[0,i].axvline(x=j, color='k', linewidth=0.5)

           if args.zoom != 0.0:
               axs[0,i].set_ylim(-args.zoom,args.zoom)

           # Ploting the Fermi level or a line at y=0.0
           if (fermi == None): # susceptibility
             axs[0,i].axhline(y=0.0, xmin=point[0], xmax=point[npoints-1], color='k', linestyle='-', linewidth=0.5)
           else: # band structure
             axs[0,i].axhline(y=fermi, xmin=point[0], xmax=point[npoints-1], color='k', linestyle='--')

           # for j in enumerate(table[:,1]):
           #   if j[1] == point[4]:
           #     print 'Value at ',name[4],': ',-1.0/table[j[0],2]

           # Plotting the results
           if (fermi == None): # susceptibility
             axs[0,i].plot(table[:,1],-1.0/table[:,2])
           else: # band structure
             axs[0,i].plot(table[:,0],table[:,1:dimbs], color='k', linewidth=1.0, linestyle='-')

         plt.tight_layout()
         plt.show()
