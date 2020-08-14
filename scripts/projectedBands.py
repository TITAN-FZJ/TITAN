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

#Ticks:
mpl.rcParams["xtick.major.width"] = 1.5
mpl.rcParams["xtick.minor.width"] = 1.0
mpl.rcParams["ytick.major.width"] = 1.5
mpl.rcParams["ytick.minor.width"] = 1.0
linewidth = 1.5
colors = [(0.07, 0.19, 0.60),(0.57, 0.05, 0.07),(0.10, 0.46, 0.13)]

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
    name_out = []
    for i in range(npoints):
      name[i], point[i] = f.readline().split()[1:]
      if name[i] == "G":
        name_out.append(r"$\Gamma$")
      else:
        name_out.append(name[i])
      # print name[i], point[i]

    Ef_line = f.readline().split()
    fermi = None
    if "Ef" in Ef_line[1]:
      fermi = float(Ef_line[2])
      # print fermi

    ##### jump garbage line
    f.readline()

    dimbs_line = f.readline().split()
    dimbs = int(dimbs_line[2])

  return npoints, name_out, point, fermi, dimbs

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
# In the next example we show how to use all options, we use \\ to break the line,
# and make it easier to read
#
# ipython  < route to script.py> <route to datafile> -- \\
#           --superconductivity \\
#           --title="Projection of \$p_z\$" \\
#           --ev \\
#           --guide \\
#           --output="cool name" \\
#           --zoom=0.1
#           --project=[[1,19],[2,20],[4,24]]
#           --separate
#
# --separate does not work together with all the rest for the moment

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="Parse bool")
  parser.add_argument("file", help="File to plot")
  parser.add_argument("--superconductivity", default=False, action="store_true" , help="Flag to plot the bands as for superconductors")
  parser.add_argument("--output", default="", action="store" , help="Output filename")
  parser.add_argument("--separate", default=False, action="store_true" , help="Plot superconductor bands in the same plot")
  parser.add_argument("--zoom", default=0.0, type=float)
  parser.add_argument("--ev", default=False, action="store_true" , help="Plot superconductor bands in the same plot")
  parser.add_argument("--noef", default=True, action="store_false" , help="Do not plot the line at Ef")
  parser.add_argument("--title", action="store", dest="title", default="")
  parser.add_argument("--guide", default=False, action="store_true" , help="Guide for band structures")
  parser.add_argument('--project', help='Array of integers', default="[[1,2,3,4,5,6,7,8,9],[10,11,12,13,14,15,16,17,18]]")
  args = parser.parse_args()

  if args.ev == True:
    ry2ev = 13.605662285137
  else:
    ry2ev = 1.0

  if args.superconductivity:
    if args.separate:
      numplots = len(sys.argv)-1
      titles = [r"Electron Bands", r"Hole Bands"]

      fig, axs = plt.subplots(1, numplots, sharey=True, squeeze=False, figsize=(6*numplots, 5))
      if args.ev == True:
        axs[0,0].set_ylabel("Energy [eV]")
      else:
        axs[0,0].set_ylabel("Energy [Ry]")

      for i in range(numplots):
        npoints, name, point, fermi = read_header(args.file)
        table = read_data(args.file)
        if(len(args.project) > 0):
          weights = read_data(args.file.replace("bandstructure","weights"))

        axs[0,i].set_title(titles[i])
        axs[0,i].set_xlim([point[0],point[npoints-1]])

        axs[0,i].set_xticks(point)
        axs[0,i].set_xticklabels(name)
        for j in point:
          axs[0,i].axvline(x=j, color='k', linewidth=0.5)

        # Ploting the Fermi level or a line at y=0.0
        if args.noef:
          if (fermi == None): # susceptibility
            axs[0,i].axhline(y=0.0, xmin=point[0], xmax=point[npoints-1], color='k', linestyle='-', linewidth=0.5)
          else: # band structure
            axs[0,i].axhline(y=0.0, xmin=point[0], xmax=point[npoints-1], color='k', linestyle='--', linewidth=0.5)
            if(i==0):
              axs[0,i].axhline(y=fermi, xmin=point[0], xmax=point[npoints-1], color='b', linestyle='--', linewidth=0.5)
            else:
              axs[0,i].axhline(y=-fermi, xmin=point[0], xmax=point[npoints-1], color='b', linestyle='--', linewidth=0.5)

        # Plotting the results
        if (fermi == None): # susceptibility
          axs[0,i].plot(table[:,1],-1.0/table[:,2])
        else: # band structure
          if(i==0):
            axs[0,i].plot(table[:,0],table[:,(table.shape[1]-1)/2+1:], color='r', linewidth=1.0, linestyle='-')
          else:
            axs[0,i].plot(table[:,0],table[:,1:(table.shape[1]-1)/2], color='k', linewidth=1.0, linestyle='-')

    else: # args.separate
      #legacy code
      numplots = 1
      titles = [r"Band Structure"]
      fig, axs = plt.subplots(1, numplots, sharey=True, squeeze=False, figsize=(6*numplots, 5))

      # Setting the units
      if args.ev == True:
        axs[0,0].set_ylabel("Energy [eV]")
      else:
        axs[0,0].set_ylabel("Energy [Ry]")

      for i in range(numplots):
        # Reads and stores the high symmetry points
        npoints, name, point, fermi, dimbs = read_header(args.file)
        # Reads the data points
        table = read_data(args.file)
        # Reads the weigths for projected band structure
        if(len(args.project) > 0):
          weights = read_data(args.file.replace("bandstructure","weights"))
        # Set custom title or default title
        if args.title != "":
          axs[0,i].set_title(args.title)
        else:
          axs[0,i].set_title(titles[i])

        axs[0,i].set_xlim([point[0],point[npoints-1]])
        # Used to zoom in a symmetrical region around zero
        if args.zoom != 0.0:
          axs[0,i].set_ylim(-args.zoom,args.zoom)

        axs[0,i].set_xticks(point)
        axs[0,i].set_xticklabels(name)

        # High symmetry points vertical lines
        for j in point:
          axs[0,i].axvline(x=j, color='k', linewidth=linewidth)

        # Ploting the Fermi level or a line at y=0.0
        if args.noef:
          if (fermi == None): # susceptibility
            axs[0,i].axhline(y=0.0, xmin=point[0], xmax=point[npoints-1], color='k', linestyle='-', linewidth=0.5)
          else: # band structure
            axs[0,i].axhline(y=0.0, color='k', linestyle='--',linewidth=0.5)

        # Wider lines for the frame
        for axis in ['top','bottom','left','right']:
          axs[0,i].spines[axis].set_linewidth(linewidth)

        # Plotting the results
        if (fermi == None): # susceptibility
          axs[0,i].plot(table[:,1],-1.0/table[:,2])
        else: # band structure
          # Print band structures as guides
          if args.guide == True:
            for j in range(1,dimbs+1):
              axs[0,i].scatter(table[:,0],table[:,j]*ry2ev,c='k',s=0.01,marker="o")

          for g, group in enumerate(eval(args.project)):
            if isinstance(group, int):
              orb = group
              for j in range(1,dimbs+1):
                weight = weights[:,orb+dimbs*(j-1)-1]
                axs[0,i].scatter(table[:,0],table[:,j]*ry2ev,c=[(0.07, 0.19, 0.40, a) for a in weight],s=10)
            else:
              for k in group:
                orb = k
                for j in range(1,dimbs+1):
                  weight = weights[:,orb+dimbs*(j-1)-1]
                  axs[0,i].scatter(table[:,0],table[:,j]*ry2ev,c=[(colors[g][0],colors[g][1],colors[g][2], a) for a in weight],s=8)

  else: # args.superconductivity:
    numplots = 1
    titles = [r"Band Structure"]
    fig, axs = plt.subplots(1, numplots, sharey=True, squeeze=False, figsize=(6*numplots, 5))
    # Setting the units
    if args.ev == True:
      axs[0,0].set_ylabel(r"E - E$_\mathrm{F}$ [eV]")
    else:
      axs[0,0].set_ylabel(r"E - E$_\mathrm{F}$ [Ry]")

    for i in range(numplots):
      # Reads and stores the high symmetry points
      npoints, name, point, fermi, dimbs = read_header(args.file)
      # Reads the data points
      table = read_data(args.file)
      # Reads the weigths for projected band structure
      if(len(args.project) > 0):
        weights = read_data(args.file.replace("bandstructure","weights"))
      # Set custom title or default title
      if args.title != "":
        axs[0,i].set_title(args.title)
      else:
        axs[0,i].set_title(titles[i])

      axs[0,i].set_xlim([point[0],point[npoints-1]])
      # Used to zoom in a symmetrical region around Ef
      if args.zoom != 0.0:
        # print(fermi*ry2ev-args.zoom,fermi*ry2ev+args.zoom)
        axs[0,i].set_ylim(-args.zoom,args.zoom)

      axs[0,i].set_xticks(point)
      axs[0,i].set_xticklabels(name)

      # High symmetry points vertical lines
      for j in point:
        axs[0,i].axvline(x=j, color='k', linewidth=linewidth)

      # Ploting the Fermi level or a line at y=0.0
      if args.noef:
        if (fermi == None): # susceptibility
          axs[0,i].axhline(y=0.0, xmin=point[0], xmax=point[npoints-1], color='k', linestyle='-', linewidth=0.5)
        else: # band structure
          axs[0,i].axhline(y=0.0, color='k', linestyle='--',linewidth=0.5)

      # Wider lines for the frame
      for axis in ['top','bottom','left','right']:
        axs[0,i].spines[axis].set_linewidth(linewidth)

      # Plotting the results
      if (fermi == None): # susceptibility
        axs[0,i].plot(table[:,1],-1.0/table[:,2])
      else: # band structure
        # Print band structures as guides
        if args.guide == True:
          for j in range(1,dimbs+1):
            axs[0,i].scatter(table[:,0],(table[:,j]-fermi)*ry2ev,c='k',s=0.01,marker="o")

        for g, group in enumerate(eval(args.project)):
          if isinstance(group, int):
            orb = group
            for j in range(1,dimbs+1):
              weight = weights[:,orb+dimbs*(j-1)-1]
              axs[0,i].scatter(table[:,0],(table[:,j]-fermi)*ry2ev,c=[(0.07, 0.19, 0.40, a) for a in weight],s=10)
          else:
            for k in group:
              orb = k
              for j in range(1,dimbs+1):
                weight = weights[:,orb+dimbs*(j-1)-1]
                axs[0,i].scatter(table[:,0],(table[:,j]-fermi)*ry2ev,c=[(colors[g][0],colors[g][1],colors[g][2], a) for a in weight],s=8)

  plt.tight_layout()
  if args.output != "":
    plt.savefig(args.output,dpi=200)
  plt.show()
