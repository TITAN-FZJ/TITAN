import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib as mpl                      # Plotting library
from matplotlib import rc                     # Improve math fonts
from mpl_toolkits.axes_grid1 import AxesGrid  # Grid plotting
import matplotlib.colors as colors            # Color selection and manipulation
import scipy.interpolate                      # Interpolation library
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
mpl.rcParams['font.family']      = 'Arial'
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


ry2ev = 1.0
# ry2ev = 13.6057 # Conversion of energy units

# plt.xkcd()



################################################################################
# Get the header from the file
################################################################################
def read_header(file):
  with open(file, "r") as f:
    Ef_line = f.readline().split()
    fermi = None
    if "Ef" in Ef_line[1]:
      fermi = float(Ef_line[2])
      # print fermi
  return fermi

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


if __name__ == "__main__":

  if(ry2ev != 1.0):
    labelx = r'$E-E_F$ [eV]'
    labely = r'LDOS [states/eV]'
  else:
    labelx = r'$E-E_F$ [Ry]'
    labely = r'LDOS [states/Ry]'

  filenameu = sys.argv[1]
  filenamed = sys.argv[2]

  # Getting fermi energy from ldos up file
  fermi = read_header(filenameu)

  # colors = brewer2mpl.get_map('Set1', 'qualitative', 5).mpl_colors
  colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#a65628']
  # colors = np.array(['k', 'g', 'r', 'b']) # The colors you wish to cycle through
  legends = np.array(['Total', 's', 'p', 'd'])

  plt.axhline(0.0, color='k', linewidth=0.75)

  datau = read_data(filenameu)
  # datau = np.loadtxt(filenameu)
  datau = datau[datau[:,0].argsort()]
  x = datau[:,0]
  for i in range(1,len(datau[0,:])):
    plt.plot((x-fermi)*ry2ev,datau[:,i]/ry2ev, color=colors[i-1], label=legends[i-1])

  datad = read_data(filenamed)
  # datad = np.loadtxt(filenamed)
  datad = datad[datad[:,0].argsort()]
  x = datad[:,0]
  for i in range(1,len(datad[0,:])):
    plt.plot((x-fermi)*ry2ev,-datad[:,i]/ry2ev, color=colors[i-1])

  a1 = max(datau[:,1]/ry2ev)
  a2 = max(datad[:,1]/ry2ev)
  max_ldos = max(a1,a2)
  ylim = [-1.1*abs(max_ldos),1.1*abs(max_ldos)]
  plt.ylim(ylim)


  plt.xlabel(labelx)
  plt.ylabel(labely)
  plt.title('LDOS')
  plt.legend()

  # if(Ef != 0.0):
  plt.axvline(0.0, color='k', linestyle='--', linewidth=0.75)

  # fig = plt.gcf()
  # fig.set_size_inches(3., 2.2)
  plt.savefig('LDOS.pdf')
  # plt.savefig(filename,transparent=True,dpi=2540, bbox_inches='tight')

  plt.show()
