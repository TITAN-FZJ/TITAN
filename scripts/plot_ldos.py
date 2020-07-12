import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib as mpl                      # Plotting library
from matplotlib import rc                     # Improve math fonts
from mpl_toolkits.axes_grid1 import AxesGrid  # Grid plotting
import matplotlib.colors as colors            # Color selection and manipulation
import scipy.interpolate                      # Interpolation library
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

parser = argparse.ArgumentParser(description="Parse bool")
parser.add_argument("files", nargs="+", help="File(s) to plot")
parser.add_argument("--superconductivity", default=False, action="store_true" , help="Flag to plot the bands as for superconductors")
parser.add_argument("--mev", default=False, action="store_true" , help="Plot superconductor bands in the same plot")
args = parser.parse_args()

if args.mev:
    ry2ev = 13.6057*1000 # Conversion of energy units


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
    labelx = r'$E-E_F$ [meV]'
    labely = r'LDOS [states/meV]'
  else:
    labelx = r'$E-E_F$ [Ry]'
    labely = r'LDOS [states/Ry]'


  numplots = len(args.files)
  if numplots > 2:
    sys.exit("Too many input files")

  # Getting fermi energy from ldos up file
  fermi = read_header(args.files[0])

  if args.superconductivity:
      # colors = brewer2mpl.get_map('Set1', 'qualitative', 5).mpl_colors
      colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#a65628', '#b98600', '#000000']
      # colors = np.array(['k', 'g', 'r', 'b']) # The colors you wish to cycle through
      legends = np.array(['Total', r'$u_s$', r'$u_p$', r'$u_d$',r'$v_s$', r'$v_p$', r'$v_d$','', '', '', '','', '', ''])
  else:
      # colors = brewer2mpl.get_map('Set1', 'qualitative', 5).mpl_colors
      colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#a65628']
      # colors = np.array(['k', 'g', 'r', 'b']) # The colors you wish to cycle through
      legends = np.array(['Total', 's', 'p', 'd', '', '', '', ''])

  fig, ax = plt.subplots(1, 1, figsize=(6, 1 + 1.5*numplots))

  plt.axhline(0.0, color='k', linewidth=0.75)

  a1 = []
  for j in range(numplots):
    data = read_data(args.files[j])
    # datau = np.loadtxt(filenameu)
    data = data[data[:,0].argsort()]
    x = [(data[k,0]-fermi)*ry2ev for k in range(len(data[:,0]))]
    signal = (-1)**j
    for i in range(1,len(data[0,:])):
      y = [signal*data[k,i]/ry2ev for k in range(len(data[:,i]))]
      plt.plot(x, y, color=colors[i-1], label=legends[j*(len(data[0,:])-1)+(i-1)])

    a1.append(max(data[:,1]/ry2ev))
  max_ldos = max(a1)

  ylim = [(-1.1*abs(max_ldos) if numplots > 1 else 0.0),1.1*abs(max_ldos)]
  plt.ylim(ylim)


  plt.xlabel(labelx)
  plt.ylabel(labely)
  plt.title('LDOS')
  plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

  # if(Ef != 0.0):
  plt.axvline(0.0, color='k', linestyle='--', linewidth=0.75)

  # fig = plt.gcf()
  # fig.set_size_inches(3., 2.2)
  plt.tight_layout()
  plt.savefig('LDOS.pdf')
  # plt.savefig(filename,transparent=True,dpi=2540, bbox_inches='tight')

  plt.show()
