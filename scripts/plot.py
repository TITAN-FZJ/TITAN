################################################################################
# Routine to plot {x,y} from a data file with many columns 
# (x is obtained from column 0)
# @author Filipe Guimarães
################################################################################
import numpy 	as np              # Numerical library
import sys                         # System library (to read arguments from command line)
import matplotlib.pyplot as plt    # Plotting library
import matplotlib as mpl           # Plotting library
from matplotlib import rc          # Improve math fonts

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
# rc('font',**{'family':'serif','serif':['Palatino']})
mpl.rcParams['text.usetex']         = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
# Default fonts
mpl.rcParams['font.size']        = 11
mpl.rcParams['font.family']      = 'Arial'
mpl.rcParams['figure.titlesize'] = 'large'
mpl.rcParams['axes.titlepad']    = 10
mpl.rcParams['lines.linewidth']  = 2
#Legends:
mpl.rcParams['legend.fontsize']     = 'medium'
mpl.rcParams['legend.fancybox']     = False
# mpl.rcParams['legend.loc']          = 'upper left'`
mpl.rcParams["font.weight"]         = "bold"
mpl.rcParams['legend.edgecolor']    = 'inherit'
mpl.rcParams["axes.labelweight"]    = "bold"
mpl.rcParams['legend.framealpha']   = None
mpl.rcParams['legend.handlelength'] = 2

# Getting filename from 1st argument of command line
filename = sys.argv[1]
# Saving information in matrix form
data=np.loadtxt(filename)		# A matrix containing all of the data in your file.
# Obtaining the number of columns in the matrix
ncol=np.size(data,axis=1)

# Defining the colors for the curves (see http://colorbrewer2.org/)
colors = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928']
#['#e41a1c','#377eb8','#4daf4a','#984ea3','#a65628']
# colors = np.array(['b', 'g', 'r', 'c', 'm', 'y', 'k'])	# The colors you wish to cycle through

#Defining legends for the curves
legends = ['curve 1', 'curve 2', 'curve 3', 'curve 4', 'curve 5', 'curve 6', 'curve 7', 'curve 8', 'curve 9']

# Loop over columns i from 1 to ncol
for i in np.arange(1,ncol):
	# Plot column i against column 0 (first column)
  plt.plot(data[:,0],data[:,i], color=colors[i-1], label=legends[i-1])

# Labels and title
plt.xlabel('x-axis label')
plt.ylabel('y-axis label')
plt.title('Title of the graph')
# Add legend
plt.legend()

plt.show()