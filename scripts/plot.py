import numpy 	as np              # Numerical library
import sys                       # System library (to read arguments from command line)
from matplotlib.pyplot import *  # Plotting library
import matplotlib.pyplot as plt  # Plotting library
# Default fonts
rcParams['font.size']        = 12
rcParams['font.family']      = 'Arial'
rcParams['figure.titlesize'] = 'large'
rcParams['lines.linewidth']  = 2
#Legends:
rcParams['legend.fontsize']  = 'medium'
rcParams['legend.fancybox'] = False
rcParams['legend.loc'] = 'upper right'
rcParams['legend.framealpha'] = None
rcParams['legend.edgecolor'] = 'inherit'
rcParams['legend.handlelength'] = 2
from matplotlib import rc        # Improve math fonts
rc('mathtext', default='regular')

# Getting filename from 1st argument of command line
filename = sys.argv[1]
# Saving information in matrix form
data=np.loadtxt(filename)		# A matrix containing all of the data in your file.
# Obtaining the number of columns in the matrix
ncol=np.size(data,axis=1)

# Defining the colors for the curves (see http://colorbrewer2.org/)
colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#a65628']
# colors = np.array(['b', 'g', 'r', 'c', 'm', 'y', 'k'])	# The colors you wish to cycle through

#Defining legends for the curves
legends = ['curve 1', 'curve 2', 'curve 3', 'curve 4', 'curve 5']

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