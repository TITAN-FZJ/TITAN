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
ry2mev = 1.0


ry2mev = 13605.7 # Conversion of energy units
fact = 10         # Interpolation factor
# ry2mev = 1.0

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
  return npoints, name_out, point, fermi

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
        data.append( [float(x) for x in line.split()] )
  # ndata = np.array(data)
  return data

################################################################################
# Interpolate the values
################################################################################
def interpolate(values):
  # Sorting values
  values.sort()

  # Separating the values for the fit
  x = np.array([v[0] for v in values])
  y = np.array([v[1] for v in values])

  # Getting number of points in x and y: nx repetitions of y[0]
  nx = sum(y == y[0])
  ny = len(x)/nx

  # Removing repetitions
  x = x[::ny]
  y = y[:ny]

  # Quantity to plot in the Colormap
  z = np.array([v[3] for v in values])

  z.shape = (nx, ny)

  # Interpolation for data on rectangular grids using a bivariate cubic spline
  f = scipy.interpolate.RectBivariateSpline(x, y, z)

  # Getting new interpolated values (tenfold points on each axis, adjust this depending on the input file)
  x = np.linspace(x.min(), x.max(), fact*nx)
  y = np.linspace(y.min(), y.max(), fact*ny)
  z = f(x, y)

  return x*ry2mev, y, -z/ry2mev

################################################################################
# set the colormap and centre the colorbar
# Obtained from:
# https://stackoverflow.com/a/50003503/3142385
################################################################################
class MidpointNormalize(colors.Normalize):
  def __init__(self, vmin, vmax, midpoint=0, clip=False):
    self.midpoint = midpoint
    colors.Normalize.__init__(self, vmin, vmax, clip)

  def __call__(self, value, clip=None):
    normalized_min = max(0.0, 1.0 / 2.0 * (1.0 - abs((float(self.midpoint) - float(self.vmin)) / (float(self.midpoint) - float(self.vmax)))))
    normalized_max = min(1.0, 1.0 / 2.0 * (1.0 + abs((float(self.vmax) - float(self.midpoint)) / (float(self.midpoint) - float(self.vmin)))))
    normalized_mid = 0.5
    x, y = [self.vmin, self.midpoint, self.vmax], [normalized_min, normalized_mid, normalized_max]
    return np.ma.masked_array(np.interp(value, x, y))


################################################################################
# Main program
################################################################################
if __name__ == "__main__":
  numplots = len(sys.argv)-1
  titles = [r'', r'']
  clabel = [r'-Im$\chi^{+-}(\omega,\mathbf{q})$']
  # titles = [r"$\#_{k}=100$M, $\eta=5\times10^{-4}$, 2nn", r"$\#_{k}=100$M, $\eta=5\times10^{-4}$, 3nn"]
  # titles = [r"$\#_{k}=10$M", r"$\#_{k}=100$k"])
  # titles = [r"$\eta=5\times10^{-3}$", r"$\eta=5\times10^{-4}$"])

  fig = plt.figure(figsize=(6*numplots, 5))

  # Create a grid with 1 row and 'numplots' column
  grid = AxesGrid(fig, 111,
                     nrows_ncols=(1, numplots),
                     axes_pad=0.2,
                     share_all=True,
                     cbar_mode='single',
                     cbar_location='right',
                     cbar_pad=0.2 )

  for i in range(numplots):
    npoints, name, point, fermi = read_header(sys.argv[i+1])
    table = read_data(sys.argv[i+1])
    y, x, z = interpolate(table)
    grid[i].set_title(titles[i])
    if i == 0:
      grid[i].set_ylabel("Energy [meV]")
    grid[i].set_xlabel("Wave vector")

    grid[i].set_xlim([point[0],point[npoints-1]])
    # grid[i].set_ylim(y.min(),y.max())
    grid[i].set_ylim(y.min(),330.0)
    # grid[i].set_ylim(1.0/(table[:,1:].min()),1.0/(table[:,1:].max()))

    grid[i].set_xticks(point)
    grid[i].set_xticklabels(name)
    for j in point:
      grid[i].axvline(x=j, color='k', linewidth=0.5)

    # Ploting the Fermi level or a line at y=0.0
    if not (fermi == None):
      grid[i].axhline(y=fermi, xmin=point[0], xmax=point[npoints-1], color='k', linestyle='--')
    else:
      grid[i].axhline(y=0.0, xmin=point[0], xmax=point[npoints-1], color='k', linestyle='-', linewidth=0.5)

    # for j in enumerate(table[:,1]):
    #   if j[1] == point[4]:
    #     print 'Value at ',name[4],': ',-1.0/table[j[0],2]

    # Plotting the results
    # im = grid[i].pcolormesh(x, y, z, cmap='jet', vmin=z.min(), vmax=z.max(), norm=MidpointNormalize(vmin=z.min(), vmax=z.max(), midpoint=0), rasterized=True)
    im = grid[i].pcolormesh(x, y, z,cmap='jet', vmin=0.0, vmax=z.max()/50.0)
    # im = grid[i].contourf(x, y, z, 800,cmap='jet', vmin=0.0, vmax=z.max()/5.0)
    #'seismic''RdBu_r'
    # grid[i].plot(table[:,1],-1.0/table[:,2])
    grid[i].set_aspect(np.diff(grid[i].get_xlim())/np.diff(grid[i].get_ylim()))
    for axis in ['top','bottom','left','right']:
      grid[i].spines[axis].set_linewidth(1.5)

    cb = grid.cbar_axes[0].colorbar(im)
    grid.cbar_axes[0].tick_params(width=2,colors='black')
    grid.cbar_axes[0].set_ylabel(clabel[i])
    # grid[i].cbar_axes[0].ticklabel_format(style='sci', scilimits=(0,0))
    grid.cbar_axes[0].yaxis.set_offset_position('left')
    # grid.cbar_axes[0].set_title(prm.cbar_titles_SOT[i] if 'SOT' in prm.output_prefix else prm.cbar_titles_Beff[i], size=12) # Title of the colorbar
    cb.solids.set_rasterized(True)
    for axis in ['top','bottom','left','right']:
      grid.cbar_axes[0].spines[axis].set_linewidth(1.5)



  # plt.tight_layout()
  plt.show()
