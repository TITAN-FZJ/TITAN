################################################################################
# Import the arguments and style
################################################################################
from head import *
import matplotlib.colors as colors            # Color selection and manipulation
import scipy.interpolate                      # Interpolation library
from mpl_toolkits.axes_grid1 import make_axes_locatable

col = 3
# kx,ky,kz,charge,spin x,spin y,spin z,morb x,morb y,morb z     

################################################################################
# Get the header from the file
################################################################################
def read_header(file):
  with open(file, "r") as f:
    header = [s for s in f.readline().replace(" ", "").replace("#", "").replace("\n", "").split(',')]
  return header

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
  nplots = len(sys.argv)-1
  
  # Create a grid with 1 row and 'numplots' column
  fig, ax = plt.subplots(1, nplots, figsize=(5*nplots,4), squeeze=False)
  plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.85, hspace=0.25,
                      wspace=0.6)

  for i in range(nplots):
    file = sys.argv[i+1]
    header = read_header(file)

    clabel = header[col]

    table = read_data(file)

    # Sorting values
    table.sort()

    # Separating the columns
    x = np.array([v[0] for v in table])
    y = np.array([v[1] for v in table])
    z = np.array([v[col] for v in table])

    im = ax[0,i].tricontourf(x, y, z, 100, cmap='bwr',vmin=z.min(), vmax=z.max(), norm=MidpointNormalize(vmin=z.min(), vmax=z.max(), midpoint=0))

    ax[0,i].set_aspect(np.diff(ax[0,i].get_xlim())/np.diff(ax[0,i].get_ylim()))
    for axis in ['top','bottom','left','right']:
      ax[0,i].spines[axis].set_linewidth(1.5)

    # Add colorbar
    cax = fig.add_axes([ax[0,i].get_position().x1+0.01,ax[0,i].get_position().y0,0.02,1.0*ax[0,i].get_position().height])
    plt.colorbar(im, cax=cax).outline.set_linewidth(1.5) # Similar to fig.colorbar(im, cax = cax)
    cax.set_ylabel(clabel)

  plt.show()
