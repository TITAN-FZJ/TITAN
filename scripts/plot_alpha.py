################################################################################
# Import the arguments and style
################################################################################
from head import *
import re

################################################################################
# Get the data from the file and save it into a matrix
################################################################################
def read_data(filename,orbitals=[]):
  # open a file using with statement
  data = []
  with open(filename,'r') as file:
    for line in file:
      # check if the current line
      # starts with "#"
      if not (line.startswith("#") or line.startswith(" #")):
        data_temp = [float(x) for x in line.split()]
        # print(data_temp)
        if len(orbitals) == 0:
          data.append(data_temp)
        else:
          data.append([data_temp[0]] + [data_temp[x] for x in orbitals])
  ndata = np.array(data)
  return ndata


if __name__ == "__main__":

  labelx = r'Broadening $\eta$'
  if(args.mev) or (args.mevlabel):
    labelx = labelx + r' [meV]'
  elif(args.ev) or (args.evlabel):
    labelx = labelx + r' [eV]'
  else:
    labelx = labelx + r' [Ry]'
  labely = r'$\alpha$'


  numcurves = len(args.files)

  # colors = brewer2mpl.get_map('Set1', 'qualitative', 5).mpl_colors
  # colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#a65628']
  colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999']
  # colors = np.array(['k', 'g', 'r', 'b']) # The colors you wish to cycle through
  # legends = np.array(['Total', 's', 'p', 'd', '', '', '', ''])


  fig, ax = plt.subplots(1, 1, figsize=(10, 6))
  ax.set_xscale('log')
  ax.set_yscale('log')

  sites = {}
  for file in args.files:
    # Getting eta from filename
    eta = float(re.findall(r'eta=(.*).dat',file)[0])

    # Reading data from file
    data = read_data(file)
    # Initializing dictionary for each site
    for i in range(len(data[:,0])):
      sites.setdefault(int(data[i,0]), {})
      sites[int(data[i,0])].setdefault('eta', [])
      if int(data[i,1]) == 1:
        sites[int(data[i,0])]['eta'].append(eta)

        sites[int(data[i,0])].setdefault('alphaxx', [])
        sites[int(data[i,0])]['alphaxx'].append(data[i,2])
      elif int(data[i,1]) == 2:
        sites[int(data[i,0])].setdefault('alphayy', [])
        sites[int(data[i,0])]['alphayy'].append(data[i,3])

  xmin = 1.0
  xmax = 0.0
  ymin = 1.0
  ymax = 0.0
  # Plotting alpha_xx
  for site in sites:
    x = sites[site]['eta']
    y = sites[site]['alphaxx']
    if all(v < 1.e-9 for v in y): 
      continue

    eta, alpha = (list(t) for t in zip(*sorted(zip(x, y))))
    ax.plot(eta, alpha, label=r"$\alpha_{xx}^{"+str(site)+"}$", color=colors[site-1], marker='o', ms='5', lw=1)
    xmin = min([xmin]+x)
    xmax = max([xmax]+x)
    ymin = min([ymin]+y)
    ymax = max([ymax]+y)

  # Plotting alpha_yy
  for site in sites:
    x = sites[site]['eta']
    y = sites[site]['alphayy']
    if all(v < 1.e-9 for v in y): 
      continue
    eta, alpha = (list(t) for t in zip(*sorted(zip(x, y))))
    ax.plot(eta, alpha, ls='--', label=r"$\alpha_{yy}^{"+str(site)+"}$", color=colors[site-1], marker='o', ms='5', lw=1)
    xmin = min([xmin]+x)
    xmax = max([xmax]+x)
    ymin = min([ymin]+y)
    ymax = max([ymax]+y)

  # Setting limits
  ax.set_xlim([0.5*xmin,1.1*xmax])
  ax.set_ylim([0.5*ymin,1.1*ymax])
  # Setting labels
  ax.set_xlabel(labelx)
  ax.set_ylabel(labely)
  ax.legend(loc='lower right',ncol=2,fontsize=10)

  plt.tight_layout()
  if args.output == "":
    plt.show()
  else:
    plt.savefig(args.output,dpi=200)
    if args.show:
      plt.show()
