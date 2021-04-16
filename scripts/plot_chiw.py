################################################################################
# Import the arguments and style
################################################################################
from head import *

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
  
  if(args.mev):
    labelx = r'$E-E_F$ [meV]'
    labely = r'$-\operatorname{Im}\chi^{+-}$ [states/meV]'
  elif(args.ev):
    labelx = r'$E-E_F$ [eV]'
    labely = r'$-\operatorname{Im}\chi^{+-}$ [states/eV]'
  else:
    labelx = r'$E-E_F$ [Ry]'
    labely = r'$-\operatorname{Im}\chi^{+-}$ [states/Ry]'

  numplots = len(args.files)
  fig, ax = plt.subplots(1, 1, figsize=(6, 5))

  # ...and their legends
  legends = []
  if args.legends != "[]":
    legends = [r"{}".format(string) for string in args.legends.split(",")]
    if len(legends) != len(numplots):
      sys.exit("Incorrect number of legends. Please provide one for each orbital. len(legends) = " + str(len(legends)) + ", len(orbitals) = " + str(len(orbitals)))

  # Getting data from input file(s)
  data = []
  for i in range(numplots):
    data.append(read_data(args.files[i]))

  # Defining the colors for the curves (see http://colorbrewer2.org/)
  colors = colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#a65628', '#b98600', '#000000']
  #['#e41a1c','#377eb8','#4daf4a','#984ea3','#a65628']
  # colors = np.array(['b', 'g', 'r', 'c', 'm', 'y', 'k'])  # The colors you wish to cycle through

  #Defining legends for the curves
  # legends = ['curve 1', 'curve 2', 'curve 3', 'curve 4', 'curve 5', 'curve 6', 'curve 7', 'curve 8', 'curve 9']

  # Loop over columns i from 1 to ncol
  for i in range(numplots):
    # Plot column i against column 0 (first column)
    plt.plot(data[i][:,0],-data[i][:,3], color=colors[i], label=(legends[i] if args.legends != "[]" else None) )

  # Labels and title
  plt.xlabel(labelx)
  plt.ylabel(labely)
  plt.title('Susceptibility')

  # Changing limits
  if args.xlim != "":
    xlim = eval(args.xlim)
    plt.xlim(xlim)

  # Add legend
  if args.legends != "[]":
    plt.legend()

  if args.output == "":
    plt.show()
  else:
    plt.savefig(args.output,dpi=200)
    if args.show:
      plt.show()
