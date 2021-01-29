################################################################################
# Import the arguments and style
################################################################################
from head import *

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

  if(args.mev):
    labelx = r'$E-E_F$ [meV]'
    labely = r'LDOS [states/meV]'
  elif(args.ev):
    labelx = r'$E-E_F$ [eV]'
    labely = r'LDOS [states/eV]'
  else:
    labelx = r'$E-E_F$ [Ry]'
    labely = r'LDOS [states/Ry]'

  if args.total is False:
    numplots = len(args.files)
    if numplots > 2:
      sys.exit("Too many input files")
  else:
    numplots = 1

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


  # Getting the orbitals to plot...
  orbitals = eval(args.orbitals)

  # ...and their legends
  if args.legends != "[]":
    legends = [r"{}".format(string) for string in args.legends.split(",")]
    if len(legends) != len(orbitals):
      sys.exit("Incorrect number of legends. Please provide one for each orbital. len(legends) = " + str(len(legends)) + ", len(orbitals) = " + str(len(orbitals)))
    if numplots == 2:
      for i in range(len(orbitals)):
        legends.append(r"")


  fig, ax = plt.subplots(1, 1, figsize=(6, 1 + 1.5*numplots))

  plt.axhline(0.0, color='k', linewidth=0.75)

  data = []
  data_int = []
  if args.total is False:
    for j in range(numplots):
      data.append(read_data(args.files[j],orbitals=orbitals))
      data[j] = data[j][data[j][:,0].argsort()]
      if args.integrated is True:
        data_int.append(np.array([integrate.cumtrapz(data[j][:,k], data[j][:,0], initial=0) for k in range(1,len(data[j][0,:])) ]).transpose())
  else:
    data.append(read_data(args.files[0],orbitals=orbitals))
    for j in range(1,len(args.files)):
      data_temp = read_data(args.files[j],orbitals=orbitals)
      for i, line in enumerate(data[0]):
        data[0][i] = np.array([data_temp[i,0]] + [data[0][i,k] + data_temp[i,k] for k in range(1,len(data_temp[i,:]))])
    data[0] = data[0][data[0][:,0].argsort()]
    if args.integrated is True:
      data_int.append(np.array([integrate.cumtrapz(data[0][:,k], data[0][:,0], initial=0) for k in range(1,len(data[0][0,:])) ]).transpose())


  a1 = []
  for j in range(numplots):
    # data = read_data(args.files[j],orbitals=orbitals)
    # datau = np.loadtxt(filenameu)
    x = [(data[j][k,0]-fermi)*ry2ev for k in range(len(data[j][:,0]))]
    signal = (-1)**j
    for i in range(1,len(data[j][0,:])):
      y = [signal*data[j][k,i]/ry2ev for k in range(len(data[j][:,i]))]
      plt.plot(x, y, color=colors[i-1], label=legends[j*(len(data[j][0,:])-1)+(i-1)])
    # Plotting integrated DOS
    if args.integrated is True:
      signal = (-1)**j
      for i in range(0,len(data_int[j][0,:])):
        y = [signal*data_int[j][k,i]/ry2ev for k in range(len(data_int[j][:,i]))]
        plt.plot(x, y, ls='--', color=colors[i], label=(legends[j*(len(data_int[j][0,:]))+i]+" int" if legends[j*(len(data_int[j][0,:]))+i] != "" else "") )


    a1.append(max(data[j][:,1]/ry2ev))
  max_ldos = max(a1)

  if args.xlim != "":
    xlim = eval(args.xlim)
    # xlim = [-eval(args.xlim),eval(args.xlim)]
    plt.xlim(xlim)

  if args.domain != "":
    plt.xlim(eval(args.domain))

  if args.ylim != "":
    ylim = [(-eval(args.ylim) if numplots > 1 else 0.0),eval(args.ylim)]
  else:
    ylim = [(-1.1*abs(max_ldos) if numplots > 1 else 0.0),1.1*abs(max_ldos)]

  plt.ylim(ylim)


  plt.xlabel(labelx)
  plt.ylabel(labely)
  if args.title == "":
    plt.title('LDOS')
  else:
    plt.title(args.title)

  plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

  if args.noef:
    plt.axvline(0.0, color='k', linestyle='--', linewidth=0.75)
  

  # fig = plt.gcf()
  # fig.set_size_inches(3., 2.2)
  plt.tight_layout()
  if args.output == "":
    plt.show()
  else:
    plt.savefig(args.output,dpi=200)
    if args.show:
      plt.show()
