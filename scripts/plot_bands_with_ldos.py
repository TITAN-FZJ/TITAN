################################################################################
# This script with all its options is to be run like this
# ipython <script_location> <band_file> <ldosu_file> <ldosd_file> -- --superconductivity --mev --title="desired title"
# the arguments starting with "--" are optional and can be ommited
################################################################################

################################################################################
# Import the arguments and style
################################################################################
from head import *

if args.superconductivity:
    # colors = brewer2mpl.get_map('Set1', 'qualitative', 5).mpl_colors
    colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#a65628', '#b98600', '#000000']
    # colors = np.array(['k', 'g', 'r', 'b']) # The colors you wish to cycle through
    legends = np.array(['Total', r'$u_s$', r'$u_p$', r'$u_d$',r'$v_s$', r'$v_p$', r'$v_d$'])
else:
    # colors = brewer2mpl.get_map('Set1', 'qualitative', 5).mpl_colors
    colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#a65628']
    # colors = np.array(['k', 'g', 'r', 'b']) # The colors you wish to cycle through
    legends = np.array(['Total', 's', 'p', 'd'])


def read_header(file):
  with open(file, "r") as f:
    count = [s for s in f.readline().split()]
    npoints = int(count[1])
    # print npoints
    name = []
    name_in = np.empty(npoints, dtype=str)
    point = np.empty(npoints, dtype=float)
    for i in range(npoints):
      name_in[i], point[i] = f.readline().split()[1:]
      if name_in[i] == "G":
        name.append(r"$\Gamma$")
      else:
        name.append(name_in[i])
      # print name[i], point[i]

    Ef_line = f.readline().split()
    fermi = None
    if "Ef" in Ef_line[1]:
      fermi = float(Ef_line[2])
      # print fermi
  return npoints, name, point, fermi

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

def read_header_ldos(file):
  with open(file, "r") as f:
    Ef_line = f.readline().split()
    fermi = None
    if "Ef" in Ef_line[1]:
      fermi = float(Ef_line[2])
      # print fermi
  return fermi

def read_data_ldos(filename):
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
  input_error = False
  if (len(args.files)==1):
    if (args.ldosu != ""):
      ldosu = list(args.ldosu)
      if (args.ldosu != ""):
        ldosd = list(args.ldosd)
        mag = True
        nplots = 3
      else:
        mag = False
        nplots = 2
    else:
      input_error = True

  elif (len(args.files)<=1 or len(args.files) >=4):
    input_error = True
  else:
    ldosu = [args.files[1]] #args.fileu

    if len(args.files)==3 :
      ldosd = [args.files[2]] #args.filed
      mag = True
    else:
      mag = False
    nplots = len(args.files)


  if input_error:
    print(f"Incorrect number of arguments: {len(args.files)}.")
    print(f"At least one band structure and one LDOS file and at maximum")
    print(f"one band structure and two LDOS files shoud be given.")
    print(f"(more than one LDOS files can be given using arguments --ldosu --ldosd)")
    exit(1)

  # Preparing labels
  if args.centeref:
    labely = r'$E-E_F$'
  else:
    labely = r'$E$'

  if(args.mev) or (args.mevlabel):
    labelx = r'$\rho(E)$ [meV$^{-1}$]'
    labely = labely + r' [meV]'
  elif(args.ev) or (args.evlabel):
    labelx = r'$\rho(E)$ [eV$^{-1}$]'
    labely = labely + r' [eV]'
  else:
    labelx = r'$\rho(E)$ [Ry$^{-1}$]'
    labely = labely + r' [Ry]'

  # Getting band structure file
  bsstruct = args.files[0]#args.fileband

  # Preparing axes
  fig, ax = plt.subplots(1,nplots, sharey=True, gridspec_kw = {'width_ratios':([1,4,1] if mag else [1,4])})
  fig.subplots_adjust(left=0.1,wspace=0.15,top=0.95, bottom=0.15, right=0.9)

  if args.title != "":
      fig.suptitle(args.title,fontsize=24)

  # Configuting first axes and labels
  ax[0].tick_params(axis='y', direction='in', left=True, right=True)
  ax[0].set_ylabel(labely, fontsize=14)
  ax[0].set_xlabel(labelx, fontsize=14)

  # Getting band structure values and fermi value from LDOS header
  npoints, name, point, fermi = read_header(bsstruct)
  table = read_data(bsstruct)
  fermi_ldos = read_header_ldos(ldosu[0])
  if args.centeref:
    shift = fermi_ldos
    fermi_ldos = 0.0
  else:
    shift = 0.0

  # Getting LDOS values and summing up if more than one file is given
  data = []
  data.append(read_data_ldos(ldosu[0]))
  for j in range(1,len(ldosu)):
    data_temp = read_data_ldos(ldosu[j])
    for i, line in enumerate(data[0]):
      data[0][i] = np.array([data_temp[i,0]] + [data[0][i,k] + data_temp[i,k] for k in range(1,len(data_temp[i,:]))])
  ndatau = data[0][data[0][:,0].argsort()]

  # If up and down LDOS files are given, get values of the second one and setup extra graph
  if mag:
    ax[2].tick_params(axis='y', direction='in', left=True, right=True)
    ax[2].set_xlabel(labelx, fontsize=14)
    data = []
    data.append(read_data_ldos(ldosd[0]))
    for j in range(1,len(ldosd)):
      data_temp = read_data_ldos(ldosd[j])
      for i, line in enumerate(data[0]):
        data[0][i] = np.array([data_temp[i,0]] + [data[0][i,k] + data_temp[i,k] for k in range(1,len(data_temp[i,:]))])
    ndatad = data[0][data[0][:,0].argsort()]

  # Preparing band structure plot
  ax[1].set_xticks(point)
  ax[1].set_xticklabels(name)
  for i in point:
    ax[1].axvline(x=i, color='k', linewidth=0.75)
  if args.noef:
    ax[1].axhline(y=fermi_ldos, color='k', linestyle='--', linewidth=0.75)

  if args.gap !=0.0:
    for i in range(len(ax)):
      ax[i].axhline(y=-args.gap, color='r', linestyle='-.', linewidth=0.75)
      ax[i].axhline(y= args.gap, color='r', linestyle='-.', linewidth=0.75)

  # Plotting band structure
  if args.superconductivity:
    ax[1].plot(table[:,0],table[:,(table.shape[1]-1)/2+1:]*ry2ev-shift*ry2ev, color='r', linewidth=1.0, linestyle='-')
    ax[1].plot(table[:,0],table[:,1:(table.shape[1]-1)/2+1]*ry2ev-shift*ry2ev, color='k', linewidth=1.0, linestyle='-')
  else:
    ax[1].plot(table[:,0],table[:,1:]*ry2ev-shift*ry2ev, color='k', linewidth=1.0, linestyle='-')
  # ymax = np.max(table[:,1:])

  ax[1].tick_params(axis='y', direction='in', left=True, right=True)

  ax[1].set_xlim([point[0],point[npoints-1]])

  a1 = max(ndatau[:,1]/ry2ev)

  # Majority spin LDOS
  x = ndatau[:,0]
  if args.onlyS:
    ax[0].plot(-ndatau[:,2]/ry2ev,(x-shift)*ry2ev,linestyle='-', color=colors[2-1], label=legends[2-1],marker=1,markersize=1)
    if args.superconductivity:
      ax[0].plot(-ndatau[:,5]/ry2ev,(x-shift)*ry2ev,linestyle='-', color=colors[5-1], label=legends[5-1],marker=1,markersize=1)
      ax[0].plot(-ndatau[:,5]/ry2ev - ndatau[:,2]/ry2ev,(x-shift)*ry2ev,linestyle='-', color=colors[1-1], label="sum",marker=1,markersize=1)
  else:
    for i in range(1,len(ndatau[0,:])):
      # ax[0].plot(-ndatau[:,i]/ry2ev,(x-shift)*ry2ev,linestyle='-', color=colors[i-1], label=legends[i-1],marker=1,markersize=1)
      ax[0].plot(-ndatau[:,i]/ry2ev,(x-shift)*ry2ev,linestyle='-', color=colors[i-1], label=legends[i-1],marker=1,markersize=1)
    for i,line in enumerate(x):
      if abs((line-shift) - fermi) == min(abs((x-shift)-fermi)):
        print(f"n_up(Ef) = {ndatau[i,1]/ry2ev}")

  # Minority spin LDOS
  if mag:
    x = ndatad[:,0]
    if args.onlyS:
      ax[2].plot(ndatad[:,2]/ry2ev,(x-shift)*ry2ev,linestyle='-', color=colors[2-1], label=legends[2-1],marker=1,markersize=1)
      if args.superconductivity:
        ax[2].plot(ndatad[:,5]/ry2ev,(x-shift)*ry2ev,linestyle='-', color=colors[5-1], label=legends[5-1],marker=1,markersize=1)
        ax[2].plot(ndatad[:,5]/ry2ev + ndatad[:,2]/ry2ev,(x-shift)*ry2ev,linestyle='-', color=colors[1-1], label="sum",marker=1,markersize=1)
    else:
      for i in range(1,len(ndatad[0,:])):
        # ax[2].plot(ndatad[:,i]/ry2ev,(x-shift)*ry2ev,linestyle='-', color=colors[i-1],marker=1,markersize=1)
        ax[2].plot(ndatad[:,i]/ry2ev,(x-shift)*ry2ev,linestyle='-', color=colors[i-1], label=legends[i-1],marker=1,markersize=1)

      for i,line in enumerate(x):
        if abs((line-shift) - fermi) == min(abs((x-shift)-fermi)):
          print(f"n_dn(Ef) = {ndatad[i,1]/ry2ev}")
    a2 = max(ndatad[:,1]/ry2ev)


  # Setting limits:
  max_ldos = max(a1,a2) if mag else a1
  xlim = 1.1*abs(max_ldos)
  ax[0].set_xlim([-xlim,0.0])
  if args.ylim != "":
      ylim = eval(args.ylim)
      ax[1].set_ylim(ylim)
  else:
      ax[1].set_ylim([(x-shift)[0]*ry2ev,(x-shift)[-1]*ry2ev])
  if mag:
    ax[2].set_xlim([0.0,xlim])

  # ax[2].axhline(y=1000.0, xmin=point[0], xmax=point[npoints-1], color='b', linestyle='--', linewidth=0.75)
  # ax[2].axhline(y=-1000.0, xmin=point[0], xmax=point[npoints-1], color='b', linestyle='--', linewidth=0.75)
  #
  # ax[0].axhline(y=1000.0, xmin=point[0], xmax=point[npoints-1], color='b', linestyle='--', linewidth=0.75)
  # ax[0].axhline(y=-1000.0, xmin=point[0], xmax=point[npoints-1], color='b', linestyle='--', linewidth=0.75)

  ax[0].set_zorder(1)
  ax[1].set_zorder(0)
  if args.noef:
    ax[0].axhline(y=fermi_ldos, xmin=point[0], xmax=point[npoints-1], color='k', linestyle='--', linewidth=0.75)
  # Getting position of the middle of the band structure frame to put legend from LDOS
  posx = ax[1].get_position().x0 + 0.5*ax[1].get_position().width
  posy = ax[1].get_position().y0 + 0.05*ax[1].get_position().height
  leg = ax[0].legend(loc="center", bbox_to_anchor=(posx,posy), bbox_transform=ax[1].transAxes)
  leg.set_zorder(2)
  ax[0].tick_params(axis='both', which='major', labelsize=14)
  ax[1].tick_params(axis='both', which='major', labelsize=14)
  if mag:
    if args.noef:
      ax[2].axhline(y=fermi_ldos, xmin=point[0], xmax=point[npoints-1], color='k', linestyle='--', linewidth=0.75)
    ax[2].tick_params(axis='both', which='major', labelsize=14)

  if args.output == "":
    plt.show()
  else:
    plt.savefig(args.output,dpi=200)
    if args.show:
      plt.show()
