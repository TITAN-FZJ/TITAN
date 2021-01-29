################################################################################
# This script with all its options is to be run like this
# ipython <script_location> <band_file> <ldosu_file> <ldosd_file> -- --superconductivity --mev --title="desired title"
# the arguments starting with "--" are optional and can be ommited
################################################################################

################################################################################
# Import the arguments and style
################################################################################
from head import *

if(ry2ev == 13.6057*1000):
  labely = r'$E-E_F$ [meV]'
elif(ry2ev == 13.6057):
  labely = r'$E-E_F$ [eV]'
else:
  labely = r'$E-E_F$ [Ry]'

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
    name = np.empty(npoints, dtype=str)
    point = np.empty(npoints, dtype=float)
    for i in range(npoints):
      name[i], point[i] = f.readline().split()[1:]
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

  bsstruct = args.files[0]#args.fileband

  ldosu = args.files[1] #args.fileu
  ldosd = args.files[2] #args.filed
  filename = args.output
  fig, ax = plt.subplots(1,3, sharey=True, gridspec_kw = {'width_ratios':[1,4,1]})
  fig.subplots_adjust(left=0.15,wspace=0.15)

  if args.title != "":
      fig.suptitle(args.title,fontsize=24)

  ax[0].tick_params(axis='y', direction='in', left=True, right=True)
  ax[0].set_ylabel(labely, fontsize=14)
  ax[2].tick_params(axis='y', direction='in', left=True, right=True)

  npoints, name, point, fermi = read_header(args.files[0])
  table = read_data(args.files[0])
  fermi_ldos = read_header_ldos(args.files[1])

  ndatau=read_data_ldos(args.files[1])
  ndatad=read_data_ldos(args.files[2])

  ndatau = ndatau[ndatau[:,0].argsort()]
  ndatad = ndatad[ndatad[:,0].argsort()]

  ax[1].set_xticks(point)
  ax[1].set_xticklabels(name)
  for i in point:
      ax[1].axvline(x=i, color='k', linewidth=0.75)
  ax[1].axhline(y=0.0, color='k', linestyle='--', linewidth=0.75)

  if args.gap !=0.0:
    for i in range(len(ax)):
      ax[i].axhline(y=-args.gap, color='r', linestyle='-.', linewidth=0.75)
      ax[i].axhline(y= args.gap, color='r', linestyle='-.', linewidth=0.75)

  if (fermi == None): # susceptibility
    ax[1].plot(table[:,1],-1.0/table[:,2])
  else: # band structure
    if args.superconductivity:
        ax[1].plot(table[:,0],table[:,(table.shape[1]-1)/2+1:]*ry2ev-fermi_ldos*ry2ev, color='r', linewidth=1.0, linestyle='-')
        ax[1].plot(table[:,0],table[:,1:(table.shape[1]-1)/2+1]*ry2ev-fermi_ldos*ry2ev, color='k', linewidth=1.0, linestyle='-')
    else:
        ax[1].plot(table[:,0],table[:,1:]*ry2ev-fermi_ldos*ry2ev, color='k', linewidth=1.0, linestyle='-')

  ax[1].tick_params(axis='y', direction='in', left=True, right=True)

  ax[1].set_xlim([point[0],point[npoints-1]])

  # Majority spin LDOS
  x = ndatau[:,0]
  if args.onlyS:
      ax[0].plot(-ndatau[:,2]/ry2ev,(x-fermi_ldos)*ry2ev,linestyle='-', color=colors[2-1], label=legends[2-1],marker=1,markersize=1)
      if args.superconductivity:
          ax[0].plot(-ndatau[:,5]/ry2ev,(x-fermi_ldos)*ry2ev,linestyle='-', color=colors[5-1], label=legends[5-1],marker=1,markersize=1)
          ax[0].plot(-ndatau[:,5]/ry2ev - ndatau[:,2]/ry2ev,(x-fermi_ldos)*ry2ev,linestyle='-', color=colors[1-1], label="sum",marker=1,markersize=1)
  else:
      for i in range(1,len(ndatau[0,:])):
          # ax[0].plot(-ndatau[:,i]/ry2ev,(x-fermi_ldos)*ry2ev,linestyle='-', color=colors[i-1], label=legends[i-1],marker=1,markersize=1)
          ax[0].plot(-ndatau[:,i]/ry2ev,(x-fermi_ldos)*ry2ev,linestyle='-', color=colors[i-1], label=legends[i-1],marker=1,markersize=1)

  # Minority spin LDOS
  x = ndatad[:,0]
  if args.onlyS:
      ax[2].plot(ndatad[:,2]/ry2ev,(x-fermi_ldos)*ry2ev,linestyle='-', color=colors[2-1], label=legends[2-1],marker=1,markersize=1)
      if args.superconductivity:
          ax[2].plot(ndatad[:,5]/ry2ev,(x-fermi_ldos)*ry2ev,linestyle='-', color=colors[5-1], label=legends[5-1],marker=1,markersize=1)
          ax[2].plot(ndatad[:,5]/ry2ev + ndatad[:,2]/ry2ev,(x-fermi_ldos)*ry2ev,linestyle='-', color=colors[1-1], label="sum",marker=1,markersize=1)
  else:
      for i in range(1,len(ndatad[0,:])):
          # ax[2].plot(ndatad[:,i]/ry2ev,(x-fermi_ldos)*ry2ev,linestyle='-', color=colors[i-1],marker=1,markersize=1)
          ax[2].plot(ndatad[:,i]/ry2ev,(x-fermi_ldos)*ry2ev,linestyle='-', color=colors[i-1], label=legends[i-1],marker=1,markersize=1)


  # Setting limits:
  a1 = max(ndatau[:,1]/ry2ev)
  a2 = max(ndatad[:,1]/ry2ev)

  max_ldos = max(a1,a2)
  xlim = 1.1*abs(max_ldos)
  ax[0].set_xlim([-xlim,0.0])
  if args.ylim != "":
      ylim = eval(args.ylim)
      ax[1].set_ylim(ylim)
  else:
      ax[1].set_ylim([1.1*(x-fermi_ldos)[0]*ry2ev,1.1*(x-fermi_ldos)[-1]*ry2ev])
  ax[2].set_xlim([0.0,xlim])

  #
  # ax[2].axhline(y=1000.0, xmin=point[0], xmax=point[npoints-1], color='b', linestyle='--', linewidth=0.75)
  # ax[2].axhline(y=-1000.0, xmin=point[0], xmax=point[npoints-1], color='b', linestyle='--', linewidth=0.75)
  #
  # ax[0].axhline(y=1000.0, xmin=point[0], xmax=point[npoints-1], color='b', linestyle='--', linewidth=0.75)
  # ax[0].axhline(y=-1000.0, xmin=point[0], xmax=point[npoints-1], color='b', linestyle='--', linewidth=0.75)

  ax[0].axhline(y=0.0, xmin=point[0], xmax=point[npoints-1], color='k', linestyle='--', linewidth=0.75)
  ax[2].axhline(y=0.0, xmin=point[0], xmax=point[npoints-1], color='k', linestyle='--', linewidth=0.75)
  posx = ax[1].get_position().x0 + 0.5*ax[1].get_position().width
  posy = ax[1].get_position().y0 + 0.9*ax[1].get_position().height
  ax[2].legend(loc="center", bbox_to_anchor=(posx,posy), bbox_transform=ax[1].transAxes)


  ax[0].tick_params(axis='both', which='major', labelsize=14)
  ax[1].tick_params(axis='both', which='major', labelsize=14)
  ax[2].tick_params(axis='both', which='major', labelsize=14)

  if args.output == "":
    plt.show()
  else:
    plt.savefig(args.output,dpi=200)
    if args.show:
      plt.show()
