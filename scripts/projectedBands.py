################################################################################
# Import the arguments and style
################################################################################
from head import *

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

    ##### jump garbage line
    f.readline()

    dimbs_line = f.readline().split()
    dimbs = int(dimbs_line[2])

  return npoints, name_out, point, fermi, dimbs

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
#
# IMPORTANT: If you want to plot the superconductor bands, then the way to
# execute this script is
#
# $ ipython < route to script.py> <route to datafile> -- --superconductivity
#
# If the user wants to have the electron and hole bands in the same plot
# then, the script has to be run as
# $ ipython < route to script.py> <route to datafile> -- --superconductivity -- together
#
# In the next example we show how to use all options, we use \\ to break the line,
# and make it easier to read
#
# ipython  < route to script.py> <route to datafile> -- \\
#           --superconductivity \\
#           --title="Projection of \$p_z\$" \\
#           --ev \\
#           --guide \\
#           --output="cool name" \\
#           --zoom=0.1
#           --project=[[1,19],[2,20],[4,24]]
#           --separate
#
# --separate does not work together with all the rest for the moment

if __name__ == "__main__":

  if args.superconductivity:
    if args.separate:
      numplots = len(sys.argv)-1
      titles = [r"Electron Bands", r"Hole Bands"]

      fig, axs = plt.subplots(1, numplots, sharey=True, squeeze=False, figsize=(6*numplots, 5))
      if args.ev == True:
        axs[0,0].set_ylabel("Energy [eV]")
      else:
        axs[0,0].set_ylabel("Energy [Ry]")

      for i in range(numplots):
        npoints, name, point, fermi = read_header(args.files[0])
        table = read_data(args.files[0])
        if(len(args.project) > 0):
          weights = read_data(args.files[1])

        axs[0,i].set_title(titles[i])
        axs[0,i].set_xlim([point[0],point[npoints-1]])

        axs[0,i].set_xticks(point)
        axs[0,i].set_xticklabels(name)
        for j in point:
          axs[0,i].axvline(x=j, color='k', linewidth=0.5)

        # Ploting the Fermi level or a line at y=0.0
        if args.noef:
          if (fermi == None): # susceptibility
            axs[0,i].axhline(y=0.0, xmin=point[0], xmax=point[npoints-1], color='k', linestyle='-', linewidth=0.5)
          else: # band structure
            axs[0,i].axhline(y=0.0, xmin=point[0], xmax=point[npoints-1], color='k', linestyle='--', linewidth=0.5)
            if(i==0):
              axs[0,i].axhline(y=fermi, xmin=point[0], xmax=point[npoints-1], color='b', linestyle='--', linewidth=0.5)
            else:
              axs[0,i].axhline(y=-fermi, xmin=point[0], xmax=point[npoints-1], color='b', linestyle='--', linewidth=0.5)

        # Plotting the results
        if (fermi == None): # susceptibility
          axs[0,i].plot(table[:,1],-1.0/table[:,2])
        else: # band structure
          if(i==0):
            axs[0,i].plot(table[:,0],table[:,(table.shape[1]-1)/2+1:], color='r', linewidth=1.0, linestyle='-')
          else:
            axs[0,i].plot(table[:,0],table[:,1:(table.shape[1]-1)/2], color='k', linewidth=1.0, linestyle='-')

    else: # args.separate
      #legacy code
      numplots = 1
      titles = [r"Band Structure"]
      fig, axs = plt.subplots(1, numplots, sharey=True, squeeze=False, figsize=(6*numplots, 5))

      # Setting the units
      if args.ev == True:
        axs[0,0].set_ylabel("Energy [eV]")
      else:
        axs[0,0].set_ylabel("Energy [Ry]")

      for i in range(numplots):
        # Reads and stores the high symmetry points
        npoints, name, point, fermi, dimbs = read_header(args.files[0])
        # Reads the data points
        table = read_data(args.files[0])
        # Reads the weigths for projected band structure
        if(len(args.project) > 0):
          weights = read_data(args.files[1])
        # Set custom title or default title
        if args.title != "":
          axs[0,i].set_title(args.title)
        else:
          axs[0,i].set_title(titles[i])

        axs[0,i].set_xlim([point[0],point[npoints-1]])
        # Used to zoom in a symmetrical region around zero
        if args.zoom != 0.0:
          axs[0,i].set_ylim(-args.zoom,args.zoom)

        axs[0,i].set_xticks(point)
        axs[0,i].set_xticklabels(name)

        # High symmetry points vertical lines
        for j in point:
          axs[0,i].axvline(x=j, color='k', linewidth=linewidth)

        # Ploting the Fermi level or a line at y=0.0
        if args.noef:
          if (fermi == None): # susceptibility
            axs[0,i].axhline(y=0.0, xmin=point[0], xmax=point[npoints-1], color='k', linestyle='-', linewidth=0.5)
          else: # band structure
            axs[0,i].axhline(y=0.0, color='k', linestyle='--',linewidth=0.5)

        # Wider lines for the frame
        for axis in ['top','bottom','left','right']:
          axs[0,i].spines[axis].set_linewidth(linewidth)

        # Plotting the results
        if (fermi == None): # susceptibility
          axs[0,i].plot(table[:,1],-1.0/table[:,2])
        else: # band structure
          # Print band structures as guides
          if args.guide == True:
            for j in range(1,dimbs+1):
              axs[0,i].scatter(table[:,0],table[:,j]*ry2ev,c='k',s=0.01,marker="o")

          for g, group in enumerate(eval(args.project)):
            if isinstance(group, int):
              orb = group
              for j in range(1,dimbs+1):
                weight = weights[:,orb+dimbs*(j-1)-1]
                axs[0,i].scatter(table[:,0],table[:,j]*ry2ev,c=[(0.07, 0.19, 0.40, a) for a in weight],s=10)
            else:
              for k in group:
                orb = k
                for j in range(1,dimbs+1):
                  weight = weights[:,orb+dimbs*(j-1)-1]
                  axs[0,i].scatter(table[:,0],table[:,j]*ry2ev,c=[(colors[g][0],colors[g][1],colors[g][2], a) for a in weight],s=8)

  else: # args.superconductivity:
    numplots = 1
    titles = [r"Band Structure"]
    fig, axs = plt.subplots(1, numplots, sharey=True, squeeze=False, figsize=(6*numplots, 5))
    # Setting the units
    if args.ev == True:
      axs[0,0].set_ylabel(r"E - E$_\mathrm{F}$ [eV]")
    else:
      axs[0,0].set_ylabel(r"E - E$_\mathrm{F}$ [Ry]")

    for i in range(numplots):
      # Reads and stores the high symmetry points
      npoints, name, point, fermi, dimbs = read_header(args.files[0])
      # Reads the data points
      table = read_data(args.files[0])
      # Reads the weigths for projected band structure
      if(len(args.project) > 0):
        weights = read_data(args.files[1])
      # Set custom title or default title
      if args.title != "":
        axs[0,i].set_title(args.title)
      else:
        axs[0,i].set_title(titles[i])

      axs[0,i].set_xlim([point[0],point[npoints-1]])
      # Used to zoom in a symmetrical region around Ef
      if args.zoom != 0.0:
        # print(fermi*ry2ev-args.zoom,fermi*ry2ev+args.zoom)
        axs[0,i].set_ylim(-args.zoom,args.zoom)

      axs[0,i].set_xticks(point)
      axs[0,i].set_xticklabels(name)

      # High symmetry points vertical lines
      for j in point:
        axs[0,i].axvline(x=j, color='k', linewidth=linewidth)

      # Ploting the Fermi level or a line at y=0.0
      if args.noef:
        if (fermi == None): # susceptibility
          axs[0,i].axhline(y=0.0, xmin=point[0], xmax=point[npoints-1], color='k', linestyle='-', linewidth=0.5)
        else: # band structure
          axs[0,i].axhline(y=0.0, color='k', linestyle='--',linewidth=0.5)

      # Wider lines for the frame
      for axis in ['top','bottom','left','right']:
        axs[0,i].spines[axis].set_linewidth(linewidth)

      # Plotting the results
      if (fermi == None): # susceptibility
        axs[0,i].plot(table[:,1],-1.0/table[:,2])
      else: # band structure
        # Print band structures as guides
        if args.guide == True:
          for j in range(1,dimbs+1):
            axs[0,i].scatter(table[:,0],(table[:,j]-fermi)*ry2ev,c='k',s=0.01,marker="o")

        for g, group in enumerate(eval(args.project)):
          if isinstance(group, int):
            orb = group
            for j in range(1,dimbs+1):
              weight = weights[:,orb+dimbs*(j-1)-1]
              axs[0,i].scatter(table[:,0],(table[:,j]-fermi)*ry2ev,c=[(0.07, 0.19, 0.40, a) for a in weight],s=10)
          else:
            for k in group:
              orb = k
              for j in range(1,dimbs+1):
                weight = weights[:,orb+dimbs*(j-1)-1]
                axs[0,i].scatter(table[:,0],(table[:,j]-fermi)*ry2ev,c=[(colors[g][0],colors[g][1],colors[g][2], a) for a in weight],s=8)

  plt.tight_layout()
  if args.output == "":
    plt.show()
  else:
    plt.savefig(args.output,dpi=200)
    if args.show:
      plt.show()
