import numpy as np
import matplotlib.pyplot as plt
#from matplotlib2tikz import save as tikz_save
import matplotlib.gridspec as gridspec
import sys
from matplotlib import rc
rc('mathtext', default='regular')

# Fonts
font = {'family': 'arial',
         'color': 'black',
         'weight':'normal',
         'size': 9,
        }

nSites=1
#plt.xkcd()

ry2ev = 13.6
if(ry2ev != 1.0):
  label = r'$E-E_F$ [eV]'
else:
  label = r'$E-E_F$ [Ry]'

# colors = brewer2mpl.get_map('Set1', 'qualitative', 5).mpl_colors
colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#a65628']
# colors = np.array(['k', 'g', 'r', 'b']) # The colors you wish to cycle through
legends = np.array(['Total', 's', 'p', 'd'])

bsstruct = sys.argv[1]

if(len(sys.argv) == 4):
  ldosu = sys.argv[2]
  ldosd = sys.argv[3]
  filename = 'BS_LDOS.pdf'
  fig, ax = plt.subplots(1,3, sharey=True, gridspec_kw = {'width_ratios':[1,4,1]})
  fig.subplots_adjust(wspace=0.15)
  ax[0].tick_params(axis='y', direction='in', left=True, right=True)
  ax[0].set_ylabel(label)
  ax[2].tick_params(axis='y', direction='in', left=True, right=True)
else:
  fig,axx = plt.subplots(1,1)
  ax = [None, axx, None]
  filename = 'BS.pdf'
  ax[1].set_ylabel(label)


fermi = 0.0

with open(bsstruct, "r") as f:
 count = np.array([int(s) for s in f.readline().split()])
 fermi = float(f.readline())
 name = np.empty([count[0]], dtype=str)
 point = np.empty([count[0]], dtype=float)
 for i in range(count[0]):
  name[i], point[i] = f.readline().split()

 table = np.loadtxt(f)
 data = []

 for i, dat in enumerate(table):
  #if( ( i % 25 == 0 )):
  for j in range(18*nSites):
   data.append([dat[0],dat[1+j]])

 table = np.array(data)

 ax[1].set_xticks(point)
 ax[1].set_xticklabels(name)
 for i in point:
  ax[1].axvline(x=i, color='k', linewidth=0.75)
 ax[1].axhline(y=0.0, xmin=point[0], xmax=point[count[0]-1], color='k', linestyle='--', linewidth=0.75)
 #ax[1].axhline(y=0.0, xmin=point[0], xmax=point[count[0]-1], color='k', linestyle='--')

 # for i in range(2,19):
 #  ax2.scatter(table[:,0],[(a-fermi) for a in table[:,i]], marker='.', c='k', s=0.2)
 ax[1].scatter(table[:,0],[(a-fermi)*ry2ev for a in table[:,1]], marker='.', c='k', s=1.0)

 if(ry2ev != 1.0):
  ax[1].set_ylim([-12.0,9])
  ax[1].set_yticks([-12, -8, -4, 0, 4, 8])
 else:
  ax[1].set_ylim([-1.0,0.7])

 ax[1].tick_params(axis='y', direction='in', left=True, right=True)

 # ax2.set_title("TITAN")
 ax[1].set_xlim([point[0],point[count[0]-1]])

if(len(sys.argv)==4):

  datau = np.loadtxt(ldosu)
  datau = datau[datau[:,0].argsort()]
  x = datau[:,0]
  for i in range(1,len(datau[0,:])):
    ax[0].plot(-datau[:,i]/ry2ev,(x-fermi)*ry2ev, color=colors[i-1], label=legends[i-1])

  datad = np.loadtxt(ldosd)
  datad = datad[datad[:,0].argsort()]
  x = datad[:,0]
  for i in range(1,len(datad[0,:])):
    ax[2].plot(datad[:,i]/ry2ev,(x-fermi)*ry2ev, color=colors[i-1])

  a1 = max(datau[:,1]/ry2ev)
  a2 = max(datad[:,1]/ry2ev)
  max_ldos = max(a1,a2)
  xlim = 1.1*abs(max_ldos)
  ax[2].set_xlim([0.0,xlim])
  ax[0].set_xlim([-xlim,0.0])

  ax[2].axhline(y=0.0, xmin=point[0], xmax=point[count[0]-1], color='k', linestyle='--', linewidth=0.75)
  ax[0].axhline(y=0.0, xmin=point[0], xmax=point[count[0]-1], color='k', linestyle='--', linewidth=0.75)

  ax[0].legend(loc=3, prop={'size': 7})

#tikz_save("plot.tex", figurewidth="10cm", figureheight="10cm")
plt.savefig(filename)

plt.show()
