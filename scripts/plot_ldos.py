import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib import rc
rc('mathtext', default='regular')

# Fonts
font = {'family': 'arial',
         'color': 'black',
         'weight':'normal',
         'size': 9,
        }

# plt.xkcd()

# Manually inserted Fermi Energy
Ef= 0.0

ry2ev = 1.0
if(ry2ev != 1.0):
  labelx = r'$E-E_F$ [eV]'
  labely = r'LDOS [states/eV]'
else:
  labelx = r'$E-E_F$ [Ry]'
  labely = r'LDOS [states/Ry]'

filenameu = sys.argv[1]
filenamed = sys.argv[2]

# colors = brewer2mpl.get_map('Set1', 'qualitative', 5).mpl_colors
colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#a65628']
# colors = np.array(['k', 'g', 'r', 'b']) # The colors you wish to cycle through
legends = np.array(['Total', 's', 'p', 'd'])

plt.axhline(0.0, color='k', linewidth=0.75)

datau = np.loadtxt(filenameu)
datau = datau[datau[:,0].argsort()]
x = datau[:,0]
for i in range(1,len(datau[0,:])):
  plt.plot((x-Ef)*ry2ev,datau[:,i]/ry2ev, color=colors[i-1], label=legends[i-1])

datad = np.loadtxt(filenamed)
datad = datad[datad[:,0].argsort()]
x = datad[:,0]
for i in range(1,len(datad[0,:])):
  plt.plot((x-Ef)*ry2ev,-datad[:,i]/ry2ev, color=colors[i-1])

a1 = max(datau[:,1]/ry2ev)
a2 = max(datad[:,1]/ry2ev)
max_ldos = max(a1,a2)
ylim = [-1.1*abs(max_ldos),1.1*abs(max_ldos)]
plt.ylim(ylim)


plt.xlabel(labelx)
plt.ylabel(labely)
plt.title('LDOS')
plt.legend()

# if(Ef != 0.0):
plt.axvline(0.0, color='k', linestyle='--', linewidth=0.75)

# fig = plt.gcf()
# fig.set_size_inches(3., 2.2)
plt.savefig('LDOS.pdf')
# plt.savefig(filename,transparent=True,dpi=2540, bbox_inches='tight')

plt.show()
