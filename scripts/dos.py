import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
from matplotlib import rc
rc('mathtext', default='regular')
from matplotlib.font_manager import FontProperties
from os.path import sep

# Fonts
font = {'family': 'arial',
         'color': 'black',
         'weight':'normal',
         'size': 9,
        }

Ef=0.776521106614553

def get_energy(col_path,atom,fermi):
	filename  = 'dos.atom' + str(atom)
	filepath  = col_path + sep + filename
	ldos      = np.loadtxt(filepath)
	x         =  (ldos[:,0]- fermi)*13.6058
	energy = np.array_split(x,2)
	return energy[0]

def get_ldos(path,atom):
	ldos=[]
	filename = 'dos.atom' + str(atom)
    	filepath = path + sep + filename
    	data     = np.loadtxt(filepath)
	x = data[:,1]/13.6058
	a1 = max(x)
	a2 = min(x)
	b = max(a1,abs(a2))
	ldos = np.array_split(x, 2)
	return ldos, b

########################################################################
# Set colorbrewer colors for the plot
########################################################################
#set2 = brewer2mpl.get_map('Set1', 'qualitative', 5).mpl_colors

set2 = ['#053061','#b2182b']
#set2 = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#a65628']

########################################################################
# Plot the LDOS for the Majority and Minority states 
########################################################################
path='.'

pos_atom=[38]
pos=[r'l$_{max}$ = 3',r'l$_{max}$ = 4']

energy = get_energy(path,pos_atom[0],Ef)
ldos, max_ldos = get_ldos(path,pos_atom[0])

plt.plot(energy, ldos[0], ls='-', c=set2[0],  lw=0.7, label=pos[0])
plt.plot(energy, ldos[1], ls='-', c=set2[0],  lw=0.7)
plt.plot([-15,15], [0,0], linestyle='-', color='black', linewidth=0.5)
plt.plot([0,0], [-15,15], linestyle='-', color='black', linewidth=0.5)

ylim=[-1.1*abs(max_ldos),1.1*abs(max_ldos)]
xlim=[min(energy),max(energy)]

plt.xlim(xlim)
plt.ylim(ylim)

plt.xlabel('E - E$_F$ (eV)', fontdict=font)
plt.ylabel('LDOS (states/eV)', fontdict=font)
plt.legend(loc='best',frameon=False, prop={'size':7}, labelspacing=0.12)
plt.tick_params(axis='x', colors='black',labelsize=8,width=1)
plt.tick_params(axis='y', colors='black',labelsize=8,width=1)

fig = matplotlib.pyplot.gcf()
fig.set_size_inches(3., 2.2)
filename= 'plots_dos/LDOS-atom' + str(pos_atom[i]) +  '.pdf'
plt.savefig(filename,transparent=True,dpi=2540, bbox_inches='tight')
plt.show()

