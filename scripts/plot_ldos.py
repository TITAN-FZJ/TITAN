import numpy as np
import matplotlib.pyplot as plt
import sys

#plt.xkcd()

Ef= 0.81848251

ry2ev = 1.0
if(ry2ev != 1.0):
  labelx = "$E-E_F$ [eV]"
  labely = 'LDOS [states/eV]'
else:
  labelx = "$E-E_F$ [Ry]"
  labely = 'LDOS [states/Ry]'

filenameu = sys.argv[1]
filenamed = sys.argv[2]

colors = np.array(['k', 'g', 'r', 'b'])	# The colors you wish to cycle through
legends = np.array(['Total', 's', 'p', 'd'])

plt.axhline(0.0, color='k')
if(Ef != 0.0):
	plt.axvline(0.0, color='k', linestyle='--')

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


plt.xlabel(labelx)
plt.ylabel(labely)
plt.title('LDOS')
plt.legend()

plt.show()

