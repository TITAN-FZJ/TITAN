import numpy as np 
import matplotlib.pyplot as plt 
import sys

#plt.xkcd()

fig, (ax1,ax2) = plt.subplots(1,2, sharey=True)
ax1.set_ylabel("Energy [Ry]")
with open(sys.argv[1], "r") as f:
 count = np.array([int(s) for s in f.readline().split()])

 fermi = float(f.readline())

 name = np.empty([count[0]], dtype=str)
 point = np.empty([count[0]], dtype=float)
 for i in range(count[0]):
  name[i], point[i] = f.readline().split()

 table = np.loadtxt(f)
 ax1.set_title("Band Structure w/o SOC")
 ax1.set_xlim([point[0],point[count[0]-1]])

 ax1.set_xticks(point)
 ax1.set_xticklabels(name)
 for i in point:
  ax1.axvline(x=i, color='k')
 ax1.axhline(y=fermi, xmin=point[0], xmax=point[count[0]-1], color='k', linestyle='--')

 for i in range(1,10):
  ax1.scatter(table[:,0],table[:,i], marker='.', c='r', s=0.1)
 for i in range(10,19):
  ax1.scatter(table[:,0],table[:,i], marker='.', c='b', s=0.1)

with open(sys.argv[2], "r") as f:
 count = np.array([int(s) for s in f.readline().split()])

 fermi = float(f.readline())

 name = np.empty([count[0]], dtype=str)
 point = np.empty([count[0]], dtype=float)
 for i in range(count[0]):
  name[i], point[i] = f.readline().split()

 table = np.loadtxt(f)
 ax2.set_title("Band structure with SOC 0.2")
 ax2.set_xlim([point[0],point[count[0]-1]])

 ax2.set_xticks(point)
 ax2.set_xticklabels(name)
 for i in point:
  ax2.axvline(x=i, color='k')
 ax2.axhline(y=fermi, xmin=point[0], xmax=point[count[0]-1], color='k', linestyle='--')

 for i in range(1,19):
  ax2.scatter(table[:,0],table[:,i], marker='.', c='k', s=0.1)


plt.show()
