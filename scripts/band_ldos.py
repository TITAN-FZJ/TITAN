import numpy as np
import matplotlib.pyplot as plt
#from matplotlib2tikz import save as tikz_save
import matplotlib.gridspec as gridspec
import sys

nSites=2
#plt.xkcd()


bsstruct = sys.argv[1]

if(len(sys.argv) == 4):
  ldosu = sys.argv[2]
  ldosd = sys.argv[3]

if(len(sys.argv) == 4):
  fig, ax = plt.subplots(1,3, sharey=True,gridspec_kw = {'width_ratios':[1,4,1]})
  fig.subplots_adjust(wspace=0)
else:
  fig,axx = plt.subplots(1,1)
  ax = [None, axx, None]

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
  ax[1].axvline(x=i, color='k')
 ax[1].axhline(y=0.0, xmin=point[0], xmax=point[count[0]-1], color='k', linestyle='--')
 #ax[1].axhline(y=0.0, xmin=point[0], xmax=point[count[0]-1], color='k', linestyle='--')

 # for i in range(2,19):
 #  ax2.scatter(table[:,0],[(a-fermi) for a in table[:,i]], marker='.', c='k', s=0.2)
 ax[1].scatter(table[:,0],[(a-fermi) for a in table[:,1]], marker='.', c='b', s=0.2, label="TITAN")
 ax[1].set_ylim([-0.8,0.7])
 # ax2.set_title("TITAN")
 ax[1].set_xlim([point[0],point[count[0]-1]])

if(len(sys.argv)==4):
  ax[0].set_ylabel("$E-E_F$ [Ry]")
else:
  ax[1].set_ylabel("$E-E_F$ [Ry]")

if(len(sys.argv)==4):
  data = [ np.loadtxt(ldosu), np.loadtxt(ldosd) ]

  dat1 = (data[0])[(data[0])[:,0].argsort()]
  dat2 = (data[1])[(data[1])[:,0].argsort()]

  x1 = [ x-fermi for x in dat1[:,0] ]
  tot1 = dat1[:,1]
  s1 = dat1[:,2]
  p1 = dat1[:,3]
  t2g1 = dat1[:,4]
  eg1 = dat1[:,5]
  d1 = [x[0]+x[1] for x in zip(dat1[:,4],dat1[:,5]) ]

  x2 = [ x-fermi for x in dat2[:,0] ]
  tot2 = [-x for x in dat2[:,1] ]
  s2 = [-x for x in dat2[:,2] ]
  p2 = [-x for x in dat2[:,3] ]
  t2g2 = [-x for x in dat2[:,4] ]
  eg2 = [-x for x in dat2[:,5] ]
  d2 = [-(x[0]+x[1]) for x in zip(dat2[:,4],dat2[:,5]) ]

  ax[2].plot(tot1,x1,color='k')
  ax[2].plot(s1,x1,color='r')
  ax[2].plot(p1,x1,color='b')
  # ax[2].plot(t2g1,x1,color='y')
  # ax[2].plot(eg1,x1,color='c')
  ax[2].plot(d1,x1,color='c')

  ax[0].plot(tot2,x2,color='k')
  ax[0].plot(s2,x2,color='r')
  ax[0].plot(p2,x2,color='b')
  # ax[0].plot(t2g2,x2,color='y')
  # ax[0].plot(eg2,x2,color='c')
  ax[0].plot(d2,x2,color='c')

  ax[2].axhline(y=0.0, xmin=point[0], xmax=point[count[0]-1], color='k', linestyle='--')
  ax[0].axhline(y=0.0, xmin=point[0], xmax=point[count[0]-1], color='k', linestyle='--')

plt.legend()
#tikz_save("plot.tex", figurewidth="10cm", figureheight="10cm")
plt.show()
