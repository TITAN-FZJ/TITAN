import json
import numpy as np
import sys

data = np.array(json.loads(sys.argv[1]))

basis = data[:3]
# basis = [
#             np.array([-0.5, 0.5, 0.5]),
#             np.array([0.5, -0.5, 0.5]),
#             np.array([0.5, 0.5, -0.5])
# ]
pln = data[3]
#pln = np.array([0.0, -1.0, 1.0])
pln = pln / np.linalg.norm(pln)

def rotateByN(a, n):
  zaxis = np.array([0.0,0.0,1.0])
  if(np.linalg.norm(np.cross(zaxis,n)) < 1e-9):
    return a

  cross = np.cross(n,zaxis)
  cross = cross / np.linalg.norm(cross)
  theta = np.arccos(np.dot(zaxis,n))

  p = [np.cos(theta*0.5), cross[0] * np.sin(theta*0.5), cross[1] * np.sin(theta*0.5), cross[2] * np.sin(theta*0.5)]
  d = [0.0, a[0], a[1], a[2]]
  q = [np.cos(theta*0.5), -cross[0] * np.sin(theta*0.5), -cross[1] * np.sin(theta*0.5), -cross[2] * np.sin(theta*0.5)]

  tmp = HamiltonianProduct(p,d)
  tmp2 = HamiltonianProduct(tmp,q)

  return tmp2[1:]


def HamiltonianProduct(a,b):
  first = a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3]
  second = a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2]
  third = a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1]
  fourth = a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0]
  return [first, second, third, fourth]

atoms = []


for m in range(-3,4):
  for n in range(-3,4):
    for l in range(-3,4):
      atoms.append((m*basis[0] + n*basis[1] + l*basis[2]))

ip = []
oop = []
for atom in atoms:
  if(np.sqrt(np.dot(atom,pln)**2) < 1.0e-8):
    ip.append(atom)
  else:
    oop.append(atom)

ip.sort(key=np.linalg.norm)
ip = ip[1:]
oop.sort(key=np.linalg.norm)

ip_smallest = [rotateByN(a,pln) for a in ip if abs(np.linalg.norm(ip[0]) - np.linalg.norm(a)) < 1e-9]
oop_smallest = [rotateByN(a,pln) for a in oop if abs(np.linalg.norm(oop[0]) - np.linalg.norm(a)) < 1e-9]
print("In plane atoms")
print(ip_smallest)
print("Out of plane atoms")
print(oop_smallest)
