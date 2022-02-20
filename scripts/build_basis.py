#!/usr/bin/env python3
"""
This program reads a 'basis' file (in the POSCAR format)
and repeat it in the direction of the Bravais vectors a1, a2 and/or a3.
The number of repetitions is chosen via the arguments --a1, --a2 and --a3,
which default to 1 for the original cell.
The output basis is put by default in the local file 'basis_out',
but can be changed via the '--output' argument.

@author: Filipe Guimaraes (f.guimaraes@fz-juelich.de)
Feb. 19, 2022

usage: build_basis.py [-h] [--a1 A1] [--a2 A2] [--a3 A3]
                      [--name NAME] [--output OUTPUT]
                      file

Basis creation script

positional arguments:
  file             Basis file (POSCAR style)

optional arguments:
  -h, --help       show this help message and exit
  --a1 A1          Number of repetitions in the direction a1 (default: 1, original cell)
  --a2 A2          Number of repetitions in the direction a2 (default: 1, original cell)
  --a3 A3          Number of repetitions in the direction a3 (default: 1, original cell)
  --name NAME      Name in the output basis (default: same as input)
  --output OUTPUT  Name of output basis file (default: ./basis_out)
  """

import numpy as np


def get_arguments():
  """
  This subroutine parse the arguments from the command line.
  """
  import argparse
  # Parse arguments
  parser = argparse.ArgumentParser(description="Basis creation script")
  parser.add_argument("file", help="Basis file (POSCAR style)")
  # parser.add_argument("--dim", default=3, type=int, help="Dimension of the system (default: dim=3)")
  parser.add_argument("--a1", default=1, type=int, help="Number of repetitions in the direction a1 (default: 1, original cell)")
  parser.add_argument("--a2", default=1, type=int, help="Number of repetitions in the direction a2 (default: 1, original cell)")
  parser.add_argument("--a3", default=1, type=int, help="Number of repetitions in the direction a3 (default: 1, original cell)")
  parser.add_argument("--name", default="", help="Name in the output basis (default: same as input)")
  parser.add_argument("--output", default="./basis_out", help="Name of output basis file (default: ./basis_out)")
  # parser.add_argument("--daemon", default=False, action="store_true" , help="Run as a 'daemon', i.e., in an infinite loop")
  return parser.parse_args()

def read_basis_file(filename):
  # Reading basis file
  with open(filename,'r') as f:
    file = f.read().split("\n")

  # Parsing the quantities
  basis = {}
  basis['name'] = file.pop(0)
  basis['a0'] = float(file.pop(0))
  basis['a1'] = np.array([float(_) for _ in file.pop(0).split()])
  basis['a2'] = np.array([float(_) for _ in file.pop(0).split()])
  basis['a3'] = np.array([float(_) for _ in file.pop(0).split()])
  basis['types'] = file.pop(0).split()
  basis['natoms'] = np.array([int(_) for _ in file.pop(0).split()])
  basis['units'] = file.pop(0)
  basis['positions'] = []
  for natoms in basis['natoms']:
    # basis['positions'].append([])
    for natom in range(natoms):
      basis['positions'].append(np.array([float(_) for _ in file.pop(0).split()]))
  return basis

def repeat_cell(basis_in,
                a1,
                a2,
                a3,
                name):
  """
  This function receives a dictionary with the original basis and cell
  and repeats it a1, a2, a3 times in the direction of each Bravais
  vector, respectively. The output system name is given in 'name'.
  """
  basis_out = {}
  basis_out['name'] = (name if name !="" else basis_in['name'])
  basis_out['a0'] = basis_in['a0']
  basis_out['a1'] = basis_in['a1']*a1
  basis_out['a2'] = basis_in['a2']*a2
  basis_out['a3'] = basis_in['a3']*a3
  basis_out['types'] = basis_in['types']
  basis_out['natoms'] = basis_in['natoms']*a1*a2*a3
  basis_out['units'] = basis_in['units']
  basis_out['positions'] = []
  for position in basis_in['positions']:
    for rep_a3 in range(a3):
      for rep_a2 in range(a2):
        for rep_a1 in range(a1):
          basis_out['positions'].append(position+basis_in['a1']*rep_a1+basis_in['a2']*rep_a2+basis_in['a3']*rep_a3)
  return basis_out

def write_basis(basis_out,
                output):
  """
  This function receives a dictionary with the repeated cell
  and outputs it in the 'output' file.
  """
  with open(output,'w') as f:
    f.write(f"{basis_out['name']}\n")
    f.write(f"{basis_out['a0']:.8f}\n")
    f.write(f"{'  '.join(f'{_:.8f}' for _ in basis_out['a1'])}\n")
    f.write(f"{'  '.join(f'{_:.8f}' for _ in basis_out['a2'])}\n")
    f.write(f"{'  '.join(f'{_:.8f}' for _ in basis_out['a3'])}\n")
    f.write(f"{'  '.join(basis_out['types'])}\n")
    f.write(f"{'  '.join(f'{_}' for _ in basis_out['natoms'])}\n")
    f.write(f"{basis_out['units']}\n")
    for position in basis_out['positions']:
      f.write(f"{'  '.join(f'{_:.8f}' for _ in position)}\n")
  return

def main():
  """
  Main function with the calls for reading, repeating and writing out.
  """
  args = get_arguments()
  # dim = args.dim

  basis_in = read_basis_file(args.file)
  
  basis_out = repeat_cell(basis_in,
                              args.a1,
                              args.a2,
                              args.a3,
                              args.name)
  write_basis(basis_out,args.output)


if __name__ == "__main__":
  main()
