# *************************************************************************************
# *                                                                                   *
# *   PAOFLOW *  Marco BUONGIORNO NARDELLI * University of North Texas 2016-2018      *
# *                                                                                   *
# *************************************************************************************
#
#  Copyright 2016-2022 - Marco BUONGIORNO NARDELLI (mbn@unt.edu) - AFLOW.ORG consortium
#
#  This file is part of AFLOW software.
#
#  AFLOW is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# *************************************************************************************

from PAOFLOW import PAOFLOW

def main():
 
  path = 'G-H-N-G-P-N-P-H'
  special_points = {'G':[0.0, 0.0, 0.0],'H':[-0.5, 0.5, 0.5],'N':[0.0, 0.5, 0.0],'P':[0.25,0.25,0.25]}

  paoflow = PAOFLOW.PAOFLOW(savedir='nb.save', smearing='m-p')
  paoflow.read_atomic_proj_QE()
  paoflow.projectability()
  paoflow.pao_hamiltonian()
  paoflow.bands(ibrav=3, nk=1500, band_path=path, high_sym_points=special_points)
  paoflow.interpolated_hamiltonian()
  paoflow.pao_eigh()
  paoflow.gradient_and_momenta()
  paoflow.adaptive_smearing()
  paoflow.dos(emin=-20., emax=40., delta=.2)
  paoflow.z2_pack()
  paoflow.finish_execution()

if __name__== '__main__':
  main()

