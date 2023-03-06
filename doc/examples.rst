.. index:: Examples

********
Examples
********

Thanks to a sporadic collaborator, TITAN includes a list of many examples, that are also used to test correctness and consistency of new implementations.
The tests are listed and detailed here.


.. index:: Example: Susceptibility, Fe bulk

Susceptibility: Iron bulk
=========================


.. index:: Graphene
.. _graphene:

Graphene bulk (2D)
==================


.. index:: Example: Response to electric field, Nanoribbon

Response: Graphene nanoribbon
=============================


.. index:: Example: Time propagation, Ni bulk

Time propagation: Nickel bulk
=============================


.. index:: Co monolayer
.. _co_monolayer:

Cobalt monolayer
================


.. index:: Example: Constraining field, Fe dimer

Constraining field: Iron dimer
==============================

.. index:: Nb monolayer

Niobium monolayer
=================


.. index:: Example: Superconductivity

Superconductivity: Niobium monolayer
====================================

.. index:: Example: Bulk Nb: Quantum Espresso, PAOFLOW and TITAN 

Bulk Nb: Quantum Espresso, PAOFLOW and TITAN
============================================

``-> tbmode = 2`` requires a file with the Hamiltonian obtained from Quantum Espresso (QE) and PAOFLOW(PF). 
In this example we detail the steps to follow to obtain these Hamiltonians and the details we must keep in 
mind while doing so.

Calculations on Quantum Espresso
--------------------------------

The files for the calculations in QE are located in the folder ``./examples/Nb_bulk_paoflow/1_quantum_espresso``.

It is necessary to get the pseudopotential file for QE, it can be found in the following URL 
`Nb pseudopotentials <http://pseudopotentials.quantum-espresso.org/legacy_tables/ps-library/nb>`_. Download the 
one called ``Nb.pbe-spn-kjpaw_psl.1.0.0.UPF`` and put it in the folder ``./examples/Nb_bulk_paoflow/1_quantum_espresso``.

The first step is to run the self consistency, the corresponding file is ``scf.in``

.. code-block:: text
    :linenos:

    &control
        calculation='scf'
        restart_mode='from_scratch',
        verbosity='high',
        prefix='nb',
        pseudo_dir = './'
        outdir='./'
    /
    &system
        ibrav=  3, celldm(1) =6.23610, nat= 1, ntyp= 1,
        ecutwfc =25.0, ecutrho=250.0,
        occupations='smearing', smearing='marzari-vanderbilt', degauss=0.05
        nspin=2, starting_magnetization(1)=0.01
    /
    &electrons
        conv_thr =  1.0d-8
        mixing_beta = 0.7
    /
    ATOMIC_SPECIES
    Nb 92.9064 Nb.pbe-spn-kjpaw_psl.1.0.0.UPF
    ATOMIC_POSITIONS
    Nb 0.00 0.00 0.00
    K_POINTS {automatic}
    48 48 48  0 0 0

We execute this file with the command ``pw.x < scf.in > scf.out``. This calculation should not take 
longer than a couple of minutes. Once it is done we must perform a non-self consistent calculation.
The correponding file is ``nscf.in``

.. code-block:: text
    :linenos:

    &control
        calculation='nscf'
        restart_mode='from_scratch',
        verbosity='high',
        prefix='nb',
        pseudo_dir = './'
        outdir='./'
    /
    &system
        ibrav=  3, celldm(1) =6.23610, nat= 1, ntyp= 1,
        ecutwfc =25.0, ecutrho=250.0,
        occupations='smearing', smearing='marzari-vanderbilt', degauss=0.05
        nspin=2, starting_magnetization(1)=0.1
        nbnd=16
    /
    &electrons
        conv_thr =  1.0d-8
        mixing_beta = 0.7
    /
    ATOMIC_SPECIES
    Nb 92.9064 Nb.pbe-spn-kjpaw_psl.1.0.0.UPF
    ATOMIC_POSITIONS
    Nb 0.00 0.00 0.00
    K_POINTS {automatic}
    12 12 12 0 0 0

We execute this file with the command ``pw.x < nscf.in > nscf.out``. Again it should not take
longer than a couple of minutes. The next necessary step is to generate the projections with ``proj.in``, in this step
we also obtain the density of states

.. code-block:: text
    :linenos:

    &projwfc
        prefix='nb'
        outdir='./'
        filpdos='./nb'
        lwrite_overlaps = .false.
        lbinary_data  = .false.
    /

To execute this file we must type ``projwf.x < proj.in > proj.out``. This is already enough to continue to PAOFLOW.
However it might be a good idea to calculate the band structure and the density of states to cross check our results.

To get the band structure we need again a non-self consistent calculation, the corresponding file is ``bands.scf``

.. code-block:: text
    :linenos:

    &control
        calculation='bands'
        restart_mode='from_scratch',
        verbosity='high',
        prefix='nb',
        pseudo_dir = './'
        outdir='./'
    /
    &system
        ibrav=  3, celldm(1) =6.23610, nat= 1, ntyp= 1,
        ecutwfc =25.0, ecutrho=250.0,
        occupations='smearing', smearing='marzari-vanderbilt', degauss=0.05
        nspin=2, starting_magnetization(1)=0.01
    /
    &electrons
        conv_thr =  1.0d-8
        mixing_beta = 0.7
    /
    ATOMIC_SPECIES
    Nb 92.9064 Nb.pbe-spn-kjpaw_psl.1.0.0.UPF
    ATOMIC_POSITIONS
    Nb 0.00 0.00 0.00

    K_POINTS {crystal_b}
    8
    0.0000000000    0.0000000000    0.0000000000 20 !G
    0.5000000000    -0.5000000000   0.5000000000 20 !H
    0.0000000000    0.0000000000    0.5000000000 20 !N
    0.0000000000    0.0000000000    0.0000000000 20 !G
    0.2500000000    0.2500000000    0.2500000000 20 !P
    0.0000000000    0.0000000000    0.5000000000 20 !N
    0.2500000000    0.2500000000    0.2500000000 20 !P
    0.5000000000    -0.5000000000   0.5000000000 20 !H

We run this file with the command ``pw.x < bands.scf.in > bands_scf.out``. When it finishes we are ready
to get the bands, now the file in charge is ``bands_scf_pp.in``

.. code-block:: text
    :linenos:

    &BANDS
        prefix = 'nb'
        outdir = './'
        filband = 'nb_bands.dat'
    /

We run this last file with the command ``bands.x < bands_scf_pp.in > bands_scf_pp.out``.

Before going to the next section, it will be important to know how many orbitals should we use on TITAN.
Since the exported Hamiltonians are not in general for 9 atomic orbitals. To know the correct number we 
must open the pseudopotential file to see which orbitals are used. In this example we have in the file a
section that looks like this 

.. code-block:: text
    :linenos:

    Valence configuration:
    nl pn  l   occ       Rcut    Rcut US       E pseu
    4S  1  0  2.00      1.100      1.150    -4.470845
    5S  2  0  2.00      1.100      1.150    -0.344216
    4P  2  1  6.00      1.100      1.500    -2.705887
    5P  3  1  0.00      1.100      1.500    -0.107294
    4D  3  2  3.00      1.300      1.700    -0.340652

Thus, the orbital configuration we will need in TITAN is ``1S 1S 3P 3P 5D``, 13 orbitals in total.

Calculations in PAOFLOW
-----------------------

This part is far less complicated than the previous one. You can either copy everything from the ``1_quantum_espresso``
to ``2_paoflow`` or you can simply add the input file for PAOFLOW in ``1_quantum_espresso``. In the latter there is a 
file called ``nb.xml`` after running the steps in the previous section, copy that file into a *new* file called ``data-file-schema.xml``.
This is the file that PAOFLOW will look for while running. The input file in this case is ``main.py``

.. code-block:: python3
    :linenos:

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

Running this file is trivial, one simply does ``ipython main.py`` and it does everything automatically. At the end there
should be a new folder called ``output`` with all the files generated by PAOFLOW: bands, DOS, and the Hamiltonian.