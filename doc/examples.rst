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

It is necessary to get the pseudopotential file for QE, it can be found in the following URL: 
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

Calculations in TITAN
---------------------

I this section many things are skipped, and we focus only on the differences or additions that one needs to make
to run this example successfully. First we need to look for the file ``z2pack_hamiltonian.dat`` in the ``output`` folder
from the last section, and copy that file into the folder ``3_titan`` with the name of ``Nb_PAOFLOW``. That file should
look like this 

.. code-block:: text
    :linenos:

    PAOFLOW Generated
       13
     2197
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ...

and so on. The ``basis`` file for bulk Nb has no changes, and it is as follows

.. code-block:: text
    :linenos:

    Nb_PAOFLOW
    6.23610
    0.5     0.5     0.5
    -0.5    0.5    0.5
    -0.5    -0.5     0.5
    Nb
    1
    L
    0.0    0.0     0.0 Nb

For the material card ``Nb`` the most important part is the list of orbitals 
``o: s s px py pz px py pz dz2 dzx dyz dx2 dxy``, as mentioned earlier the pseudopotential we used for Nb on 
QE contains 13 states, 2 semicore states s and p, and the 9 usual valence orbitals, thus we list the orbitals 
accordingly on TITAN. The rest of the interactions are set to zero since they must already be included in the 
Hamiltonian.

.. code-block:: text
    :linenos:

    Nb
    o: s s px py pz px py pz dz2 dzx dyz dx2 dxy
    0.0         !Fermi level
    0.0 0.0     !Un and Um in eV
    0.0 0.0     !lambda_SOI for p and d orbitals
    0.0 0.0 0.0 !superconducting lambda

The ``kbands`` file also remains unchanged

.. code-block:: text
    :linenos:

    G       0.0     0.0     0.0
    H       -0.5    0.5     0.5
    N       0.0     0.5     0.0
    P       0.25    0.25   0.25

Finally the most important part is the ``input`` file. In this example we calculate the band structure, hence we 
pick ``-> itype = 4``, and since we are using the inputs from DFT we choose ``-> tbmode = 2``. We do not want 
to run self consistent calculations in this case, otherwise some bands would move and everything ends up in a 
disaster. To avoid this we fix the Fermi level by adding the runoption ``fixEf`` to ``-> Options =`` and we fix
the Fermi level at zero as it comes from PAOFLOW with ``-> Ef = 0.d0``. All the fields must be off, hence we 
make ``-> FIELD = F``, ``-> SOC = F``, and ``-> superCond = F``. Finally, as mentioned earlier, we request TITAN 
*not to perform a self-consistent calculation*; we manage this by skipping the self-consistency with ``-> skipsc = T``.

.. code-block:: text
    :linenos:

    !*******************************************************************************!
    !                           CHOOSE WHAT TO CALCULATE                            !
    -> itype = 4
    ! Slater-Koster (1), DFT (2)
    -> tbmode = 2
    !===============================================================================!
    !                                OUTPUT FILE                                    !
    -> output = output/bands
    !===============================================================================!
    !                         OPTIONAL RUNNING VARIABLES                            !
    -> Options = createfolders positions ontheflysc eigenstates nojac fixEf
    -> Ef = 0.d0
    !===============================================================================!
    !                               SYSTEM VARIABLES                                !
    ! Lattice and surface direction:                                                !
    -> nn_stages = 1
    ! Small imaginary part of the Green function:
    -> eta = 1e-3
    !===============================================================================!
    !
    -> sysdim = 3
    !===============================================================================!
    !                          STATIC MAGNETIC ZEEMAN FIELD                         !
    !                         (in the spin reference system)                        !
    -> FIELD = F  ! Magnetic field on (T) or off (F) (choose one form below)
    ! Spherical coordinates (default, angles in units of pi)                        !
    -> hwa   = 5.788381e-5
    -> hwt   = 0.0E+00
    -> hwp   = 0.0E+00
    ! Cartesian coordinates                                                         !
    -> hwx = 0.0E+00
    -> hwy = 0.0E+00
    -> hwz = 0.0E+00
    !===============================================================================!
    !                             SPIN ORBIT COUPLING                               !
    -> SOC = F
    -> socscale = 1.00d0
    !===============================================================================!
    !                      REAL-TIME PROPAGATION  (itype=11)                        !
    ! Intensity of transverse magnetic field:                                       !
    -> hw1 = 5.788381e-6
    ! Frequency of transverse magnetic field:                                       !
    -> hw = 11.576762e-5
    ! Total time:                                                                   !
    -> integration_time = 5.e6
    ! Step size:                                                                    !
    -> step = 120
    !===============================================================================!
    !                   DIRECTION OF IN-PLANE ELECTRIC FIELD                        !
    -> ebasis = cartesian
    -> dirEfield = 1.0 0.0 0.0
    !===============================================================================!
    !                        IN-PLANE CURRENTS TO CALCULATE                         !
    -> n0sc1 = 1   ! First neighbor
    -> n0sc2 = 6   ! Last neighbor
    !
    !                    CURRENT AND DISTURBANCE RENORMALIZATION                    !
    !          (only used when currents are calculated - itype=7 and 8)             !
    -> renorm = T        ! Turn renormalization on (T) or off (F)
    -> renormnb = 1      ! Reference neighbor (where the charge current will be 1)
    !===============================================================================!
    !                           INTEGRATION VARIABLES                               !
    ! Approximate number of k-points: (nkpt > 0)                                    !
    -> nkpt = 10000
    ! Number of parts to divide energy integral in complex plane:                   !
    -> parts = 2
    -> n1gl = 128
    !  Number of parts to divide energy integral in real axis:                      !
    -> parts3 = 1
    -> n3gl = 128
    !===============================================================================!
    !                       SUPERCONDUCTING VARIABLES                               !
    -> superCond = F
    !===============================================================================!
    !                              PLOTTING VARIABLES                               !
    ! Energy range and number of points:                                            !
    -> emin =  -30.0
    -> emax =  15.0
    -> nEner = 800
    !-> skip_steps = 0
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    ! Wave vector path and number of points:                                        !
    -> band = G H N G P N P H
    -> nQvec = 500
    -> qbasis = bravais
    !===============================================================================!
    !                              SELF-CONSISTENCY                                 !
    !        (Finding center of the bands for the given number of particles)        !
    -> skipsc = T ! skip self-consistency calculation when it finds previous results
    ! File to use as starting point:
    !-> scfile = CoSCTSOC.dat
    !===============================================================================!
    !-> fermi_layer = 1 ! Which layer should the Fermi-Level be taken from
    !*******************************************************************************!

This should be enough to get the band structure with TITAN from PAOFLOW inputs. This procedure for getting the inputs
for TITAN is the same, regardless of the system to use. 