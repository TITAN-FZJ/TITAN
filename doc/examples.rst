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