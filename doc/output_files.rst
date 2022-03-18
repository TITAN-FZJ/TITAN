.. index:: Output Files

************
Output Files
************

These are the basic files that TITAN outputs during any type of calculation.

.. tabs::

    .. tab:: ``-> output``

        This file is defined in the :guilabel:`input` file via the parameter ``-> output``.
        It contains the main output of TITAN, with the basic information and timings.
        This file can be used to follow the calculation, using for example:

        .. code-block:: bash

            tail -f output_file


    .. tab:: :guilabel:`parameter.in`

        In this file, all the parameters read from the :guilabel:`input` :ref:`file <input>` are dumped.
        It can be used to check if all the parameters were correctly read, as well as a simpler and smaller way to
        store information on a run (since the :guilabel:`input` file may have more text, commented variables, etc.)


    .. tab:: :guilabel:`Atoms`

        This file contains the generated lattice of the system.
        It is used when ``-> tbmode = 1`` (SK parameters), as for the other case the positions of the atoms are read from the input.
        This file is

        .. code-block:: text

            # Lattice constant:
            a0
            # Lattice vectors:
            a1   |  in the same units as a0 is given
            a2   |  i.e., the values shown here are obtained 
            a3   |  as a1, a2, a3 from basis multiplied by a0
            # Basis Atoms:
            r1   |  position of the atoms in the same unit as a0,
            r2   |  already taking into account the units
            ...  |  given in the basis file
            =========================================================================
            #---------------------------- Neighbor Atoms ----------------------------
            # Atom:            1
            # Stage:           1
            n11_1   |  nearest neighbor positions of atom 1
            n11_2   |  
            ...     |
            # Stage:            2
            n12_1   |  second nearest neighbor positions of atom 1
            n12_2   |  
            ...     |
            ...     |
            # Atom:            2
            # Stage:            1
            n21_1   |  nearest neighbor positions of atom 2
            n21_2   |  
            ...     |
            # Stage:            2
            n22_1   |  second nearest neighbor positions of atom 2
            n22_2   |  
            ...     |
            ...     |

        The number of stages is chosen via the parameter ``-> nn_stages`` in :guilabel:`input`.


Apart from these files, TITAN will generate output files that are put on the ``${TITAN}/results`` folder.
An example of its structure is given below, for a system with 4 atoms in the unit cell (``4Sites``):

.. code-block:: text

    results/
    ├── FSOC                      | Spin-orbit coupling = False
    │   └── selfconsistency
    └── TSOC                      | Spin-orbit coupling = True
        ├── 4Sites                | 4 sites in the unit cell
        │   ├── A                 | Gilbert damping alpha...
        │   │   ├── Slope         | ...from the slope of the susceptibility
        │   │   └── TCM           | ...from the torque correlation methods
        │   ├── Beff              | Effective magnetic fields (when electric field is applied)
        │   ├── BS                | Band Structure
        │   ├── CC                | Charge current (currently not working)
        │   ├── CD                | Charge disturbance
        │   ├── FS                | Fermi Surface / Iso-Energy
        │   ├── HF                | Hartree-Fock susceptibilities
        │   ├── Jij               | Coupling tensor
        │   ├── LC                | Orbital angular momentum current (currently not working)
        │   ├── LD                | Orbital angular momentum disturbance
        │   ├── LDOS              | Local Density of States
        │   ├── RPA               | RPA susceptibilities
        │   ├── SC                | Spin current (currently not working)
        │   ├── SD                | Spin disturbance
        │   ├── SHA               | Spin Hall Angle (currently not working)
        │   ├── SOT               | Spin Orbit Torques
        │   └── time_propagation  | Real-time propagation observables
        └── selfconsistency

It contains first a separation into calculations with (``T``) and without (``F``) SOC.
Then, a folder with the number of sites in the unit cell is used to store all the different types of results.
For a given calculation (i.e., spin-orbit status and number of sites), the complete folder structure is always generated.
In this case, only the calculation with ``-> SOC = T`` was run.

A second folder is located inside the SOC separation, called ``selfconsistency``.
There, the self-consistency files (one for each atom in the unit cell, and also depending on a few parameters) are stored.
TITAN always tries to read the files for the same parameters from these folders to use as an `initial guess`.

.. note::
    When the effective electronic interaction ``Um`` or ``Un`` is non-zero,
    the self-consistency needs an initial density (that is, for the pure tight-binding hamiltonian).
    In this case, the values are calculated for each atom in the unit cell, and files ``initialrho`` are put into ``results/FSOC/selfconsistency``.