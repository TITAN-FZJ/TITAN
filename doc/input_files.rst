.. index:: Input Files

***********
Input Files
***********


Main files
==========

Here we list and explain the basic files that are required in every calculation of TITAN.
Examples of these files can be found :doc:`here<examples>`.

.. index:: input
.. _input:

:guilabel:`input`
-----------------

In the :guilabel:`input` file, the main options and parameters are located.

.. note::
    This file is processed in the subroutine ``get_parameters`` of :guilabel:`mod_io.F90`.

.. index:: basis
.. _basis:

:guilabel:`basis`
-----------------

This file contains the information about the lattice: the lattive parameter, bravais vectors, atom types, and position of the atoms in the unit cell.
This file is written in the `POSCAR` format.

.. note::
    This file is processed in the subroutine ``read_basis`` of :guilabel:`mod_polyBasis.F90`.


The basis file for Fe bulk, for example, is:

.. code-block:: text
    :linenos:

    Fe bulk               ! Name
    5.30                  ! Lattice parameter
    0.5	    0.5	   -0.5   ! Bravais vector a1
    0.5	   -0.5	    0.5   ! Bravais vector a2
    -0.5    0.5	    0.5   ! Bravais vector a3
    Fe                    ! Elements
    1                     ! Number of atoms of each element
    L                     ! Units of position vectors
    0.0   0.0   0.0       ! Position of the atoms

| The first line contains a name, for easy identification.  
| The second is the lattice parameter :math:`a_0`. It can be given in any units (Angstrons or atomic units), which will change the output quantities that depend on length units.  
| From the third to fifth lines are the Bravais vectors :math:`\mathbf{a}_1`, :math:`\mathbf{a}_2`, :math:`\mathbf{a}_3` (the above case is for a `bcc` lattice).  
| The sixth lines contain the different elements that compose the material, with the number of each listed on the seventh line, in the same order.  
| Line number 8 contains the units where the positions are given; it can be given in:

* ``cartesian``, where the positions :math:`x,y,z` in each line are only multiplied by the lattice parameter. This means that the real position of the atoms are:

    :math:`\mathbf{r} = (x \mathbf{\hat{x}}+ y \mathbf{\hat{y}}+ z \mathbf{\hat{z}})a_0`

* ``bravais`` or ``lattice``, where the numbers are multiplied also by the Bravais vectors, such that

    :math:`\mathbf{r} = (x \mathbf{a}_1+ y \mathbf{a}_2+ z \mathbf{a}_3)a_0`

The positions of all the basis atoms must be then listed starting on line 9.

.. important::
    The number of position lines should be the same number as the sum of the numbers in line 7.


.. index:: Parameter file

parameter file
--------------

Each element listed in line 6 of the :ref:`basis file<basis>` require a parameter file with the same name as the element is given.
This file contains basis information such as SOC strength, superconductivity :math:`\lambda`, effective electron-electron interaction :math:`U_n` and :math:`U_m`.
For Slater-Koster parameters (i.e., ``-> tbmode = 1``), it also contains the basis information of the original system as well as the SK two-center integrals.

.. note::
    This file is processed in the subroutine ``readElementFile`` of :guilabel:`TightBinding.F90`.

``-> tbmode = 1``
#################

Below is an example of a parameter file :guilabel:`Fe` for the Fe bulk system.

.. code-block:: text
    :linenos:

    Fe bcc bulk
    5.30
    0.5     0.5     -0.5
    0.5     -0.5    0.5
    -0.5    0.5     0.5
    Fe
    1
    L
    0.0    0.0     0.0
    3                  ! Dimension of the system
    0.7415             ! Original Fermi energy
    0.58  0.37  7.05   ! Original occupations s, p, d
    0.000  0.07353     ! Effective electron-electron interaction Un, Um
    0.000  0.004       ! SOC strength for p and d orbitals
    0.0 0.0 0.0        ! Superconductivity parameters
    3                  ! Number of stages of neighbors
    s on-site	 1.2017709017      1        Fe
    p on-site	 1.8725119829      2        Fe
    t2g on-site	 0.6881678104      3        Fe
    eg on-site	 0.6643740535      4        Fe
    sss 1st nn	-0.1394413859      5        Fe
    pps 1st nn	 0.2681021988      6        Fe
    ppp 1st nn	 0.0297146384      7        Fe
    dds 1st nn	-0.0508569255      8        Fe
    ddp 1st nn	 0.0309574008      9        Fe
    ddd 1st nn	-0.0030320531     10        Fe
    sps 1st nn	 0.1777951121     11        Fe
    sds 1st nn	-0.0678095073     12        Fe
    pds 1st nn	-0.0930757448     13        Fe
    pdp 1st nn	 0.0208929181     14        Fe
    sss 2nd nn	-0.0314096436     15        Fe
    pps 2nd nn	 0.1884829849     16        Fe
    ppp 2nd nn	 0.0390681326     17        Fe
    dds 2nd nn	-0.0312470067     18        Fe
    ddp 2nd nn	 0.0061819027     19        Fe
    ddd 2nd nn	 0.0007075703     20        Fe
    sps 2nd nn	 0.0735426247     21        Fe
    sds 2nd nn	-0.0388437621     22        Fe
    pds 2nd nn	-0.0602805056     23        Fe
    pdp 2nd nn	-0.0038276755     24        Fe
    sss 3rd nn	 0.0181787629     25        Fe
    pps 3rd nn	-0.0444739647     26        Fe
    ppp 3rd nn	 0.0164096598     27        Fe
    dds 3rd nn	 0.0016750902     28        Fe
    ddp 3rd nn	 0.0003651654     29        Fe
    ddd 3rd nn	-0.0005600667     30        Fe
    sps 3rd nn	-0.0256738886     31        Fe
    sds 3rd nn	-0.0003887220     32        Fe
    pds 3rd nn	 0.0006854520     33        Fe
    pdp 3rd nn	-0.0028157220     34        Fe

When using the SK parameters, the parameter file first contains initially the POSCAR file of the original system, from where the parameters were obtained.
In the case above, it is contained between lines 1 and 9, and is the same as the :guilabel:`basis` file, as the system to be investigated is the same as the parameters were obtained.
The parameter file then lists the following quantities:

* Dimension of the system: The value can be different from the ``-> isysdim`` on the :guilabel:`input` file. For example, the parameters can be obtained from Fe bulk (3D) and used in a layered system, where ``-> isysdim = 2``.
* Fermi energy: The value of the Fermi energy from the system where the parameters were obtained.
* `s`, `p` and `d` occupations: The number of electrons on the original system (NOTE: This may be deprecated. It was used as :math:`n_0` in the :math:`U_n (n-n_0)` term of the Hamiltonian).
* Effective electron-electron interaction :math:`U_n`, :math:`U_m`: charge and magnetic parts of the intra-atomic electron-electron interaction, *in the same units as the hopping parameters*
* SOC strength for p and d orbitals: :math:`\lambda_{SOC}` for `p` and `d` orbitals (the spherical `s` orbital does not contribute to the SOC term)
* Superconductivity :math:`\lambda`: This can be given in different ways, in the following order or priority:

  * A **single** value for all the orbitals;
  * **One value per orbital** of the given atom type;
  * **Three values**, one per general orbital type (`s`, `p` and `d`);
  * **Nine values**, one per specific orbital type (`s`, `px`, `py`, `pz`, `dxy`, `dyz`, `dzx`, `dx2`, `dz2`).

* Number of neighbor stages: This determines how many parameters will be read on the lines below.

After these parameters are given, for each neightbor stage, TITAN reads 10 values of each two center integrals: `sss`, `pps`, `ppp`, `dds`, `ddp`, `ddd`, `sps`, `sds`, `pds`, `pdp`, one per line.
These values can be obtained, for example in the Handbook of Papaconstantopoulos [Papa]_.
They are given in this format to make it easy for copying and pasting from the book.

.. note::
    There used to be a `database website <http://esd.cos.gmu.edu/tb/>`_ containing the SK parameters from the Handbook of Papaconstantopoulos, but the website has changed. 
    `Another website <http://cmasc.gmu.edu/esd/>`_ seems to have the original database, but it crashes when requesting a parameter.


``-> tbmode = 2``
#################

When using hamiltonians from DFT (either PAOFLOW or Wannier), the basis and the SK parameters are not needed since the hamiltonian is directly read (and it doesn't need to be built from the parameters or the hoppings for given neighbors to be constructed).
It is also assumed that the system is the same as the one to be investigated in TITAN.

Below is an example of a parameter file :guilabel:`Co` for the Co monolayer system.

.. code-block:: text
    :linenos:

    Co             ! Name
    0.0            ! Original Fermi energy
    0.000 1.0      ! Effective electron-electron interaction Un, Um
    0.000 0.085    ! SOC strength for p and d orbitals
    0.0 0.0 0.0    ! Superconductivity parameters

The parameter file is simplified in this case, and only contains:

* Name of the system
* Fermi energy: The value of the Fermi energy from the system where the parameters were obtained.
* Effective electron-electron interaction :math:`U_n`, :math:`U_m`: charge and magnetic parts of the intra-atomic electron-electron interaction, *in the same units as the hopping parameters*
* SOC strength for p and d orbitals: :math:`\lambda_{SOC}` for `p` and `d` orbitals (the spherical `s` orbital does not contribute to the SOC term)
* Superconductivity :math:`\lambda`: This can be given in different ways, in the following order or priority:

  * A **single** value for all the orbitals;
  * **One value per orbital** of the given atom type;
  * **Three values**, one per general orbital type (`s`, `p` and `d`);
  * **Nine values**, one per specific orbital type (`s`, `px`, `py`, `pz`, `dxy`, `dyz`, `dzx`, `dx2`, `dz2`).



Secondary files
===============

Other files that may be necessary depending on the calculation:


.. index:: hamiltonian_file

Hamiltonian file
----------------

The Hamiltonian obtained from DFT when ``-> tbmode = 2``


.. index:: kbands

:guilabel:`kbands`
------------------

The file :guilabel:`kbands` includes the definition of k-points to be used in the band structure (``itype=4``) and q-dependent susceptibility calculation (``itype=7``).
The units is defined on the :guilabel:`input` file via the parameter ``-> qbasis``.



.. index:: initialmag

:guilabel:`initialmag`
----------------------

TITAN initial magnetic moments for the atoms in the unit cell is `2` along the `z-`direction.
This can be changed by adding ``-> magbasis`` in the :guilabel:`input` file (possible values are ``cartesian`` or ``spherical``), and adding the values in a file called :guilabel:`initialmag`.
It must contain 3 values per line, with the number of lines given by the number of atoms in the unit cell (following the same order as the positions).




