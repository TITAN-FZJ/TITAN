.. index:: Input Files

***********
Input Files
***********


Main files
==========

The basic files needed by TITAN are the following:

.. index:: input

:guilabel:`input`
-----------------

In the :guilabel:`input` file, the main options and parameters are located.


.. index:: basis

:guilabel:`basis`
-----------------

The :guilabel:`basis` contains the information about the lattice: the lattive parameter, bravais vectors, atom types, and position of the atoms in the unit cell.
This file is written in shte `POSCAR` format.

.. index:: Parameter file

parameter file
--------------

Each element listed in the :guilabel:`basis` file require an parameter file with the same name as the element.


Secondary files
===============

Other files that may be necessary depending on the calculation:

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




