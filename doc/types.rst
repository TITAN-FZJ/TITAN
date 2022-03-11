.. index:: Types of calculations

*********************
Types of calculations
*********************

TITAN can calculate different types of properties of the systems.
The type of calculation is chosen via the parameter ``-> itype`` in the :guilabel:`input` file.
Values vary from 1 to 12, with the following meaning:


.. index:: Self-consistency

Self-consistency only (``itype = 1``)
=====================================

The first step on TITAN is always the :doc:`self-consistency<self_consistency>`.
In the special case that ``Un=Um=0`` and ``fixEf``, the system does not need a self-consistency.

If a self-consistency file for the given parameters already exists, it will be automatically read and used as the initial guess.
To skip a new convergence altogether, the parameter ``-> skipsc`` can be set to ``T``.

.. toctree::
    self_consistency

.. index:: LDOS, Local Density of States

Local Density of States (``itype = 2``)
=======================================


.. index:: LDOS and Coupling

Local Density of States and Coupling (``itype = 3``)
====================================================


.. index:: Band Structure

Band Structure (``itype = 4``)
==============================


.. index:: Fermi Surface

Iso-Energy Quantities (``itype = 5``)
=====================================

Default: Fermi-Surface, can loop


.. index:: Coupling tensor, MEI, Magnetic Exchange Interaction, DMI, Dzyaloshinskii-Moriya Interaction

Coupling tensor (``itype = 6``)
===============================

Includes Magnetic Exchange Interaction (MEI) and Dzyaloshinskii-Moriya Interaction (DMI).


.. index:: Response function

Susceptibility (``itype = 7``)
==============================


.. index:: Response functions

Response functions (``itype = 8``)
==================================


.. index:: Response functions (Magnetic Field)

Response functions x magnetic field (``itype = 9``)
===================================================

Different than the previous one, this type outputs the response functions as a function of the applied static magnetic field.


.. index:: Gilbert damping

Gilbert Damping (``itype = 10``)
================================

Calculation of Gilbert Damping using the slope of the inverse susceptibility or via the Torque correlation models


.. index:: Time propagation

Time Propagation (``itype = 11``)
=================================

Calculation of Gilbert Damping using the slope of the inverse susceptibility or via the Torque correlation models



.. index:: Coupling (Real Space)

Coupling (Real Space) (``itype = 12``)
======================================

Coupling tensor in real space





