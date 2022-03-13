.. index:: Introduction
.. sectionauthor:: Filipe Guimarães <f.guimaraes@fz-juelich.de>

************
Introduction
************

TITAN is a simulation code that focus on Time-dependent transport and angular momentum properties of nanostructures.
It is based on a multi-orbital tight-binding model, where the hamiltonian that can be constructed via Slater-Koster parameters [SK]_ [Papa]_ or via hamiltonians from codes such as `PAOFLOW <http://aflowlib.org/src/paoflow/>`_ or `Wannier90 <http://www.wannier.org>`_.
The electron-electron interaction can be taken into account (apart from the single-particle contribution inside the tight-binding hamiltonian) via
an effective Hubbard-like hamiltonian [Hubbard]_.
The spin-orbit interaction may also be included via a :math:`\lambda \mathbf{L}\cdot \mathbf{S}` term.

.. attention::
    The basis for the :math:`\mathbf{L}` matrix is assumed to be the local atomic orbitals (i.e., *s*, *px*, *py*, *pz*, *dxy*, *dyz*, *dzx*, *dx2* and *dz2*).


A static magnetic field :math:`\mathbf{B}` can be applied in any direction via a term :math:`g \mu_B \mathbf{B}\cdot \mathbf{S}` in the hamiltonian.
Response functions are obtained in Kubo's linear response approach [Kubo]_.
