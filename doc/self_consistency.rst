.. index:: Self-consistency

*********************
Self-consistency only
*********************

The first step on TITAN is always the self-consistency.

.. math::

  n^{\text{out}}_{i\mu}(\{n^{\text{in},d},\mathbf{m}^{\text{in},d}\},E_\text{F}) - n^{\text{in}}_{i\mu} &= 0 \\
  \mathbf{m}^{\text{out},d}_i(\{n^{\text{in},d},\mathbf{m}^{\text{in},d}\},E_\text{F}) - \mathbf{m}^{\text{in},d}_i &= 0 \\
  \sum_{i\mu} n^\text{out}_{i\mu}(\{n^\text{in},\mathbf{m}^\text{in}\},E_\text{F}) - N_\text{el} &= 0

The equations are written in this way because it uses the non-linear root finder ``dnsqe`` from `Slatec <https://www.netlib.org/slatec/src/>`_ (which is included, together with its dependencies, as a module in :guilabel:`mod_dnsqe.F90`).
They are processed in the following way: 

* Initial values :math:`n^{\text{in}}_\mu, \mathbf{m}^{\text{in}}, \Delta_\text{SC}^{\text{in}}` are given (the charge must be orbital resolved due to the form of the :doc:`Hamiltonian <hamiltonian>`);
* With those values, the new :doc:`expectation values <expectation_values>` are obtained;
* The sum of the occupations are calculated (when ``fixEf`` is not used);
* If the output values are not the same as the input (within a tolerance that may be changed with the :ref:`input parameter <input>` ``-> magtol``), ``dnsqe`` will generate new "initial" values :math:`n^{\text{in}}_\mu, \mathbf{m}^{\text{in}}, \Delta_\text{SC}^{\text{in}}` and the process restarts.

The last equation represents *total charge neutrality*. 
This equation determines the Fermi level :math:`E_\text{F}` such that the total initial number of electrons :math:`N_\text{el}` is kept constant.

.. note::
  In the special case that ``Un=Um=0`` and ``fixEf``, the system does not need to perform the self-consistency.
  It will only calculate the ground-state expectation values :math:`n, \mathbf{m}, \Delta_\text{SC}` for every orbital and site in the unit cell, and write to the output self-consistency files, together with the Fermi energy :math:`E_\text{F}`, to be read when running the code with the same parameters again.

If a self-consistency file for the given parameters already exists, it will be automatically read and used as the initial guess.

.. tip::
  To skip a new convergence altogether, the parameter ``-> skipsc`` can be set to ``T``.
