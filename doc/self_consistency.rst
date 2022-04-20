.. index:: Self-consistency

*********************
Self-consistency only
*********************

The first step on TITAN is always the :doc:`self-consistency<self_consistency>`.
In the special case that ``Un=Um=0`` and ``fixEf``, the system does not need a self-consistency.

If a self-consistency file for the given parameters already exists, it will be automatically read and used as the initial guess.
To skip a new convergence altogether, the parameter ``-> skipsc`` can be set to ``T``.

.. math::

    n^{\text{out}}_{i\mu}(\{n^{\text{in},d},\mathbf{m}^{\text{in},d}\},E_\text{F}) - n^{\text{in}}_{i\mu} &= 0 \\
    \mathbf{m}^{\text{out},d}_i(\{n^{\text{in},d},\mathbf{m}^{\text{in},d}\},E_\text{F}) - \mathbf{m}^{\text{in},d}_i &= 0 \\
    \sum_{i\mu} n^\text{out}_{i\mu}(\{n^\text{in},\mathbf{m}^\text{in}\},E_\text{F}) - N_\text{el} &= 0

The equations are written in this way because it uses the non-linear root finder ``dnsqe`` from `Slatec <https://www.netlib.org/slatec/src/>`_ (which is included, together with its dependencies, as a module in :guilabel:`mod_dnsqe.F90`).


