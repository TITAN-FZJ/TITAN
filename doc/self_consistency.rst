.. index:: Self-consistency

*********************
Self-consistency only
*********************

The first step on TITAN is always the :doc:`self-consistency<self_consistency>`.
In the special case that ``Un=Um=0`` and ``fixEf``, the system does not need a self-consistency.

If a self-consistency file for the given parameters already exists, it will be automatically read and used as the initial guess.
To skip a new convergence altogether, the parameter ``-> skipsc`` can be set to ``T``.




