.. index:: Response functions

******************
Response functions
******************

For ``-> itype = 8``, TITAN calculates frequency-dependent response to electric fields.

.. attention::
    Only the local responses (local spin and orbital moments, torques) are currently calculated.
    Currents are non-local responses in the form of :math:`\delta\langle O_{ij} \rangle`, and this is
    currently not working with the generalized cell.
    The main problem is the relation between Hartree-Fock and RPA responses.
