.. index:: Response functions (Magnetic Field)

***********************************
Response functions x magnetic field
***********************************

For ``-> itype = 9``, the magnetic-field-dependent response functions to electric fields are calculated
Different than the :doc:`frequency-dependent<response>`, this type outputs the response functions as a function of the applied static magnetic field.

.. attention::
    Only the local responses (local spin and orbital moments, torques) are currently calculated.
    Currents are non-local responses in the form of :math:`\delta\langle O_{ij} \rangle`, and this is
    currently not working with the generalized cell.
    The main problem is the relation between Hartree-Fock and RPA responses.
