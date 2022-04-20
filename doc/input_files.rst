.. index:: Input Files

***********
Input Files
***********


Main files
==========

Here we list and explain the basic files that are required in every calculation of TITAN.
Examples of these files can be found :doc:`here<examples>`.

.. tip::
    For these files, comments can be added with ``!`` (i.e., text on the right of it is not processed). 


.. index:: input
.. _input:

:guilabel:`input`
-----------------

In the :guilabel:`input` file, the main options and parameters are located.

.. note::
  This file is processed in the subroutine ``get_parameters`` of :guilabel:`mod_io.F90`.

.. attention::
  Inputs should be added one per line, and must start with ``->`` (see subroutine ``read_file`` in :guilabel:`mod_input.F90`).

The following is an exhaustive list of all possible input parameters that may be read from :guilabel:`input`,
including if they are required (and for which ``itype``) or optional.

* ``-> tbmode``: Which type of tight-binding hamiltonian to use:

  1. Slater-Koster Hamiltonian: Uses the SK parameters read from the :ref:`parameter files<parameter_file>` to build the Hamiltonian. Hopping between elements are built by default with a geometric average :math:`\bar{t}_{12} = \text{sgn}(t_1 + t_2) \sqrt{|t_1||t_2|}` (can be changed with the ``simplemix`` run option, see below). The lattice is generated using the information on the :guilabel:`basis` :ref:`file <basis>` and can be controlled with the variables:

    * ``-> nn_stages`` (default: :math:`2`): Number of nearest neighbor shells/stages to consider
    * ``-> relTol`` (default: :math:`0.05`): Relative tolerance to group neighboring distances in the same shell/stage.

  2. DFT Hamiltonian: Hamiltonian read from direct DFT output via programs such as `PAOFLOW <http://aflowlib.org/src/paoflow/>`_ or `Wannier90 <http://www.wannier.org>`_. The :ref:`Hamiltonian file <hamiltonian_file>` contains information on the amount of orbitals, neighboring shells and the hoppings.

  .. tip::
    The string ``_dfttype=S`` is added to the self-consistency files for ``-> tbmode=1`` and ``_dfttype=D`` for ``-> tbmode=2``.

* ``-> itype`` (**required**): Which calculation will be executed. More information can be found on the :doc:`Types of calculations<types>` page.
* ``-> nkpt`` (**required**): Number of k-points in the Brillouin Zone (see note below about calculations with Green functions). It can be given as (in order of priority):

  * (sysdim) numbers, that will be multiplied to get the total number of kpoints;
  * One number (nkpt), that will be splited into :math:`(\text{nkpt})^{\frac{1}{d}}` per dimension.

  .. attention::
    Some quantities are calculated by integrating Green functions in the imaginary complex plane. In that case,
    the larger the imaginary part :math:`y` of the energy argument :math:`E=E_\text{F}+\mathrm{i}y`, the smoother are the functions to be integrated. Therefore, the 
    BZ mesh is adapted depending on :math:`y` (that we call adaptive mesh) as :math:`\frac{\text{nkpt}}{(\frac{y}{y_0})^\sqrt{d}}`, with :math:`y_0` being the imaginary part of the 
    energy value closes to the real axis (see ``-> eta`` below).

* ``-> minimumBZmesh`` (default: :math:`1000``): Minimum number of kpoints on the smallest BZ mesh (see explanation above). If larger than ``nkpt``, this value it used.
* ``-> eta`` (**required**): Small imaginary part (broadening) of the energy argument :math:`E=E_\text{F}+\mathrm{i}\eta` (in calculations using Green function) or the energy equivalent of the temperature (in calculations using eigenstates) :math:`\beta = \frac{1}{k_\text{B}T} \equiv \frac{1}{\pi\eta}`, where :math:`k_\text{B}` is the Boltzmann constant and the factor :math:`\frac{1}{\pi}` was added to cancel the existing one on the first Matsubara pole (:math:`\hbar\omega _{n}={\frac {(2n+1)\pi }{\beta }}`, for :math:`n=0`) and make it closer to the limit of the integration of the Green functions (usually from :math:`\eta` to :math:`\infty` in the imaginary energy axis).
* ``-> etap`` (only for response calculations, ``-> itype=7-9``, default: same as ``eta``): Same as above, but as the responses are composed of a product of Green functions, there's the possibility to use different broadenings for each. May be useful to make the meaning of the broadening the same as the one added in the susceptibilities calculated with eigenstates (in the discussion of Gilbert damping, for example).
* ``-> sysdim`` (default: :math:`3`): Dimension :math:`d` of the system to be studied (i.e., the one given in :guilabel:`basis`).
* ``-> orbitals`` (default: `s px py pz dxy dyz dzx dx2 dz2`): Which orbitals to use in the hamiltonian for all the atoms in the basis. The default are the 1 s, 3 p and 5 d. For type- or atom-dependent orbitals, use the option in the :ref:`parameter file <parameter_file>`.
* ``-> fermi_layer`` (default: :math:`1``): which site/layer to use the (initial) Fermi energy from.
* ``-> Ef`` (default: `none`): Value of Fermi energy to overwrite the original one. By default, the Fermi level is then modified to keep the total number of electrons constant. To fix the Fermi level, use the run option ``fixEf`` (see below).
* ``-> output`` (**required**): Main output file of TITAN, where information about all the steps are logged (see more info on :doc:`Output files<output_files>`).

  .. tip::
    If the string ``#nkpt`` or ``#eta`` is used on the output filename, it is substituted to the respective parameters.
    The file can also be in a subfolder (e.g., ``-> output = output/outfile``), which will be created during runtime.
    Existing files are overwritten.

* ``-> suffix`` (default: `none`): A string to be added in all results files. This is useful to avoid files to be overwritten (when the file does not already contain changes in a given parameter).
* ``-> SOC`` (default: ``T``): Option to turn the spin-orbit Coupling on ``T`` or off ``F``. This option is used in the folder structure as ``TSOC`` or ``FSOC``, respectively (see :doc:`Output Files<output_files>`). There is also a run option for linear SOC, whose folder is ``LSOC`` (see below).
* ``-> socscale`` (only valid for ``-> SOC = T`` ): A real number to scale the strength of the SOC :math:`\lambda_\text{SOC}` in all the atoms.
* ``-> skipsc`` (default: ``F``): Option to skip the self-consistency when it was already done (i.e., the files already exist).
* ``-> maxiter`` (default: :math:`99999999``): Maximum number of steps in the self-consistency.
* ``-> scfile`` (default: `none`): Filename to use as initial input for the self-consistency. Useful to use an already converged calculation as input for other parameters (instead of copying and renaming the file).
* ``-> addelectrons`` (default: :math:`0.0`): Number of electrons to add (positive) or remove (negative) from the system. It does not affect the calculations when the runoption ``fixEf`` is used.
* ``-> Options`` (optional): This parameter can have multiple entries (separated with white spaces) that trigger different behaviour of the code. These are all the recognized keywords - note that some are not implemented, others need to be fixed and yet others are deprecated:

  * ``ry2ev``: Transform energy units from Rydbert to electron-Volt;
  * ``tesla``: Flag to indicate the program that the magnetic field parameters should be read as in units of Tesla;
  * ``verbose``: Flag to indicate the program to log more information (not used in the code);
  * ``debug``: This flag can be used to run codes in the ``if(ldebug)`` region in :guilabel:`main.F90`;
  * ``addresults``: Avoid creating (and overwriting) results output files. Files must already exist.
  * ``createfiles``: Only create files with headers (useful when planning to run calculations with ``addresults`` afterwards).
  * ``createfolders`` (**recommended**): Create :ref:`all the sub-folders structure <results-structure>` for the results. Doesn't affect if it already exists (nothing is overwritten).
  * ``kpoints``: Creates an :doc:`output file<output_files>` :guilabel:`kpoints` with all the kpoints, that can be plotted, for example, to check if the BZ is correct.
  * ``positions``: Creates an :doc:`output file<output_files>` :guilabel:`Atoms` containing, for each atom in the basis, its neighboring atoms (types and positions) separated by stages/shells.
  * ``lineargfsoc`` (not implemented in all types of calculations): Includes SOC in a linear approximation within the Green functions, i.e., :math:`G = G_0+G_0H_\text{SOC}G_0`, where :math:`G` is the Green function with SOC and :math:`G_0` without.
  * ``linearsoc`` (should work for at least ``-> itype=7-9``): Calculate responses up to linear term in :math:`\lambda_\text{SOC}`.

  .. note::
    When ``lineargfsoc`` or ``linearsoc`` is used, the string ``_linearsoc`` is appended to the results files.
    Also, in the case of ``linearsoc``, the SOC folder where the :doc:`output files<output_files>` are put is called ``LSOC``.

  * ``nojac`` (always active when using ``eigenstates``): Does not use analytical jacobian subroutine, calculates it numerically.
  * ``hfresponses``: Calculate only HF responses (no RPA ones).
  * ``ontheflysc`` (**recommended**): Write every self-consistency step into files (that can be restated later with ``-> skipsc = F``).
  * ``rotatemag``: Rotate the magnetization to the direction of the external magnetic field before the self-consistency (useful when ``-> SOC = F``, as the final magnetization should point in this direction).
  * ``nolb``: Do not add the orbital Zeeman term :math:`H_\text{ZL} = g_L\mu_B\sum_i \mathbf{L}_i\cdot\mathbf{B}`. 
  * ``nodiag`` (``-> itype=7-9``): Do not diagonalize the susceptibility matrix (which is done by default). Useful when the susceptibility matrix become singular.
  * ``writeonscreen`` (``-> itype=8-9``): Also writes some results on the output file (useful for some checks).
  * ``sortfiles``: Only sort files (no calculation is done). Useful when the calculation was done in parallel and finished before sorting the files (therefore, the results may be in random order).
  * ``simplemix`` (``-> tbmode=1``): Calculate hoppings between elements as simple average between each of the values, i.e., :math:`\bar{t}_{12} = \frac{t_1 + t_2}{2}`, instead of the default geometric average :math:`\bar{t}_{12} = \text{sgn}(t_1 + t_2) \sqrt{|t_1||t_2|}`.
  * ``fixEf``: Do not include the total charge neutrality equation :math:`\sum_i n_i(\{n,\mathbf{m}\},E_\text{F}) = N` on the self-consistency. This means the Fermi energy is fixed. Useful to use in conjunction to the ``-> Ef`` input to define the Fermi energy.
  * ``eigenstates``: Calculate the expectation values using the eigenstates and the Fermi distribution :math:`f(E)`, i.e., :math:`\langle \mathcal{O}\rangle = \sum_n f(E_\text{n})\langle n | \mathcal{O} | n \rangle`.

    .. attention::
      The run option ``eigenstates`` is always activated for the real time propgation (``-> itype = 11``), as the implementation relies on the evolution of the coefficients :math:`c_n(t)`.

  * ``printfieldonly`` (``-> tbmode=11``): Use this option to print the time-dependent magnetic and/or electric field that is applied on real time propagations. This is useful to check the field before running the full calculation.

    .. note:: The following ``-> Options`` are either not working, not implemented or obsolete:\

      * ``sha`` (``-> itype=8-9``): Calculate the Spin Hall Angle by doing: :math:`\Theta_\text{SH} = \frac{I^\text{S}_{\perp}}{I^\text{C}_{\|}}`. Must be fixed after current calculations are working (see :doc:`Known Issues<known_issues>`).
      * ``lgtv``: Only calculate total longitudinal and transverse currents. Must be fixed after current calculations are working (see :doc:`Known Issues<known_issues>`).
      * ``checkjac``: Option to check if jacobian is correctly calculated (comparing numerical with analytical implementation). Although the subroutine ``check_jacobian`` is in place in :guilabel:`mod_self_consistency.F90`, it is not called anywhere.
      * ``forceoccupation`` (not implemented yet): Option to input the occupations and find the center of the bands, instead of doing total charge neutrality.

* ``-> n0sc1``, ``-> n0sc2`` (only for current calculations, ``-> itype=8-9``, required but may be obsolete): Initial and final neighbor to calculate the current. This option is not working, since the currents are not being calculated (see :doc:`Known Issues<known_issues>`), but it may also be obsolete since the neighbors are now defined differently. Must be checked when reactivating the current.
* ``-> magtol`` (default: :math:`10^{-12}`): Tolerance to use in the self-consistency.
* ``-> magbasis`` (required to read :guilabel:`inputmag`, see :ref:`below <initialmag>`): This option specifies in which basis to interpret the coordinates of the initial magnetization. Can be ``cartesian``, ``spherical`` or ``bravais`` (maybe ``neighbor`` could be implemented later on).
* ``-> ebasis`` (only valid for ``-> itype=8-9``, default: ``spherical``): Basis to interpret the ``dirEfield`` coordinates. Can be ``cartesian``, ``spherical`` or ``bravais`` (maybe ``neighbor`` could be implemented later on).
* ``-> dirEfield`` (only valid for ``-> itype=8-9``, **required**): Coordinates for the direction of the electric field (Only the direction matters - if it is not a unit vector, it is transformed to one). How many coordinates depends on the ``-> ebasis``:

  * ``cartesian``: `x`, `y` and `z` coordinates;
  * ``spherical``: :math:`\theta` and :math:`\phi` coordinates (in degrees);
  * ``bravais``: :math:`n_1, n_2, n_3`, such that :math:`\mathbf{E} = n_1 \mathbf{a}_1+ n_2 \mathbf{a}_2+ n_3 \mathbf{a}_3`

* ``-> FIELD`` (default: ``F``): Activate ``T`` or deactivate ``F`` static magnetic field. When active, the field can be given as:

  * Spherical coordinates: ``-> hwa``, ``-> hwt`` and ``-> hwp`` (the later two in degrees). Each of these components accept two forms of input:
    * A single number, for a fixed value of the respective magnetic field coordinate to be applied;
    * Three numbers: initial, final, and number of points :math:`\text{npts}`. This will create a loop from the initial to the final values (including them), with :math:`\text{npts}+1` points.
  * Cartesian coordinates: Although not recommended, the code also accepts cartesian coordinates via ``-> hwx``, ``-> hwy`` and ``-> hwz``. To use this, the variable ``-> hwa`` must not be set. **The cartesian components cannot be looped.**

  .. tip::
    When looping over magnetic fields, the variable ``-> skip_steps_hw`` can be used to skip a given number of steps in the beginning of the calculation. This option is useful to skip points that were already calculated.

    The magnetic field can also be modified by other parameters. They are:  

      * ``-> hwscale``: Scale the magnetic field by a real value;
      * ``-> hwtrotate``: Rotate the angle :math:`\theta` (in degrees).
      * ``-> hwprotate``: Rotate the angle :math:`\phi` (in degrees).

    These parameters receive one value per atom in the unit cell (Fortran also accepts repetitions, e.g., ``3*0.5`` for 3 times the factor ``0.5``).

* ``-> constraining_field``: Activate calculations with contraining fields (i.e., applied fields that results in a given magnetization). When activated, it requires ``-> magbasis`` and the :guilabel:`initialmag` with the desired magnetization directions for each site in the unit cell. It also requires:

  * ``-> constr_type``: 

    1. Transverse constraining field: Contrain only the direction of the magnetization (i.e., a field is applied to keep the direction of :math:`\mathbf{m}`` fixed, not its length). An extra parameter ``-> cmix`` can be given (default: :math:`10^{-2}`), to set the mixing parameter between the steps (lower values take longer to converge, but high values may not converge at all);
    2. Full constraning field: The magnetic moment is completely fixed by the application of a constraining field.

* ``-> superCond`` (default: ``F``): Activate ``T`` or deactivate ``F`` superconductivity. When active, the hamiltonian has blocks for electrons and holes. The coupling is defined by the parameter :math:`\lambda_\text{SC}` that is defined in the :ref:`parameter file<parameter_file>`.

  .. tip::
    When ``-> superCond = T`` but :math:`\lambda_\text{SC} = 0` in the parameter file, the results must be the same as ``-> superCond = F``.

Energy Integration:

* ``-> parts``: Number of parts to split the imaginary part from :math:`\eta` to :math:`\infty` of the :doc:`energy integration </energy_integration>`. The parts are splitted exponentially after the transformation to an integral from `0` to `1` (e.g., the calculation is splitted in the intervals `(0,0.01)`, `(0.01,0.1)` and `(0.1,1.0)` for ``-> parts=3``). Each part is integrated using Gauss-Legendre quadrature points.
* ``-> n1gl``: Number of Gauss-Legendre points in each part (defined by ``-> parts``) of the imaginary-axis energy integration.
* ``-> parts3`` (``-> itype=7-9``): Number of parts to split the real part of the :doc:`energy integration </energy_integration>`. Different than the imaginary part, here the integration is splitted in equal pieces (i.e., linear split).
* ``-> n3gl`` (``-> itype=7-9``): Number of Gauss-Legendre points in each part (defined by ``-> parts3``) of the real-axis energy integration.

  .. tip::
    The string ``_parts=X_parts3=Y`` is added to the file name of the results files for the responses (``-> itype=7-9``).

Energy loop:

* ``-> emin`` (``-> itype=2-3,7-9``, required): Starting energy value for energy loops.
* ``-> emax`` (``-> itype=2-3,7-9``, required): Ending energy value for energy loops.
* ``-> nEner`` (``-> itype=2-3,7-9``, required): Number of energy points in the loop (the real number of calculated points is :math:`\text{nEner}+1` to include the edges.)
* ``-> skip_steps`` (default: :math:`0`): Number of energy points to skip on the beginning of the loop (useful to skip points that were already calculated).

Wave-vector loop:

* ``-> band`` (**required** for ``-> itype=4``, default: ``G``(:math:`\Gamma`) for ``-> itype=6-9`` ): Wave-vector points defined in the :guilabel:`kbands` :ref:`file <kbands>` to make a path in the reciprocal space.
* ``-> nQvec`` (``-> itype=2-3,7-9``, default: :math:`0`): Number of Q vectors in the path defined above. The real number of calculated points is :math:`\text{nQvec}+1` to include the edges. When only one letter (point) is given in ``-> band``, only one point is used in the loop.
* ``-> qbasis`` (default: ``bravais``): Which basis to interpret the coordinates in the :guilabel:`kbands` file. Can be ``cartesian``, ``spherical`` or ``bravais`` (which is actually in units of the reciprocal vectors).

* ``-> parField`` (default: :math:`1`): Number of fields to group together in the MPI field communicator. It is (supposedly) used when looping over magnetic fields. Must be tested.
* ``-> parFreq`` (default: :math:`1`): Number of frequencies to group together in the MPI frequency communicator. It is (supposedly) used when looping over frequencies. Must be tested.

* ``-> renorm`` (default: ``F``): Option to renormalize the currents dividing them by the current flowing to one of the neighbors (chosen via ``-> renormnb``). Currently not working (see :doc:`Known Issues<known_issues>`) and probably obsolete.



Specific for some ``itypes``:

For ``-> itype=5``:

* ``fs_energy``: Energy to calculate the iso-energy observables. Can be given as:
  * A single number, the energy to calculate the observables;
  * Three numbers: initial energy :math:`E_i`, final energy :math:`E_f`, and number of points :math:`\text{npts}`. This will create a loop from :math:`E_i` to :math:`E_f` (including those energies), with :math:`\text{npts}+1` points.

For ``-> itype=11``:

For more information on the method and the variables defined here, see :doc:`real time propagation methods<real_time_propagation>`.

* ``-> integration_time`` (**required**): Maximum time.
* ``-> step`` (default: :math:`\text{integration_time}/10^{4}`): Initial time step :math:`\Delta t`. Note that the step will be adjusted at each calculation to avoid large errors.
* ``-> sc_tol`` (default: :math:`10^{-2}`): Time-propagation self-consistency tolerance.
* ``-> abs_tol`` (default: :math:`10^{-3}`): Absolute tolerance for...
* ``-> rel_tol`` (default: :math:`10^{-3}`): Relative tolerance for...
* ``-> safe_factor`` (default: :math:`0.9`): Safety factor for...

* ``-> electric`` (default: ``F``): Turn on or off the application of a time-dependent electric field. Default is an oscillatory electric field.
  
  * ``-> pulse_e`` (default: ``F``): Parameter to define if the electric field is pulsed (``T``) or oscillatory (``F``).

    If ``-> pulse_e = F`` (oscillatory field), the following quantities are required:

    * ``-> hE_0`` (**required**): Amplitude of the oscillation :math:`\hbar E_0`.
    * ``-> hw_e`` (**required**): Frequency/energy of the oscillation :math:`\hbar \omega`.
    * Polarization (see box :ref:`Defining the polarization <polarization>` below).

    If ``-> pulse_e = T``:

    * ``-> npulse_e`` (default: :math:`1`): Number of pulses to be applied.

    Each pulse is defined by a few quantities (which should contain one value per ``npulse_e``):

    * ``-> hE_0`` (**required**): Intensity of the pulse :math:`\hbar E_0`.
    * ``-> hw_e`` (**required**): Internal frequency/energy of the pulse :math:`\hbar \omega`.
    * ``-> tau_e`` (**required**): Width of the pulse :math:`\tau_\text{E}`.
    * ``-> delay_e`` (default: :math:`t=\tau_\text{E}/2``): By default, each pulse starts when the previous ends. This can be changed with this variable (if given, must be one value per pulse).
    * Polarization (see box :ref:`Defining the polarization <polarization>` below).

* ``-> magnetic`` (default: ``F``): Turn on or off the application of a time-dependent magnetic field. Default is an oscillatory magnetic field.
  
  * ``-> pulse_m`` (default: ``F``): Parameter to define if the magnetic field is pulsed (``T``) or oscillatory (``F``).

    If ``-> pulse_m = F`` (oscillatory field), the following quantities are required:

    * ``-> hw1_m`` (**required**): Amplitude of the oscillation :math:`\hbar omega_1`.
    * ``-> hw_e`` (**required**): Frequency/energy of the oscillation :math:`\hbar \omega`.
    * Polarization (see box :ref:`Defining the polarization <polarization>` below).

    If ``-> pulse_m = T``:

    * ``-> npulse_m`` (default: :math:`1`): Number of pulses to be applied.

    Each pulse is defined by a few quantities (which should contain one value per ``npulse_m``):

    * ``-> hw1_m`` (**required**): Intensity of the pulse :math:`\hbar \omega_1`.
    * ``-> hw_m`` (**required**): Internal frequency/energy of the pulse :math:`\hbar \omega`.
    * ``-> tau_m`` (**required**): Width of the pulse :math:`\tau_\text{M}`.
    * ``-> delay_m`` (default: :math:`t=\tau_\text{M}/2``): By default, each pulse starts when the previous ends. This can be changed with this variable (if given, must be one value per pulse).
    * Polarization (see box :ref:`Defining the polarization <polarization>` below).

  .. note::
    When a time-dependent electric (magnetic) field is used, the string ``_efield`` (``_magfield``) is attached to the results filename. Besides that, if the fields are pulsed, the string ``pulse`` is also included.

.. _polarization: 

  .. admonition:: Defining the polarization
    :class: danger

    The polarization for the field is required. It can be given in two ways:

    * ``-> polarization_e`` (for electric field) and ``-> polarization_m`` (for magnetic field): This can be either one of:

      *  ``x``: :math:`\cos(\omega t)\mathbf{\hat{x}}`;
      *  ``y``: :math:`\cos(\omega t)\mathbf{\hat{y}}`;
      *  ``z``: :math:`\cos(\omega t)\mathbf{\hat{z}}`;
      *  ``p`` (:math:`+`): :math:`\cos(\omega t)\mathbf{\hat{x}}+\mathrm{i}\sin(\omega t)\mathbf{\hat{y}}`;
      *  ``m`` (:math:`-`): :math:`\cos(\omega t)\mathbf{\hat{x}}-\mathrm{i}\sin(\omega t)\mathbf{\hat{y}}`.

    * Since the options above are restricted, that input can be omitted and instead each polatization can be given by two vectors of 3 coordinates each. For the electric field, they are: ``-> polarization_vec_ip_e`` (in-phase, i.e., proportional to :math:`\cos(\omega t)`) and ``-> polarization_vec_op_e`` (out-of-phase, i.e., proportional to :math:`\sin(\omega t)`). The input is equivalent for magnetic fields, but with the parameters ``-> polarization_vec_ip_m`` and ``-> polarization_vec_op_m``. For example, the ``p`` polarization above is equivalent to an in-phase component ``1.0 0.0 0.0`` and an out-of-phase ``0.0 1.0 0.0``.

    **Note that for pulses, it is expected one polarization (either 1 letter or 3 numbers) for each pulse.**

For ``-> itype=12``:

These are the quantities related to the calculation of the magnetic coupling tensor in real space. 
For more details, see :doc:`magnetic tensor methods<magnetic_tensor>`.

* ``-> cluster_layers`` (default: :math:`2`): TODO: CHECK THE CODE, THERE'S NO DEFAULT ON MOD_IO
* ``-> nqpt`` (default: same as ``-> nkpt``): Number of wave vectors to perform the Fourier transform from :math:`J(\mathbf{q})` to :math:`J(\mathbf{R})`.


.. note::
  All the parameters read from :guilabel:`input` are logged into the :guilabel:`parameter.in` :doc:`output file<output_files>`.






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

* ``cartesian``, where the positions :math:`x, y, z` in each line are only multiplied by the lattice parameter. This means that the real position of the atoms are:

    :math:`\mathbf{r} = (x \mathbf{\hat{x}}+ y \mathbf{\hat{y}}+ z \mathbf{\hat{z}})a_0`

* ``bravais`` or ``lattice``, where the numbers :math:`n_1, n_2, n_3` are multiplied also by the Bravais vectors, such that

    :math:`\mathbf{r} = (n_1 \mathbf{a}_1+ n_2 \mathbf{a}_2+ n_3 \mathbf{a}_3)a_0`

The positions of all the basis atoms must be then listed starting on line 9.

.. important::
    The number of position lines should be the same number as the sum of the numbers in line 7.


.. index:: Parameter file
.. _parameter_file:

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
* SOC strength for p and d orbitals: :math:`\lambda_{\text{SOC}}` for `p` and `d` orbitals (the spherical `s` orbital does not contribute to the SOC term)
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
* SOC strength for p and d orbitals: :math:`\lambda_{\text{SOC}}` for `p` and `d` orbitals (the spherical `s` orbital does not contribute to the SOC term)
* Superconductivity :math:`\lambda`: This can be given in different ways, in the following order or priority:

  * A **single** value for all the orbitals;
  * **One value per orbital** of the given atom type;
  * **Three values**, one per general orbital type (`s`, `p` and `d`);
  * **Nine values**, one per specific orbital type (`s`, `px`, `py`, `pz`, `dxy`, `dyz`, `dzx`, `dx2`, `dz2`).



Secondary files
===============

Other files that may be necessary depending on the calculation:


.. index:: hamiltonian_file
.. _hamiltonian_file:

Hamiltonian file
----------------

The Hamiltonian obtained from DFT when ``-> tbmode = 2``


.. index:: kbands
.. _kbands:

:guilabel:`kbands`
------------------

The file :guilabel:`kbands` includes the definition of k-points to be used in the band structure (``itype=4``) and q-dependent susceptibility calculation (``itype=7``).
The units is defined on the :guilabel:`input` file via the parameter ``-> qbasis``.



.. index:: initialmag
.. _initialmag:

:guilabel:`initialmag`
----------------------

TITAN initial magnetic moments for the atoms in the unit cell is `2.0` along the `z-`direction.
This can be changed by adding ``-> magbasis`` in the :guilabel:`input` file (possible values are ``cartesian`` or ``spherical``), and adding the values in a file called :guilabel:`initialmag`.
It must contain 3 values per line, with the number of lines given by the number of atoms in the unit cell (following the same order as the positions).

* ``cartesian``: `\mathbf{m}_i^x`, `\mathbf{m}_i^y` and `\mathbf{m}_i^z` coordinates;
* ``spherical``: :math:`|\mathbf{m}_i|`, :math:`\mathbf{m}_i^\theta` and :math:`\mathbf{m}_i^\phi` coordinates;



