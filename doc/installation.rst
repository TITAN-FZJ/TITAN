.. index:: Installation
.. sectionauthor:: Filipe Guimarães <f.guimaraes@fz-juelich.de>
    
************
Installation
************

TITAN is written in `Modern Fortran <https://fortran-lang.org/>`_, and parallelized using MPI and openMP.
An initial version using accelerators is also implemented using CUDA Fortran and the `cuSolver library to solve the eigenvalue problem <https://docs.nvidia.com/cuda/cusolver/index.html#cuSolverDN-lt-t-gt-syevd>`_.
TITAN's compilation uses CMake.

It was mainly configured and optimized to be used on the supercomputers of the Jülich Supercomputing Centre, so the compilations steps will be exemplified there.
To compile TITAN, first :doc:`download<download>` the source code to your local computer or supercomputer of choice.

To use the intel compiler:

.. code-block:: bash

    module load Intel IntelMPI CMake

or to use GFortran:

.. code-block:: bash

    module load GCC ParaStationMPI CMake

(Could also be done using ``OpenMPI`` instead of ``ParaStationMPI``.)

The compilation should be done via the ``${TITAN}/build`` folder.
To avoid incompatibility with previous compilations, it is recommended to clean this folder before starting a new (different) compilation.
So, to compile the code, the following should be done:

.. code-block:: bash

    cd build
    rm -r *
    cmake ../ -DPLATFORM=<system> [-DCOMPILER=<compiler> -DDEBUG=ON -DSUFFIX=<suffix> -DPREP=ON]
    make -j 32

The following options are possible:

* ``-DPLATFORM=<system>``: The ``<system>`` is defined for the Intel compiler, since the flags were different depending on the system. The following options are possible:

    * ``jurecadc``, ``juwelsbooster`` or ``jusuf``
    * ``juwels``
    * ``jurecabooster``
    * ``jusuf``
    * ``rwth`` for the Claix supercomputers
    * ``iff`` for the PGI-1 local clusters
    * ``osx`` for MacOS

    Not all of them have different flags, but they can be adapted separately in the ``CMakeLists.txt`` file.
* ``-DCOMPILER=<compiler>``: CMake looks for `mpiifort` and `mpif90` (in this order) as compiler. The `COMPILER` flag can be used in case something is wrongly recognized. Possible values are: ``intel``, ``gfortran`` (or ``gnu``), or ``nvfortran``. For the last two, the executable filename is appended with ``_gnu`` or ``_nv``, respectively.
* ``-DDEBUG=ON``: The debug flags are not used by default. Use this option to activate them. When this is used, the executable is appended with ``_debug``.
* ``-DSUFFIX=<suffix>``: This variable can be used to add an extra suffix to the executable.
* ``-DPREP=ON``: Use this option to compile with `SCORE-P` instrumentation. This needs the module ``Score-P``, and the results can be analised with `Scalasca <https://www.scalasca.org/>`_ (via the module `Scalasca`).

The generated binary can be found in ``${TITAN}/bin/``.