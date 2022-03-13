.. index:: Installation
.. sectionauthor:: Filipe Guimarães <f.guimaraes@fz-juelich.de>
    
************
Installation
************

TITAN is written in FORTRAN, and parallelized using MPI and openMP.
An initial version using accelerators is also implemented using CUDA Fortran.
Its compilation uses CMake.


It was mainly configured and optimized to be used on the supercomputers of the Jülich Supercomputing Centre, so the compilations steps will be exemplified.
To compile TITAN, first :doc:`download<download>` the source code to your local computer or supercomputer of choice.

To use the intel compiler:

.. code-block:: bash

    module load Intel IntelMPI CMake

or to use GFortran:

.. code-block:: bash

    module load GCC ParaStationMPI CMake

Then, compile the code using

.. code-block:: bash

    cmake ../ -DPLATFORM=system [-DCOMPILER=compiler -DDEBUG=ON]
    make -j 24


The ``system`` is defined for the Intel compiler, and can be ``jurecadc``, ``juwels``, ``jusuf``, or ``iff``.
TITAN was tested for the compilers: ``ifort`` (default), ``gfortran`` (gnu) and ``nvfortran``.
The generated binary can be found in ``${TITAN}/bin/``, and its filename contains the system, compiler (if not ifort) and ``debug`` if the respective flag was used.