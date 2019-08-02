<h1> TITAN </h1>

**TI**me-dependent 
**T**ransport and
**A**ngular momentum of
**N**anostructures

A paper containing a detailed description of the method used in the code is in preparation. In the meantime, we shall greatly appreciate if scientific work done using **TITAN** will contain an explicit acknowledgement and the following reference:<br>
F. S. M. Guimarães et al., Sci. Rep. 7, 3686 (2017)<br>
https://www.nature.com/articles/s41598-017-03924-1

This program calculates the electric and spin excitations for bulk
and thin films. It can describe Ferromagnetic Resonance (FMR) experiments
and intrinsic Spin and Orbital Momentum Hall Effects (including Anomalous 
and Planar Hall effects), for example.
A multi-orbital tight-binding model is used, with the el-el interaction 
described by a Hubbard-like hamiltonian [2].  The spin-orbit interaction
is present, and a static magnetic field can be applied in any direction.
Excitations are obtained in Kubo's linear response approach [1].

An applied AC electric field (described by a time-dependent potential
vector) that couples to the current density is applied, and the following
responses are calculated (for each site in the unit cell):
 - local spin accumulations;
 - local orbital momentum accumulations in each site;
 - spin-orbit, exchange and external torques;
 - effective magnetic fields;
 - spin and charge currents*;
 - orbital momentum currents*;
 - DC components of pumped spin currents*.<br>
(* currents need reimplementation for the new generalized cell)

The full 4x4 magnetic susceptibility matrix, including transverse and 
longitudinal blocks, is also evaluated.

The LDOS in each site, band structure, Fermi surface, and full exchange 
coupling tensor (including DMI and anisotropic terms are also implemented.

The program is currently written for 9 orbitals per site 
(1 <i>s</i>,3 <i>p</i> and 5 <i>d</i>), for a given set of 
nearest neighbours hopping parameters (as described by Slater and Koster [3]).
The generalized response functions are obtained as a function of the frequency
within the Random-Phase Approximation, which are calculated from the mean-field counterpart.
These ones are evaluated in terms of the monoeletronic Green functions. 
The integration in k-space (2D or 3D) is calculated using a uniformly spaced k-mesh.
The integration over energy is done in the complex plane - part in the imaginary axis
and part in the real axis (for &#969;&#8800;0).
The code is parallelized using openMP and MPI.

Plotting scripts (in Python) are provided in the 'scripts' folder.

[1] R. Kubo, Statistical-Mechanical Theory of Irreversible Processes. I. General Theory and Simple Applications to Magnetic and Conduction Problems, J. Phys. Soc. Jpn. 12, 570 (1957).

[2] J. Hubbard, Electron Correlations in Narrow Energy Bands, Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences 276, 238 (1963).

[3] J. C. Slater and G. F. Koster, Simplified LCAO Method for the Periodic Potential Problem, Phys. Rev. 94, 1498 (1954).

OBS: Legacy version uses DFT parameters as described in

[4] O. K. Andersen and O. Jepsen, Phys. Rev. Lett. 53, 2571 (1984).

```
!*******************************************************************************!
!                           CHOOSE WHAT TO CALCULATE                            !
!  Options:                                                                     !
!   0 - Test before SC                                                          !
!   1 - Self consistency only                                                   !
!   2 - Test after SC                                                           !
!   3 - LDOS and exchange interactions as a function of energy                  !
!   4 - Band structure                                                          !
!   5 - Charge and spin density at Fermi surface                                !
!   6 - Exhange interactions and anisotropies (full tensor Jij)                 !
!   7 - Local susceptibility as a function of energy                            !
!   8 - Parallel currents, disturbances, susceptibilities, torques and          !
!       effective fields as a function of energy                                !
!   9 - dc-limit as a function of magnetic field (abs, theta or phi)            !
!       (use 'emin' as frequency)                                               !
!  10 - Calculation of Gilbert Damping using Torque Correlation models          !
!  11 - Real-time propagation (NEW!)                                            !
!===============================================================================!
!                         OPTIONAL RUNNING VARIABLES                            !
!   noUonall      - Set U=0 on all layers                                       !
!   noUonNM       - Set U=0 on NM layers                                        !
!   GSL           - Calculate ground state of Orbital Angular Momentum L        !
!   verbose       - Write intermediate steps and progress ars                   !
!   ry2ev         - Convert energy units to eV                                  !
!   tesla         - Read values of static magnetic field in tesla               !
!   debug         - Turn on to use this optional variable for debugging         !
!   createfiles   - Just create files (WARNING: overwrite existing files!)      !
!   addresults    - Add results to existing files, instead of creating new ones !
!                   (files MUST exist)                                          !
!   slatec        - Use non-linear root finder dnsqe from Slatec library        !
!   kpoints       - Write kpoints and weights into files                        !
!   lineargfsoc   - Adds linear order SOC in the GF                             !
!   linearsoc     - Adds linear order SOC in the susceptibility                 !
!   nojac         - Perform self-consistency without using the jacobian         !
!   hfresponses   - Use Hartree-Fock response functions                         !
!   ontheflysc    - Writes intermediate self-consistency steps to file          !
!   rotatemag     - Rotate the magnetization to the direction of the field      !
!   nolb          - Remove Orbital Zeeman from the hamiltonian                  !
!   nodiag        - Don't diagonalize susceptibilities                          !
!   sortfiles     - Only sort files (no calculation is done)                    !
!   writeonscreen - Also write some results on output files                     !
!   sha           - Only calculate the SHA for input parameters                 !
!   lgtv          - Only calculate total longitudinal and transverse currents   !
!===============================================================================!
!                           INTEGRATION VARIABLES                               !
! k-points are obtained using Cunningham special points for BCC110 and FCC100   !
! Ref.: S. Cunningham, Phys. Rev. B 10, 4988 (1974)                             !
!                                                                               !
! Energy integration use Gauss-Legendre quadrature                              !
!===============================================================================!
!                                BAND STRUCTURE                                 !
! K-points must be listed in file 'kbands' and the path is given in the input.  !
! e.g.:                                                                         !
-> band = G N M G
! where the file 'kbands' contains                                              !
G	0.0	0.0	0.0
N	0.0	0.5	0.0
M 0.5 0.5 0.0
! This can also be used to calculate the q-dependent susceptibility             !
!*******************************************************************************!
! To stop the program:                                                          !
! - after an energy step, just create a file called 'stop'                      !
! (e.g., in terminal: touch stop)                                               !
! - after any layer/field calculation, just create a file called 'stopout'      !
! (e.g., in terminal: touch stopout)                                            !
!*******************************************************************************!
```

## Compiling the code

First clone the repository on your local computer or supercomputer of choice.
```
git clone https://iffgit.fz-juelich.de/titan/TITAN.git
cd TITAN/build
``` 
Depending on the system you are compiling on choose one of the following options.
### Jureca 
```
module load Intel IntelMPI NAG/Mark26 CMake
cmake ../ -DPLATFORM=jureca 
make -j12
```

### Juqueen 
```
module load nag 
cmake ../ -DPLATFORM=juqueen
make -j16
```

### Linux 
```
cmake ../ 
make -j4
```

The binary can be found in `TITAN/bin/`
## Running the code 

After compiling you can start running the code.
An example calculation can be found in `TITAN/example`

### Jureca 

```
#!/bin/bash -x 
#SBATCH --nodes=1
#SBATCH -J JobName
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH -e error-%j
#SBATCH -o output-%j
#SBATCH --mail-user=your@mail
#SBATCH --mail-type=ALL
#SBATCH --partition=batch 
#SBATCH --time=1:00:00

module load Intel IntelMPI NAG/Mark26

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=1

export OMP_STACKSIZE=150m
export KMP_AFFINITY=scatter

srun ~/TITAN/bin/titan_jureca.exe --exports=NAG_KUSARI_FILE --exports=OMP_NUM_THREADS --exports=MKL_NUM_THREADS --exports=KMP_AFFINITY
``` 

### Juqueen
```
# @ job_name = JobName
# @ error = error-$(job_name).$(jobid).err
# @ output = output-$(job_name).$(jobid).out
# @ environment = COPY_ALL
# @ wall_clock_limit = 24:00:00
# @ notification = always
# @ notify_user = your@mail
# @ job_type = bluegene
## @ bg_connectivity = TORUS
# @ bg_size = 1024
# @ queue

# IMPORTANT RELATIONS
#
# np = number of points*256 (2*64 + 2*64)
NUMPROCS=1024
# bg_size*ranks_per_node = np (or bg_size = np/ranks_per_node)
# threads_per_cpu*ranks_per_node = 64
NUMTHREADS=64
RANKSPERNODE=1
#

module load nag

runjob --np ${NUMPROCS} --exp-env NAG_KUSARI_FILE --envs OMP_NUM_THREADS=${NUMTHREADS} --ranks-per-node ${RANKSPERNODE} --args " " --exe ~/TITAN/bin/titan_juqueen.exe
```

### Linux 
``` 
export OMP_NUM_THREADS=8
mpirun -np 1 ~/TITAN/bin/titan.exe --exports=OMP_NUM_THREADS
```
