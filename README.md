<h1> TITAN </h1>

**T**ime-dependent description of
**I**tinerant electrons:
**T**ransport and
**A**ngular momentum properties of
**N**anostructures

This program calculates the intrinsic Spin and Orbital Momentum Hall
Effects (including the Anomalous and Planar Hall effects) for
ultrathin films composed by (heavy + magnetic) transition metals
in a system where the spin-orbit interaction is present. We use
a multi-orbital tight-binding model in Kubo's linear response approach.
The el-el interaction is described by a Hubbard-like hamiltonian.
The program is written for 9 orbitals per site (1 <i>s</i>,3 <i>p</i> and 5 <i>d</i>),
including first and second nearest neighbours hopping matrices.

We apply an AC electric field (described by a time-dependent potential
vector) that couples to the current density, and we calculate the
response of:
 - local spin accumulations in each layer;
 - local orbital momentum accumulations in each layer;
 - spin and charge currents that flow parallel to the layers;
 - orbital momentum currents that flow parallel to the layers;
 - spin-orbit, exchange and external torques;
 - effective magnetic fields;
 - DC components of pumped spin currents.

We can also calculate LDOS in each layer, band structure, Fermi surface,
and full exchange coupling tensor (including DMI and anisotropic terms).

We calculate the generalized response functions as a function of the frequency
within Random Phase Approximation, writing them in terms
of the mean field counterparts. These ones are written in terms of the
monoeletronic Green functions. The integration in k<sub>//</sub> is calculated
using the generation of 2D points by Cunningham.
It is parallelized using openMP. The number of points in the energy
can be set in the input file - usually, it uses a set of 128 points
in the imaginary axis and 64 points in the real axis (for &#969;&#8800;0)
and it is parallelized with MPI.
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
!                             LATTICE DEFINITIONS                               !
!         For spin quantization and direction of in-plane electric field        !
!---------------------------------- BCC (110) ----------------------------------!
! L - long axis (x - short axis; y - out-of-plane)                              !
! S - short axis (x - long axis; y - out-of-plane)                              !
! P - perpendicular out-of plane (x - short axis; y - long axis)                !
!                                        ^(only for quantization)               !
! O - other (specify in main program)  -> (only for external field)             !
! or choose the direction of one of the neighbors                               !
!                           L direction                                         !
!                               |                                               !
!                            1  |  4                                            !
!                             \ | /                                             !
!                   ------5-----|-----6-------  S direction                     !
!                             / | \                                             !
!                            2  |  3                                            !
!                               |                                               !
! Layers grow in the (-0.5,0.5,0.0) direction (in the lattice referential)      !
!---------------------------------- FCC (100) ----------------------------------!
! In the following, choose the direction of one of the neighbors                !
!                (1-4) - first n.n. ; (5-8) second n.n.                         !
!              (use 9 for out-of-plane z, x // dirEfield)                       !
!                               6                                               !
!                               |                                               !
!                            2  |  1                                            !
!                             \ | /                                             !
!                   ------7-----|-----5-------                                  !
!                             / | \                                             !
!                            3  |  4                                            !
!                               |                                               !
!                               8                                               !
! Layers grow in the +z direction (in the lattice referential)                  !
!===============================================================================!
!                           INTEGRATION VARIABLES                               !
! k-points are obtained using Cunningham special points for BCC110 and FCC100   !
! Ref.: S. Cunningham, Phys. Rev. B 10, 4988 (1974)                             !
!                                                                               !
! Energy integration use Gauss-Legendre quadrature                              !
!===============================================================================!
!                                BAND STRUCTURE                                 !
! Choose one of the following high symmetry direction in k space:               !
!---------------------------------- BCC (110) ----------------------------------!
! GH - HP - PN - NG  (or the opposite HG - PH - NP - GN)                        !
!---------------------------------- FCC (100) ----------------------------------!
! GM - MX - XG  (or the opposite GX - XM - MG)                                  !
!===============================================================================!
!                                DFT PARAMETERS                                 !
! Ref.:	O. K. Andersen and O. Jepsen, Phys. Rev. Lett. 53, 2571 (1984).         !
!*******************************************************************************!
! To stop the program:                                                          !
! - after an energy step, just create a file called 'stop'                      !
! (e.g., in terminal: touch stop)                                               !
! - after any layer/field calculation, just create a file called 'stopout'      !
! (e.g., in terminal: touch stopout)                                            !
!*******************************************************************************!
```
