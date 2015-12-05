<h1> PROGRAM DHE (Dynamical Hall Effects) </h1>

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
 - local spin densities in each layer;
 - local orbital momentum densities in each layer;
 - spin currents that flow parallel to the layers;
 - orbital momentum currents that flow parallel to the layers;
 - charge currents that flow parallel to the layers;
 - effective magnetic fields.

We can also calculate LDOS in each layer, band structure, Fermi surface,
and full exchange coupling tensor (including DMI and anisotropic terms).

We calculate the generalized response functions as a function of the frequency
within Random Phase Approximation, writing them in terms
of the mean field counterparts. These ones are written in terms of the
monoeletronic Green functions. The integration in k<sub>//</sub> is calculated
using the generation of 2D points by Cunningham.
It is parallelized using openMP. The number of points in the energy
can be set in the input file (usually, it uses a set of 128 points
in the imaginary axis and 64 points in the real axis (for &#969;&#8800;0)
and it is parallelized with MPI.
