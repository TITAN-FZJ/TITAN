# Input File Documentation

All commands are written with a preceeding `->` and followed by an equal sign `=` which separates command and value, e.g.
```
  -> itype = 1
```
## Mandatory commands
### Type of Calculation: `itype`
1 - Only Self-Consistency

2 - (Debugging)

3 - Local Density of States and Coupling

4 - Band Structure

5 - Fermi Surface 

6 - Coupling

7 - RPA & HF Susceptibilities and Gilbert Damping from Slope

8 - Calculate everything (except 2-6)

9 - Calculate everything in DC limit (except 2-6)

10 - Calculate Gilbert from Torque Correlation Model with SO and XC Torques

### Name of output log file: `output`
```fortran
character(len=200) :: output
```
### Number of next neighbor stages: `nn_stages`
```fortran
integer :: nn_stages
```
### Broadening of Greens function: `eta`
```fortran
real(double) :: eta
```

### Status of static magnetic field: `field`
```fortran
logical :: field
```
In case the static magnetic field is turned on, more parameter become mandatory:
#### Static magnetic field amplitude: `hwa`
```fortran
real(double), dimension(:) :: hwa
```
If `dimension(1)` use only a single amplitude