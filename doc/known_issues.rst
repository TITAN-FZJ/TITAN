.. index:: known_issues

************
Known issues
************

Here we list the known issues on TITAN.

Currents responses
==================

In the past, TITAN had hard-coded unit cells, and could calculate response functions
of the type :math:`\chi_{ijll}` or :math:`\chi_{iikl}`. 
With the change to the generalized lattice, the responses should also be generalized to :math:`\chi_{ijkl}` (in principle),
which represents current from site `i` to site `j` due to a current flowing between sites `k` and `l`.
However, the relation between RPA and HF, given in a simple form by :math:`\chi = \chi_{HF} - \chi_{HF}.U.\chi` can be solved for the first
two cases, but become complicated for the most general 4-indices one.
Ideally, the 

There are some possible work-arounds:

* Implement the calculation for the same 3-indices, but that may be difficult to identify the cases and build the matrix;
* Implement the HF resposes only, which could give already valuable information on some cases where :math:`U=0`.

