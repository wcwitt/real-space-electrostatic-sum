[![Build Status](https://travis-ci.com/wcwitt/real-space-electrostatic-sum.svg?branch=master)](https://travis-ci.com/wcwitt/real-space-electrostatic-sum)
# real-space-electrostatic-sum

Implementation of the real-space electrostatic sum outlined in [Pickard, *Phys. Rev. Mat.* **2**, 013806, 2018](https://doi.org/10.1103/PhysRevMaterials.2.013806).

Potentially faster than the ubiquitous Ewald sum found in many electronic structure codes and elsewhere.

Repository contains:

* a [Fortran module](src/real_space_electrostatic_sum.f90) with the main routines;
* a [C-style interface](src/real_space_electrostatic_sum_c.f90);
* a [Python wrapper](python/real_space_electrostatic_sum.py) built with ctypes;
* a [Jupyter notebook](https://nbviewer.jupyter.org/github/wcwitt/real-space-electrostatic-sum/blob/master/python/benchmarking.ipynb) with examples and benchmarking.

Now with forces and stresses!

Warning: The "loc" array expects Fortran ordering, so beware when calling from C/C++.
