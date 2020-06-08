[![Build Status](https://travis-ci.com/wcwitt/real-space-electrostatic-sum.svg?branch=master)](https://travis-ci.com/wcwitt/real-space-electrostatic-sum)
# real-space-electrostatic-sum

Implementation of the real-space electrostatic sum outlined in [Pickard, *Phys. Rev. Mat.* **2**, 013806, 2018](https://doi.org/10.1103/PhysRevMaterials.2.013806). Includes force and stress routines.

Potentially faster than the ubiquitous Ewald sum found in many electronic structure codes and elsewhere.

Repository contains:

* a [Fortran module](source/real_space_electrostatic_sum.f90) with the main routines;
* a [C-style interface](source/c_real_space_electrostatic_sum.f90);
* a [Python wrapper](python/real_space_electrostatic_sum.py) built with ctypes;
* a [Jupyter notebook](https://nbviewer.jupyter.org/github/wcwitt/real-space-electrostatic-sum/blob/master/python/benchmarking.ipynb) with examples and benchmarking;
* a [set of unit tests](https://github.com/wcwitt/real-space-electrostatic-sum/blob/master/test/test.py) for the energy, force, and stress routines.

To build and test:

```
mkdir build && cd build
cmake ..
make
cd ../test
python test.py
```
