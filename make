#!/bin/bash

if [ -d module ]; then rm -r module; fi
mkdir module
if [ -d library ]; then rm -r library; fi
mkdir library

# flags
flags="-g -Og -fimplicit-none -fcheck=all -fbacktrace -Wall -Wextra -Wconversion -pedantic"

# compile into a shared library
gfortran -shared -fPIC $flags -Jmodule \
    source/real_space_electrostatic_sum.f90 source/c_real_space_electrostatic_sum.f90 \
    -o library/real_space_electrostatic_sum.so

