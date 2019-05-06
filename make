#!/bin/bash

if [ -d obj ]; then rm -r obj; fi
mkdir obj
if [ -d mod ]; then rm -r mod; fi
mkdir mod
if [ -d lib ]; then rm -r lib; fi
mkdir lib

# flags
flags="-g -Og -fimplicit-none -fcheck=all -fbacktrace -Wall -Wextra -Wconversion -pedantic"

# compile into a shared library
gfortran -shared -fPIC $flags -Jmod \
    src/real_space_electrostatic_sum.f90 src/real_space_electrostatic_sum_c.f90 \
    -o lib/real_space_electrostatic_sum.so

