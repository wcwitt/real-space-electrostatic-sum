cmake_minimum_required(VERSION 3.1)
project(real_space_electrostatic_sum LANGUAGES Fortran C)

add_library(real_space_electrostatic_sum 
    SHARED
    source/real_space_electrostatic_sum.f90
    source/c_real_space_electrostatic_sum.f90)

target_include_directories(real_space_electrostatic_sum PUBLIC include)
