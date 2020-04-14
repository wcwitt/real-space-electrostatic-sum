# MIT License
# 
# Copyright (c) 2019-2020 William C. Witt
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import os
import ctypes as ct
import numpy as np

#______________________________________________________________________________
#                                                             ctypes setup

# load library
lib = ct.cdll.LoadLibrary(os.path.dirname(os.path.abspath(__file__)) 
                            + '/../library/real_space_electrostatic_sum.so')

# set argtypes and restype for 'c_real_space_electrostatic_sum_energy'
lib.c_real_space_electrostatic_sum_energy.argtypes = [
        ct.POINTER(ct.c_double), # a1
        ct.POINTER(ct.c_double), # a2
        ct.POINTER(ct.c_double), # a3
        ct.POINTER(ct.c_int),    # n
        ct.POINTER(ct.c_double), # rx
        ct.POINTER(ct.c_double), # ry
        ct.POINTER(ct.c_double), # rz
        ct.POINTER(ct.c_double), # z
        ct.POINTER(ct.c_double), # rc
        ct.POINTER(ct.c_double), # rd
        ct.POINTER(ct.c_double)] # e
lib.c_real_space_electrostatic_sum_energy.restype = None

# set argtypes and restype for 'c_real_space_electrostatic_sum_force'
lib.c_real_space_electrostatic_sum_force.argtypes = [
        ct.POINTER(ct.c_double), # a1
        ct.POINTER(ct.c_double), # a2
        ct.POINTER(ct.c_double), # a3
        ct.POINTER(ct.c_int),    # n
        ct.POINTER(ct.c_double), # rx
        ct.POINTER(ct.c_double), # ry
        ct.POINTER(ct.c_double), # rz
        ct.POINTER(ct.c_double), # z
        ct.POINTER(ct.c_double), # rc
        ct.POINTER(ct.c_double), # rd
        ct.POINTER(ct.c_double), # fx
        ct.POINTER(ct.c_double), # fy
        ct.POINTER(ct.c_double)] # fz
lib.c_real_space_electrostatic_sum_force.restype = None

# set argtypes and restype for 'c_real_space_electrostatic_sum_stress'
lib.c_real_space_electrostatic_sum_stress.argtypes = [
        ct.POINTER(ct.c_double), # a1
        ct.POINTER(ct.c_double), # a2
        ct.POINTER(ct.c_double), # a3
        ct.POINTER(ct.c_int),    # n
        ct.POINTER(ct.c_double), # rx
        ct.POINTER(ct.c_double), # ry
        ct.POINTER(ct.c_double), # rz
        ct.POINTER(ct.c_double), # z
        ct.POINTER(ct.c_double), # rc
        ct.POINTER(ct.c_double), # rd
        ct.POINTER(ct.c_double)] # s
lib.c_real_space_electrostatic_sum_stress.restype = None

#______________________________________________________________________________
#                                                                   energy

def energy(a1, a2, a3, n, rx, ry, rz, z, rc, rd):

    # create c variables (except for numpy arrays)
    n_c = ct.c_int(n)
    rc_c = ct.c_double(rc)
    rd_c = ct.c_double(rd)
    e_c = ct.c_double()

    # ensure numpy arrays are stored as expected
    a1_c = np.require(a1, dtype='float64', requirements=['C','A'])
    a2_c = np.require(a2, dtype='float64', requirements=['C','A'])
    a3_c = np.require(a3, dtype='float64', requirements=['C','A'])
    rx_c = np.require(rx, dtype='float64', requirements=['C','A'])
    ry_c = np.require(ry, dtype='float64', requirements=['C','A'])
    rz_c = np.require(rz, dtype='float64', requirements=['C','A'])
    z_c = np.require(z, dtype='float64', requirements=['C','A'])

    # call library function
    lib.c_real_space_electrostatic_sum_energy(
            a1_c.ctypes.data_as(ct.POINTER(ct.c_double)),
            a2_c.ctypes.data_as(ct.POINTER(ct.c_double)),
            a3_c.ctypes.data_as(ct.POINTER(ct.c_double)),
            ct.byref(n_c),
            rx_c.ctypes.data_as(ct.POINTER(ct.c_double)),
            ry_c.ctypes.data_as(ct.POINTER(ct.c_double)),
            rz_c.ctypes.data_as(ct.POINTER(ct.c_double)),
            z_c.ctypes.data_as(ct.POINTER(ct.c_double)),
            ct.byref(rc_c),
            ct.byref(rd_c),
            ct.byref(e_c))

    # return the energy
    return e_c.value

#______________________________________________________________________________
#                                                                    force

def force(a1, a2, a3, n, rx, ry, rz, z, rc, rd):

    # create c variables (except for numpy arrays)
    n_c = ct.c_int(n)
    rc_c = ct.c_double(rc)
    rd_c = ct.c_double(rd)

    # ensure numpy arrays are stored as expected
    a1_c = np.require(a1, dtype='float64', requirements=['C','A'])
    a2_c = np.require(a2, dtype='float64', requirements=['C','A'])
    a3_c = np.require(a3, dtype='float64', requirements=['C','A'])
    rx_c = np.require(rx, dtype='float64', requirements=['C','A'])
    ry_c = np.require(ry, dtype='float64', requirements=['C','A'])
    rz_c = np.require(rz, dtype='float64', requirements=['C','A'])
    z_c = np.require(z, dtype='float64', requirements=['C','A'])

    # create numpy arrays for forces
    fx = np.require(np.zeros(n, dtype='float64'), requirements=['C','A'])
    fy = np.require(np.zeros(n, dtype='float64'), requirements=['C','A'])
    fz = np.require(np.zeros(n, dtype='float64'), requirements=['C','A'])

    # call library function
    lib.c_real_space_electrostatic_sum_force(
            a1_c.ctypes.data_as(ct.POINTER(ct.c_double)),
            a2_c.ctypes.data_as(ct.POINTER(ct.c_double)),
            a3_c.ctypes.data_as(ct.POINTER(ct.c_double)),
            ct.byref(n_c),
            rx_c.ctypes.data_as(ct.POINTER(ct.c_double)),
            ry_c.ctypes.data_as(ct.POINTER(ct.c_double)),
            rz_c.ctypes.data_as(ct.POINTER(ct.c_double)),
            z_c.ctypes.data_as(ct.POINTER(ct.c_double)),
            ct.byref(rc_c),
            ct.byref(rd_c),
            fx.ctypes.data_as(ct.POINTER(ct.c_double)),
            fy.ctypes.data_as(ct.POINTER(ct.c_double)),
            fz.ctypes.data_as(ct.POINTER(ct.c_double)))

    # return the forces
    return fx, fy, fz

#______________________________________________________________________________
#                                                                   stress

def stress(a1, a2, a3, n, rx, ry, rz, z, rc, rd):

    # create c variables (except for numpy arrays)
    n_c = ct.c_int(n)
    rc_c = ct.c_double(rc)
    rd_c = ct.c_double(rd)

    # ensure numpy arrays are stored as expected
    a1_c = np.require(a1, dtype='float64', requirements=['C','A'])
    a2_c = np.require(a2, dtype='float64', requirements=['C','A'])
    a3_c = np.require(a3, dtype='float64', requirements=['C','A'])
    rx_c = np.require(rx, dtype='float64', requirements=['C','A'])
    ry_c = np.require(ry, dtype='float64', requirements=['C','A'])
    rz_c = np.require(rz, dtype='float64', requirements=['C','A'])
    z_c = np.require(z, dtype='float64', requirements=['C','A'])

    # create numpy array for stress
    stress = np.require(np.zeros(6, dtype='float64'), requirements=['C','A'])

    # call library function
    lib.c_real_space_electrostatic_sum_stress(
            a1_c.ctypes.data_as(ct.POINTER(ct.c_double)),
            a2_c.ctypes.data_as(ct.POINTER(ct.c_double)),
            a3_c.ctypes.data_as(ct.POINTER(ct.c_double)),
            ct.byref(n_c),
            rx_c.ctypes.data_as(ct.POINTER(ct.c_double)),
            ry_c.ctypes.data_as(ct.POINTER(ct.c_double)),
            rz_c.ctypes.data_as(ct.POINTER(ct.c_double)),
            z_c.ctypes.data_as(ct.POINTER(ct.c_double)),
            ct.byref(rc_c),
            ct.byref(rd_c),
            stress.ctypes.data_as(ct.POINTER(ct.c_double)))

    # return the stress
    return stress
