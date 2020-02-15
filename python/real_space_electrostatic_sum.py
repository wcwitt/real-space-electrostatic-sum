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

# set argtypes and restype for 'c_energy'
lib.c_energy.argtypes = [ct.POINTER(ct.c_double), # a1
                         ct.POINTER(ct.c_double), # a2
                         ct.POINTER(ct.c_double), # a3
                         ct.POINTER(ct.c_int),    # num
                         ct.POINTER(ct.c_double), # loc
                         ct.POINTER(ct.c_double), # z
                         ct.POINTER(ct.c_double), # rc
                         ct.POINTER(ct.c_double), # rd
                         ct.POINTER(ct.c_double)] # ene
lib.c_energy.restype = None

# set argtypes and restype for 'c_force'
lib.c_force.argtypes = [ct.POINTER(ct.c_double), # a1
                        ct.POINTER(ct.c_double), # a2
                        ct.POINTER(ct.c_double), # a3
                        ct.POINTER(ct.c_int),    # num
                        ct.POINTER(ct.c_double), # loc
                        ct.POINTER(ct.c_double), # z
                        ct.POINTER(ct.c_double), # rc
                        ct.POINTER(ct.c_double), # rd
                        ct.POINTER(ct.c_double), # fx
                        ct.POINTER(ct.c_double), # fy
                        ct.POINTER(ct.c_double)] # fz
lib.c_force.restype = None

# set argtypes and restype for 'c_stress'
lib.c_stress.argtypes = [ct.POINTER(ct.c_double), # a1
                         ct.POINTER(ct.c_double), # a2
                         ct.POINTER(ct.c_double), # a3
                         ct.POINTER(ct.c_int),    # num
                         ct.POINTER(ct.c_double), # loc
                         ct.POINTER(ct.c_double), # z
                         ct.POINTER(ct.c_double), # rc
                         ct.POINTER(ct.c_double), # rd
                         ct.POINTER(ct.c_double)] # s
lib.c_stress.restype = None

#______________________________________________________________________________
#                                                                   energy

def energy(a1, a2, a3, num, loc, z, rc, rd):

    # create c variables (except for numpy arrays)
    num_c = ct.c_int(num)
    rc_c = ct.c_double(rc)
    rd_c = ct.c_double(rd)
    ene_c = ct.c_double()

    # ensure numpy arrays are stored as expected
    a1_c = np.require(a1, dtype='float64', requirements=['F_CONTIGUOUS', 'ALIGNED'])
    a2_c = np.require(a2, dtype='float64', requirements=['F_CONTIGUOUS', 'ALIGNED'])
    a3_c = np.require(a3, dtype='float64', requirements=['F_CONTIGUOUS', 'ALIGNED'])
    loc_c = np.require(loc, dtype='float64', requirements=['F_CONTIGUOUS', 'ALIGNED'])
    z_c = np.require(z, dtype='float64', requirements=['F_CONTIGUOUS', 'ALIGNED'])

    # call library function
    lib.c_energy(a1_c.ctypes.data_as(ct.POINTER(ct.c_double)),
                 a2_c.ctypes.data_as(ct.POINTER(ct.c_double)),
                 a3_c.ctypes.data_as(ct.POINTER(ct.c_double)),
                 ct.byref(num_c),
                 loc_c.ctypes.data_as(ct.POINTER(ct.c_double)),
                 z_c.ctypes.data_as(ct.POINTER(ct.c_double)),
                 ct.byref(rc_c),
                 ct.byref(rd_c),
                 ct.byref(ene_c))

    # return the energy
    return ene_c.value

#______________________________________________________________________________
#                                                                    force

def force(a1, a2, a3, num, loc, z, rc, rd):

    # create c variables (except for numpy arrays)
    num_c = ct.c_int(num)
    rc_c = ct.c_double(rc)
    rd_c = ct.c_double(rd)

    # ensure numpy arrays are stored as expected
    a1_c = np.require(a1, dtype='float64', requirements=['F_CONTIGUOUS', 'ALIGNED'])
    a2_c = np.require(a2, dtype='float64', requirements=['F_CONTIGUOUS', 'ALIGNED'])
    a3_c = np.require(a3, dtype='float64', requirements=['F_CONTIGUOUS', 'ALIGNED'])
    loc_c = np.require(loc, dtype='float64', requirements=['F_CONTIGUOUS', 'ALIGNED'])
    z_c = np.require(z, dtype='float64', requirements=['F_CONTIGUOUS', 'ALIGNED'])

    # create numpy arrays for forces
    fx = np.require(np.zeros(num, dtype='float64'), requirements=['F_CONTIGUOUS', 'ALIGNED'])
    fy = np.require(np.zeros(num, dtype='float64'), requirements=['F_CONTIGUOUS', 'ALIGNED'])
    fz = np.require(np.zeros(num, dtype='float64'), requirements=['F_CONTIGUOUS', 'ALIGNED'])

    # call library function
    lib.c_force(a1_c.ctypes.data_as(ct.POINTER(ct.c_double)),
                a2_c.ctypes.data_as(ct.POINTER(ct.c_double)),
                a3_c.ctypes.data_as(ct.POINTER(ct.c_double)),
                ct.byref(num_c),
                loc_c.ctypes.data_as(ct.POINTER(ct.c_double)),
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

def stress(a1, a2, a3, num, loc, z, rc, rd):

    # create c variables (except for numpy arrays)
    num_c = ct.c_int(num)
    rc_c = ct.c_double(rc)
    rd_c = ct.c_double(rd)

    # ensure numpy arrays are stored as expected
    a1_c = np.require(a1, dtype='float64', requirements=['F_CONTIGUOUS', 'ALIGNED'])
    a2_c = np.require(a2, dtype='float64', requirements=['F_CONTIGUOUS', 'ALIGNED'])
    a3_c = np.require(a3, dtype='float64', requirements=['F_CONTIGUOUS', 'ALIGNED'])
    loc_c = np.require(loc, dtype='float64', requirements=['F_CONTIGUOUS', 'ALIGNED'])
    z_c = np.require(z, dtype='float64', requirements=['F_CONTIGUOUS', 'ALIGNED'])

    # create numpy array for stress
    stress = np.require(np.zeros(6, dtype='float64'), requirements=['F_CONTIGUOUS', 'ALIGNED'])

    # call library function
    lib.c_stress(a1_c.ctypes.data_as(ct.POINTER(ct.c_double)),
                 a2_c.ctypes.data_as(ct.POINTER(ct.c_double)),
                 a3_c.ctypes.data_as(ct.POINTER(ct.c_double)),
                 ct.byref(num_c),
                 loc_c.ctypes.data_as(ct.POINTER(ct.c_double)),
                 z_c.ctypes.data_as(ct.POINTER(ct.c_double)),
                 ct.byref(rc_c),
                 ct.byref(rd_c),
                 stress.ctypes.data_as(ct.POINTER(ct.c_double)))

    # return the stress
    return stress
