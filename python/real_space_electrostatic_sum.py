# MIT License
# 
# Copyright (c) 2019 William C. Witt
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

#!/usr/bin/env/python

import ctypes as ct
import numpy as np

#______________________________________________________________________________
#                                                             ctypes setup

# load library
lib = ct.cdll.LoadLibrary('../lib/real_space_electrostatic_sum.so')

# set argtypes and restype for 'energy_c'
lib.energy_c.argtypes = [ct.POINTER(ct.c_double), # a_x
                         ct.POINTER(ct.c_double), # a_y
                         ct.POINTER(ct.c_double), # a_z
                         ct.POINTER(ct.c_int),    # num
                         ct.POINTER(ct.c_double), # loc
                         ct.POINTER(ct.c_double), # chg
                         ct.POINTER(ct.c_double), # r_c
                         ct.POINTER(ct.c_double), # r_d
                         ct.POINTER(ct.c_double)] # ene
lib.energy_c.restype = None

#______________________________________________________________________________
#                                                                   energy

def energy(a_x, a_y, a_z, num, loc, chg, r_c, r_d):

    # create c variables (except for numpy arrays)
    num_c = ct.c_int(num)
    r_c_c = ct.c_double(r_c)
    r_d_c = ct.c_double(r_d)
    ene_c = ct.c_double()

    # ensure numpy arrays are stored as expected
    a_x_c = np.require(a_x, dtype='float64', requirements=['F_CONTIGUOUS', 'ALIGNED'])
    a_y_c = np.require(a_y, dtype='float64', requirements=['F_CONTIGUOUS', 'ALIGNED'])
    a_z_c = np.require(a_z, dtype='float64', requirements=['F_CONTIGUOUS', 'ALIGNED'])
    loc_c = np.require(loc, dtype='float64', requirements=['F_CONTIGUOUS', 'ALIGNED'])
    chg_c = np.require(chg, dtype='float64', requirements=['F_CONTIGUOUS', 'ALIGNED'])

    # call library function
    lib.energy_c(a_x_c.ctypes.data_as(ct.POINTER(ct.c_double)),
                 a_y_c.ctypes.data_as(ct.POINTER(ct.c_double)),
                 a_z_c.ctypes.data_as(ct.POINTER(ct.c_double)),
                 ct.byref(num_c),
                 loc_c.ctypes.data_as(ct.POINTER(ct.c_double)),
                 chg_c.ctypes.data_as(ct.POINTER(ct.c_double)),
                 ct.byref(r_c_c),
                 ct.byref(r_d_c),
                 ct.byref(ene_c))

    # extract the result
    return ene_c.value
