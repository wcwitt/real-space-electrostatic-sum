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

import numpy as np
import os
import sys
import unittest

sys.path.append(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), '../python/'))
import real_space_electrostatic_sum

class TestRealSpaceElectrostaticSum(unittest.TestCase):

    def test_energy(self):

        # Al
        a_1 = np.array([5.41141973394663, 0.00000000000000, 0.00000000000000])
        a_2 = np.array([2.70570986697332, 4.68642696013821, 0.00000000000000])
        a_3 = np.array([2.70570986697332, 1.56214232004608, 4.41840571073226])
        loc = np.zeros([1,3])
        chg = 3.0 * np.ones(loc.shape[0])
        h_max = 4.42
        r_d_hat = 2.0
        ene = real_space_electrostatic_sum.energy(
                a_1, a_2, a_3, loc.shape[0], loc[:,0], loc[:,1], loc[:,2], chg,
                3.0*r_d_hat**2*h_max, r_d_hat*h_max)
        ewald = -2.69595457432924945
        self.assertAlmostEqual(ene, ewald, places=9)

        # Si
        a_1 = np.array([7.25654832321381, 0.00000000000000, 0.00000000000000])
        a_2 = np.array([3.62827416160690, 6.28435519169252, 0.00000000000000])
        a_3 = np.array([3.62827416160690, 2.09478506389751, 5.92494689524090])
        loc = np.array([[0.0,  0.0,  0.0],
                        [0.25, 0.25, 0.25]])
        loc = (np.vstack((a_1, a_2, a_3)).T).dot(loc.T).T # to cartesian
        chg = 4.0 * np.ones(loc.shape[0])
        h_max = 5.92
        r_d_hat = 2.0
        ene_per_ion = real_space_electrostatic_sum.energy(
                a_1, a_2, a_3, loc.shape[0], loc[:,0], loc[:,1], loc[:,2], chg,
                3.0*r_d_hat**2*h_max, r_d_hat*h_max) / loc.shape[0]
        ewald_per_ion = -8.39857465282205418/loc.shape[0]
        self.assertAlmostEqual(ene_per_ion, ewald_per_ion, places=9)

        # SiO2
        a_1 = np.array([ 9.28422445623683, 0.00000000000000, 0.00000000000000])
        a_2 = np.array([-4.64211222811842, 8.04037423353787, 0.00000000000000])
        a_3 = np.array([ 0.00000000000000, 0.00000000000000, 10.2139697101486])
        loc = np.array([[0.41500, 0.27200, 0.21300],
                        [0.72800, 0.14300, 0.54633],
                        [0.85700, 0.58500, 0.87967],
                        [0.27200, 0.41500, 0.78700],
                        [0.14300, 0.72800, 0.45367],
                        [0.58500, 0.85700, 0.12033],
                        [0.46500, 0.00000, 0.33333],
                        [0.00000, 0.46500, 0.66667],
                        [0.53500, 0.53500, 0.00000]])
        loc = (np.vstack((a_1, a_2, a_3)).T).dot(loc.T).T # to cartesian
        chg = 6.0 * np.ones(loc.shape[0]) # most are O
        chg[6:] = 4.0                     # three are Si
        h_max = 10.21
        r_d_hat = 2.0
        ene_per_ion = real_space_electrostatic_sum.energy(
                a_1, a_2, a_3, loc.shape[0], loc[:,0], loc[:,1], loc[:,2], chg,
                3.0*r_d_hat**2*h_max, r_d_hat*h_max) / loc.shape[0]
        ewald_per_ion = -69.48809871723248932 / loc.shape[0]
        self.assertAlmostEqual(ene_per_ion, ewald_per_ion, places=9)

        # Al2SiO5
        a_1 = np.array([14.7289033699982, 0.00000000000000, 0.00000000000000])
        a_2 = np.array([0.00000000000000, 14.9260018049230, 0.00000000000000])
        a_3 = np.array([0.00000000000000, 0.00000000000000, 10.5049875335275])
        loc = np.array([[0.23030, 0.13430, 0.23900],
                        [0.76970, 0.86570, 0.23900],
                        [0.26970, 0.63430, 0.26100],
                        [0.73030, 0.36570, 0.26100],
                        [0.76970, 0.86570, 0.76100],
                        [0.23030, 0.13430, 0.76100],
                        [0.73030, 0.36570, 0.73900],
                        [0.26970, 0.63430, 0.73900],
                        [0.00000, 0.00000, 0.24220],
                        [0.50000, 0.50000, 0.25780],
                        [0.00000, 0.00000, 0.75780],
                        [0.50000, 0.50000, 0.74220],
                        [0.37080, 0.13870, 0.50000],
                        [0.42320, 0.36270, 0.50000],
                        [0.62920, 0.86130, 0.50000],
                        [0.57680, 0.63730, 0.50000],
                        [0.12920, 0.63870, 0.00000],
                        [0.07680, 0.86270, 0.00000],
                        [0.87080, 0.36130, 0.00000],
                        [0.92320, 0.13730, 0.00000],
                        [0.24620, 0.25290, 0.00000],
                        [0.42400, 0.36290, 0.00000],
                        [0.10380, 0.40130, 0.00000],
                        [0.75380, 0.74710, 0.00000],
                        [0.57600, 0.63710, 0.00000],
                        [0.89620, 0.59870, 0.00000],
                        [0.25380, 0.75290, 0.50000],
                        [0.07600, 0.86290, 0.50000],
                        [0.39620, 0.90130, 0.50000],
                        [0.74620, 0.24710, 0.50000],
                        [0.92400, 0.13710, 0.50000],
                        [0.60380, 0.09870, 0.50000]])
        loc = (np.vstack((a_1, a_2, a_3)).T).dot(loc.T).T # to cartesian
        chg = 6.0 * np.ones(loc.shape[0]) # most are O
        chg[8:13]  = 3.0                  # eight are Al
        chg[14] = 3.0
        chg[16] = 3.0
        chg[18] = 3.0
        chg[20] = 4.0                     # four are Si
        chg[23] = 4.0
        chg[26] = 4.0
        chg[29] = 4.0
        h_max = 14.93
        r_d_hat = 2.0
        ene_per_ion = real_space_electrostatic_sum.energy(
                a_1, a_2, a_3, loc.shape[0], loc[:,0], loc[:,1], loc[:,2], chg,
                3.0*r_d_hat**2*h_max, r_d_hat*h_max) / loc.shape[0]
        ewald_per_ion = -244.05500850908111943 / loc.shape[0]
        self.assertAlmostEqual(ene_per_ion, ewald_per_ion, places=9)

    def test_force(self):

        # SiO2
        a_1 = np.array([ 9.28422445623683, 0.00000000000000, 0.00000000000000])
        a_2 = np.array([-4.64211222811842, 8.04037423353787, 0.00000000000000])
        a_3 = np.array([ 0.00000000000000, 0.00000000000000, 10.2139697101486])
        loc = np.array([[0.41500, 0.27200, 0.21300],
                        [0.72800, 0.14300, 0.54633],
                        [0.85700, 0.58500, 0.87967],
                        [0.27200, 0.41500, 0.78700],
                        [0.14300, 0.72800, 0.45367],
                        [0.58500, 0.85700, 0.12033],
                        [0.46500, 0.00000, 0.33333],
                        [0.00000, 0.46500, 0.66667],
                        [0.53500, 0.53500, 0.00000]])
        loc = (np.vstack((a_1, a_2, a_3)).T).dot(loc.T).T # to cartesian
        chg = 6.0 * np.ones(loc.shape[0]) # most are O
        chg[6:] = 4.0                     # three are Si
        h_max = 10.21
        r_d_hat = 2.0

        # compute forces
        fx, fy, fz = real_space_electrostatic_sum.force(
                a_1, a_2, a_3, loc.shape[0], loc[:,0], loc[:,1], loc[:,2], chg,
                3.0*r_d_hat**2*h_max, r_d_hat*h_max)
        f = np.vstack((fx, fy, fz)).T

        # compare with finite-difference-derived forces
        d = 1e-4
        a = np.vstack((a_1, a_2, a_3)).T
        ai = np.linalg.inv(a)
        for i in range(loc.shape[0]): # loop over atoms
            for j in range(3):        # loop over x/y/z

                # add tiny amount to coordinate
                loc[i,j] += d
                x = ai.dot(loc[i,:]) # get scaled position
                x = x - np.floor(x)  # ensure position is within cell
                loc[i,:] = a.dot(x)  # recover cartesian position
                ene_p = real_space_electrostatic_sum.energy(
                        a_1, a_2, a_3, loc.shape[0], loc[:,0], loc[:,1], loc[:,2], chg,
                        3.0*r_d_hat**2*h_max, r_d_hat*h_max)

                # subtract tiny amount from coordinate 
                loc[i,j] -= 2*d
                x = ai.dot(loc[i,:]) # get scaled position
                x = x - np.floor(x)  # ensure position is within cell
                loc[i,:] = a.dot(x)  # recover cartesian position
                ene_m = real_space_electrostatic_sum.energy(
                        a_1, a_2, a_3, loc.shape[0], loc[:,0], loc[:,1], loc[:,2], chg,
                        3.0*r_d_hat**2*h_max, r_d_hat*h_max)

                # compute approximate force and compare
                f_finite_difference = -(ene_p-ene_m)/(2*d)
                self.assertAlmostEqual(f[i,j], f_finite_difference, places=6)

                # restore coordinate to original
                loc[i,j] += d
                x = ai.dot(loc[i,:]) # get scaled position
                x = x - np.floor(x)  # ensure position is within cell
                loc[i,:] = a.dot(x)  # recover cartesian position

    def test_stress(self):

        # define initial strain (applied to generate nontrivial stresses)
        t = np.eye(3) + np.array([[-0.25,  0.35, -0.15],
                                  [ 0.35,  0.15,  0.25],
                                  [-0.15,  0.25, -0.20]])

        # strained SiO2
        a_1 = t.dot([ 9.28422445623683, 0.00000000000000, 0.00000000000000])
        a_2 = t.dot([-4.64211222811842, 8.04037423353787, 0.00000000000000])
        a_3 = t.dot([ 0.00000000000000, 0.00000000000000, 10.2139697101486])
        loc = np.array([[0.41500, 0.27200, 0.21300],
                        [0.72800, 0.14300, 0.54633],
                        [0.85700, 0.58500, 0.87967],
                        [0.27200, 0.41500, 0.78700],
                        [0.14300, 0.72800, 0.45367],
                        [0.58500, 0.85700, 0.12033],
                        [0.46500, 0.00000, 0.33333],
                        [0.00000, 0.46500, 0.66667],
                        [0.53500, 0.53500, 0.00000]])
        loc = (np.vstack((a_1, a_2, a_3)).T).dot(loc.T).T # to cartesian
        chg = 6.0 * np.ones(loc.shape[0]) # most are O
        chg[6:] = 4.0                     # three are Si
        h_max = 10.21
        r_d_hat = 2.0

        # compute stress
        stress = real_space_electrostatic_sum.stress(
                a_1, a_2, a_3, loc.shape[0], loc[:,0], loc[:,1], loc[:,2], chg,
                3.0*r_d_hat**2*h_max, r_d_hat*h_max)

        # compute stress numerically
        d = 1e-4
        for i in range(6):
            eps = np.zeros([3,3])

            # nudge up by tiny amount
            if i==0 or i==1 or i==2:
                eps[i,i] = d
            elif i==3:
                eps[1,2] = d
            elif i==4:
                eps[0,2] = d
            elif i==5:
                eps[0,1] = d
            a_e = (np.eye(3) + eps).dot(np.vstack((a_1, a_2, a_3)).T)
            loc_e = (np.eye(3) + eps).dot(loc.T).T
            ene_p = real_space_electrostatic_sum.energy(
                    a_e[:,0], a_e[:,1], a_e[:,2], loc.shape[0], loc_e[:,0], loc_e[:,1], loc_e[:,2], chg,
                    3.0*r_d_hat**2*h_max, r_d_hat*h_max)

            # nudge down by tiny amount
            eps = -eps
            a_e = (np.eye(3) + eps).dot(np.vstack((a_1, a_2, a_3)).T)
            loc_e = (np.eye(3) + eps).dot(loc.T).T
            ene_m = real_space_electrostatic_sum.energy(
                    a_e[:,0], a_e[:,1], a_e[:,2], loc.shape[0], loc_e[:,0], loc_e[:,1], loc_e[:,2], chg,
                    3.0*r_d_hat**2*h_max, r_d_hat*h_max)

            # compute approximate stress and compare
            volume = np.linalg.det(np.vstack((a_1, a_2, a_3)).T)
            s_finite_difference = 1.0/volume*(ene_p-ene_m)/(2*d)
            self.assertAlmostEqual(stress[i], s_finite_difference, places=8)

if __name__ == '__main__':
    unittest.main()
