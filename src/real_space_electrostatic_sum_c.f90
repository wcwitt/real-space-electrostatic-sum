! MIT License
! 
! Copyright (c) 2019 William C. Witt
! 
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
! 
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

module real_space_electrostatic_sum_c

    use iso_c_binding, only: c_double, c_int
    use real_space_electrostatic_sum

    implicit none

contains

subroutine energy_c(a_x, a_y, a_z, num, loc, chg, r_c, r_d, ene) bind(c, name='energy_c')
!______________________________________________________________________________
!
    implicit none

    real(c_double),  intent(in)   ::  a_x(3), a_y(3), a_z(3)
    integer(c_int),  intent(in)   ::  num
    real(c_double),  intent(in)   ::  loc(num,3)
    real(c_double),  intent(in)   ::  chg(num)
    real(c_double),  intent(in)   ::  r_c
    real(c_double),  intent(in)   ::  r_d
    real(c_double),  intent(out)  ::  ene
!______________________________________________________________________________
!
    call energy(a_x, a_y, a_z, num, loc, chg, r_c, r_d, ene)

end subroutine

end module
