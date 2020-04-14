! MIT License
! 
! Copyright (c) 2019-2020 William C. Witt
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

module c_real_space_electrostatic_sum

    use iso_c_binding, only: c_double, c_int
    use real_space_electrostatic_sum

    implicit none

contains

subroutine c_real_space_electrostatic_sum_energy(&
        a1, a2, a3, n, rx, ry, rz, z, rc, rd, e) bind(c)
!______________________________________________________________________________
!
    implicit none

    real(c_double), intent(in)   ::  a1(3), a2(3), a3(3)
    integer(c_int), intent(in)   ::  n
    real(c_double), intent(in)   ::  rx(n), ry(n), rz(n)
    real(c_double), intent(in)   ::  z(n)
    real(c_double), intent(in)   ::  rc
    real(c_double), intent(in)   ::  rd
    real(c_double), intent(out)  ::  e
!______________________________________________________________________________
!
    call energy(a1, a2, a3, n, rx, ry, rz, z, rc, rd, e)

end subroutine

subroutine c_real_space_electrostatic_sum_force(&
        a1, a2, a3, n, rx, ry, rz, z, rc, rd, fx, fy, fz) bind(c)
!______________________________________________________________________________
!
    implicit none

    real(c_double), intent(in)   ::  a1(3), a2(3), a3(3)
    integer(c_int), intent(in)   ::  n
    real(c_double), intent(in)   ::  rx(n), ry(n), rz(n)
    real(c_double), intent(in)   ::  z(n)
    real(c_double), intent(in)   ::  rc
    real(c_double), intent(in)   ::  rd
    real(c_double), intent(out)  ::  fx(n), fy(n), fz(n)
!______________________________________________________________________________
!
    call force(a1, a2, a3, n, rx, ry, rz, z, rc, rd, fx, fy, fz)

end subroutine

subroutine c_real_space_electrostatic_sum_stress(&
        a1, a2, a3, n, rx, ry, rz, z, rc, rd, s) bind(c)
!______________________________________________________________________________
!
    implicit none

    real(c_double), intent(in)   ::  a1(3), a2(3), a3(3)
    integer(c_int), intent(in)   ::  n
    real(c_double), intent(in)   ::  rx(n), ry(n), rz(n)
    real(c_double), intent(in)   ::  z(n)
    real(c_double), intent(in)   ::  rc
    real(c_double), intent(in)   ::  rd
    real(c_double), intent(out)  ::  s(6)
!______________________________________________________________________________
!
    call stress(a1, a2, a3, n, rx, ry, rz, z, rc, rd, s)

end subroutine

end module
