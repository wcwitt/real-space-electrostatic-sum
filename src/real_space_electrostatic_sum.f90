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

module real_space_electrostatic_sum

    implicit none

    integer,  parameter  ::  dp = 8
    real(dp), parameter  ::  pi = 3.14159265358979323846264338327950288419_dp
    real(dp), parameter  ::  sqrt_pi = sqrt(pi)
    real(dp), parameter  ::  one_third = 1.0_dp/3.0_dp

contains

subroutine energy(a_x, a_y, a_z, num, loc, chg, r_c, r_d, ene)
!______________________________________________________________________________
!
    implicit none

    real(dp), intent(in)   ::  a_x(3), a_y(3), a_z(3)
    integer,  intent(in)   ::  num
    real(dp), intent(in)   ::  loc(num,3)
    real(dp), intent(in)   ::  chg(num)
    real(dp), intent(in)   ::  r_c
    real(dp), intent(in)   ::  r_d
    real(dp), intent(out)  ::  ene

    real(dp) ::  vol, rho, e_i, q_i, r_ij, a(3,3), b_t(3,3), &
                 d_100, d_010, d_001, origin_j(3), xyz_ij(3), r_a
    integer  ::  i, j, shift_x, shift_y, shift_z, &
                 shift_x_max, shift_y_max, shift_z_max
!______________________________________________________________________________
!
    ! compute cell volume and average density
    vol = a_x(1) * (a_y(2)*a_z(3) - a_y(3)*a_z(2)) + &
          a_x(2) * (a_y(3)*a_z(1) - a_y(1)*a_z(3)) + &
          a_x(3) * (a_y(1)*a_z(2) - a_y(2)*a_z(1))
    vol = abs(vol) ! for left-handed coordinate systems
    rho = sum(chg) / vol

    ! compute reciprocal lattice vectors
    a(:,1) = a_x;  a(:,2) = a_y;  a(:,3) = a_z
    call invert_3x3(a, b_t)  ! note: b_t still missing factor of 2*pi

    ! compute distances between planes (accounts for missing 2*pi in b_t)
    d_100 = 1.0_dp / sqrt(b_t(1,1)*b_t(1,1) + b_t(1,2)*b_t(1,2) + b_t(1,3)*b_t(1,3))
    d_010 = 1.0_dp / sqrt(b_t(2,1)*b_t(2,1) + b_t(2,2)*b_t(2,2) + b_t(2,3)*b_t(2,3))
    d_001 = 1.0_dp / sqrt(b_t(3,1)*b_t(3,1) + b_t(3,2)*b_t(3,2) + b_t(3,3)*b_t(3,3))

    ! compute the number of cells to include along each direction
    shift_x_max = ceiling(r_c / d_100)
    shift_y_max = ceiling(r_c / d_010)
    shift_z_max = ceiling(r_c / d_001)

    ! loop over atoms in cell
    ene = 0.0_dp
    do i = 1, num
    
        ! loop to compute the sum in Eq. (13)
        e_i = 0.0_dp
        q_i = chg(i)  ! b/c the i==j part of the sum is skipped below
        do shift_z = -shift_z_max, shift_z_max
        do shift_y = -shift_y_max, shift_y_max
        do shift_x = -shift_x_max, shift_x_max

            ! get origin of shifted cell and loop over atoms in that cell
            origin_j = shift_x*a_x + shift_y*a_y + shift_z*a_z            
            do j = 1, num

                ! exclude i==j
                if (i==j .and. shift_x==0 .and. shift_y==0 .and. shift_z==0) cycle

                ! compute distance between atoms
                xyz_ij = loc(i,:) - (origin_j + loc(j,:))
                r_ij = sqrt(sum(xyz_ij * xyz_ij))

                ! only proceed if r_ij < r_c
                if (r_ij > r_c) cycle
                e_i = e_i + chg(j) * erfc(r_ij/r_d) / r_ij  ! mult. by chg(i) below
                q_i = q_i + chg(j)

            end do  ! j

        end do  ! shift_x
        end do  ! shift_y
        end do  ! shift_z
        e_i = 0.5_dp * chg(i) * e_i
    
        ! compute the correction terms in Eq. (19)
        r_a = (3.0_dp*q_i/(4.0_dp*pi*rho))**one_third  ! adaptive cutoff
        e_i = e_i - pi * chg(i) * rho * r_a * r_a  &
                  + pi * chg(i) * rho * (r_a*r_a - r_d*r_d/2.0_dp) * erf(r_a/r_d)  &
                  + sqrt_pi * chg(i) * rho * r_a * r_d * exp(-r_a*r_a/(r_d*r_d)) &
                  - 1.0/(sqrt_pi * r_d) * chg(i) * chg(i)

        ! increment the total energy
        ene = ene + e_i

    end do  ! i

end subroutine

subroutine invert_3x3(a, b)
!______________________________________________________________________________
!
!   inverts a 3x3 matrix. not hyper-optimized.
!______________________________________________________________________________
!
    implicit none

    real(dp), intent(in)   ::  a(3,3)
    real(dp), intent(out)  ::  b(3,3)

    real(dp)  ::  det
!______________________________________________________________________________
!
    ! compute cofactor matrix
    b(1,1) = a(2,2)*a(3,3) - a(3,2)*a(2,3)
    b(2,1) = a(2,3)*a(3,1) - a(3,3)*a(2,1)
    b(3,1) = a(2,1)*a(3,2) - a(3,1)*a(2,2)
    b(1,2) = a(1,3)*a(3,2) - a(3,3)*a(1,2)
    b(2,2) = a(1,1)*a(3,3) - a(3,1)*a(1,3)
    b(3,2) = a(1,2)*a(3,1) - a(3,2)*a(1,1)
    b(1,3) = a(1,2)*a(2,3) - a(2,2)*a(1,3)
    b(2,3) = a(1,3)*a(2,1) - a(2,3)*a(1,1)
    b(3,3) = a(1,1)*a(2,2) - a(2,1)*a(1,2)

    ! compute determinant and divide
    det = a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(2,3))
    det = det - a(1,2)*(a(2,1)*a(3,3) - a(3,1)*a(2,3))
    det = det + a(1,3)*(a(2,1)*a(3,2) - a(3,1)*a(2,2))
    b = b / det

end subroutine

end module
