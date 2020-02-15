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

module real_space_electrostatic_sum

    implicit none

    integer,  parameter  ::  dp = 8
    real(dp), parameter  ::  pi = 3.14159265358979323846264338327950288419_dp
    real(dp), parameter  ::  sqrt_pi = sqrt(pi)
    real(dp), parameter  ::  one_third = 1.0_dp/3.0_dp

contains

subroutine energy(a1, a2, a3, num, loc, z, rc, rd, ene)
!______________________________________________________________________________
!
    implicit none

    real(dp), intent(in)   ::  a1(3), a2(3), a3(3)
    integer,  intent(in)   ::  num
    real(dp), intent(in)   ::  loc(num,3)
    real(dp), intent(in)   ::  z(num)
    real(dp), intent(in)   ::  rc
    real(dp), intent(in)   ::  rd
    real(dp), intent(out)  ::  ene

    real(dp) ::  vol, rho, ei, qi, rij, a(3,3), bt(3,3), &
                 d_100, d_010, d_001, origin_j(3), xyz_ij(3), ra
    integer  ::  i, j, shift1, shift2, shift3, &
                 shift1max, shift2max, shift3max
!______________________________________________________________________________
!
    ! compute cell volume and average density
    vol = a1(1) * (a2(2)*a3(3) - a2(3)*a3(2)) + &
          a1(2) * (a2(3)*a3(1) - a2(1)*a3(3)) + &
          a1(3) * (a2(1)*a3(2) - a2(2)*a3(1))
    vol = abs(vol) ! for left-handed coordinate systems
    rho = sum(z) / vol

    ! compute reciprocal lattice vectors
    a(:,1) = a1;  a(:,2) = a2;  a(:,3) = a3
    call invert_3x3(a, bt)  ! note: bt still missing factor of 2*pi

    ! compute distances between planes (accounts for missing 2*pi in bt)
    d_100 = 1.0_dp / sqrt(sum(bt(1,:) * bt(1,:)))
    d_010 = 1.0_dp / sqrt(sum(bt(2,:) * bt(2,:)))
    d_001 = 1.0_dp / sqrt(sum(bt(3,:) * bt(3,:)))

    ! compute the number of cells to include along each direction
    shift1max = ceiling(rc / d_100)
    shift2max = ceiling(rc / d_010)
    shift3max = ceiling(rc / d_001)

    ! prepare for loop over ions in cell
    ene = 0.0_dp

    ! loop over ions in cell
    do i = 1, num
    
        ! prepare for loop over neighboring ions
        ei = 0.0_dp
        qi = z(i)  ! b/c the i==j part of the sum is skipped below

        ! loop over cells
        do shift3 = -shift3max, shift3max
        do shift2 = -shift2max, shift2max
        do shift1 = -shift1max, shift1max

            ! get origin of shifted cell and loop over ions in that cell
            origin_j = shift1*a1 + shift2*a2 + shift3*a3            
            do j = 1, num

                ! exclude i==j
                if (i==j .and. shift1==0 .and. shift2==0 .and. shift3==0) cycle

                ! compute distance between points
                xyz_ij = loc(i,:) - (origin_j + loc(j,:))
                rij = sqrt(sum(xyz_ij * xyz_ij))

                ! only proceed if rij < rc
                if (rij > rc) cycle

                ! update energy and charge
                ei = ei + z(j) * erfc(rij/rd) / rij
                qi = qi + z(j)

            end do  ! j

        end do  ! shift1
        end do  ! shift2
        end do  ! shift3

        ! apply 1/2 z(i) factor to energy
        ei = 0.5_dp * z(i) * ei
    
        ! add correction terms
        ra = (3.0_dp * qi / (4.0_dp * pi * rho))**one_third  ! adaptive cutoff
        ei = ei - pi * z(i) * rho * ra * ra  &
              + pi * z(i) * rho * (ra*ra - rd*rd/2.0_dp) * erf(ra/rd)  &
              + sqrt_pi * z(i) * rho * ra * rd * exp(-ra*ra/(rd*rd)) &
              - 1.0/(sqrt_pi * rd) * z(i) * z(i)

        ! increment the total energy
        ene = ene + ei

    end do  ! i

end subroutine

subroutine force(a1, a2, a3, num, loc, z, rc, rd, fx, fy, fz)
!______________________________________________________________________________
!
    implicit none

    real(dp), intent(in)   ::  a1(3), a2(3), a3(3)
    integer,  intent(in)   ::  num
    real(dp), intent(in)   ::  loc(num,3)
    real(dp), intent(in)   ::  z(num)
    real(dp), intent(in)   ::  rc
    real(dp), intent(in)   ::  rd
    real(dp), intent(out)  ::  fx(num), fy(num), fz(num)

    real(dp) ::  rij, rijrij, rij_rd, a(3,3), bt(3,3), &
                 d_100, d_010, d_001, origin_j(3), xyz_ij(3), t
    integer  ::  i, j, shift1, shift2, shift3, &
                 shift1max, shift2max, shift3max
!______________________________________________________________________________
!
    ! compute reciprocal lattice vectors
    a(:,1) = a1;  a(:,2) = a2;  a(:,3) = a3
    call invert_3x3(a, bt)  ! note: bt still missing factor of 2*pi

    ! compute distances between planes (accounts for missing 2*pi in bt)
    d_100 = 1.0_dp / sqrt(sum(bt(1,:) * bt(1,:)))
    d_010 = 1.0_dp / sqrt(sum(bt(2,:) * bt(2,:)))
    d_001 = 1.0_dp / sqrt(sum(bt(3,:) * bt(3,:)))

    ! compute the number of cells to include along each direction
    shift1max = ceiling(rc / d_100)
    shift2max = ceiling(rc / d_010)
    shift3max = ceiling(rc / d_001)

    ! prepare for loop over ions in cell
    fx = 0.0_dp
    fy = 0.0_dp
    fz = 0.0_dp

    ! loop over ions in cell
    do i = 1, num
    
        ! loop over cells
        do shift3 = -shift3max, shift3max
        do shift2 = -shift2max, shift2max
        do shift1 = -shift1max, shift1max

            ! get origin of shifted cell and loop over ions in that cell
            origin_j = shift1*a1 + shift2*a2 + shift3*a3            
            do j = 1, num

                ! exclude i==j
                if (i==j .and. shift1==0 .and. shift2==0 .and. shift3==0) cycle

                ! compute distance between points
                xyz_ij = loc(i,:) - (origin_j + loc(j,:))
                rijrij = sum(xyz_ij * xyz_ij)
                rij = sqrt(rijrij)
                rij_rd = rij / rd

                ! only proceed if rij < rc
                if (rij > rc) cycle

                ! compute common term
                t = z(j) * & 
                    (2.0_dp / sqrt_pi * rij_rd * exp(-rij_rd * rij_rd) &
                        + erfc(rij_rd)) / (rij * rijrij)

                ! add contributions to forces
                fx(i) = fx(i) + t * xyz_ij(1)
                fy(i) = fy(i) + t * xyz_ij(2)
                fz(i) = fz(i) + t * xyz_ij(3)

            end do  ! j

        end do  ! shift1
        end do  ! shift2
        end do  ! shift3

        ! apply z(i) factor
        fx(i) = z(i) * fx(i)
        fy(i) = z(i) * fy(i)
        fz(i) = z(i) * fz(i)

    end do  ! i

end subroutine

subroutine stress(a1, a2, a3, num, loc, z, rc, rd, s)
!______________________________________________________________________________
!
    implicit none

    real(dp), intent(in)   ::  a1(3), a2(3), a3(3)
    integer,  intent(in)   ::  num
    real(dp), intent(in)   ::  loc(num,3)
    real(dp), intent(in)   ::  z(num)
    real(dp), intent(in)   ::  rc
    real(dp), intent(in)   ::  rd
    real(dp), intent(out)  ::  s(6)

    real(dp) ::  vol, rho, qi, rij, rijrij, rij_rd, a(3,3), bt(3,3), &
                 d_100, d_010, d_001, origin_j(3), xyz_ij(3), t, si(6), &
                 ra, ra_rd
    integer  ::  i, j, shift1, shift2, shift3, &
                 shift1max, shift2max, shift3max
!______________________________________________________________________________
!
    ! compute cell volume and average density
    vol = a1(1) * (a2(2)*a3(3) - a2(3)*a3(2)) + &
          a1(2) * (a2(3)*a3(1) - a2(1)*a3(3)) + &
          a1(3) * (a2(1)*a3(2) - a2(2)*a3(1))
    vol = abs(vol) ! for left-handed coordinate systems
    rho = sum(z) / vol

    ! compute reciprocal lattice vectors
    a(:,1) = a1;  a(:,2) = a2;  a(:,3) = a3
    call invert_3x3(a, bt)  ! note: bt still missing factor of 2*pi

    ! compute distances between planes (accounts for missing 2*pi in bt)
    d_100 = 1.0_dp / sqrt(sum(bt(1,:) * bt(1,:)))
    d_010 = 1.0_dp / sqrt(sum(bt(2,:) * bt(2,:)))
    d_001 = 1.0_dp / sqrt(sum(bt(3,:) * bt(3,:)))

    ! compute the number of cells to include along each direction
    shift1max = ceiling(rc / d_100)
    shift2max = ceiling(rc / d_010)
    shift3max = ceiling(rc / d_001)

    ! prepare for loop over ions in cell
    s = 0.0_dp

    ! loop over ions in cell
    do i = 1, num
    
        ! prepare for loop over neighboring ions
        si = 0.0_dp
        qi = z(i)  ! b/c the i==j part of the sum is skipped below
        
        ! loop over cells
        do shift3 = -shift3max, shift3max
        do shift2 = -shift2max, shift2max
        do shift1 = -shift1max, shift1max

            ! get origin of shifted cell and loop over ions in that cell
            origin_j = shift1*a1 + shift2*a2 + shift3*a3            
            do j = 1, num

                ! exclude i==j
                if (i==j .and. shift1==0 .and. shift2==0 .and. shift3==0) cycle

                ! compute distance between points
                xyz_ij = loc(i,:) - (origin_j + loc(j,:))
                rijrij = sum(xyz_ij * xyz_ij)
                rij = sqrt(rijrij)
                rij_rd = rij / rd

                ! only proceed if rij < rc
                if (rij > rc) cycle

                ! compute common term
                t = z(j) * & 
                    (2.0_dp / sqrt_pi * rij_rd * exp(-rij_rd * rij_rd) &
                        + erfc(rij_rd)) / (rij * rijrij)

                ! add contributions to stresses and charge
                si(1:3) = si(1:3) + t * xyz_ij(1:3) * xyz_ij(1:3)
                si(4) = si(4) + t * xyz_ij(2) * xyz_ij(3)
                si(5) = si(5) + t * xyz_ij(1) * xyz_ij(3)
                si(6) = si(6) + t * xyz_ij(1) * xyz_ij(2)
                qi = qi + z(j)

            end do  ! j

        end do  ! shift1
        end do  ! shift2
        end do  ! shift3

        ! apply factor of -z(i) / (2 * volume)
        si = -z(i) / (2.0_dp * vol) * si
    
        ! add remaining terms
        ra = (3.0_dp * qi / (4.0_dp * pi * rho))**one_third  ! adaptive cutoff
        ra_rd = ra / rd
        t = pi / (2.0_dp * vol) * rd * rd * rho * z(i) &
            * (1.0_dp - 2.0/sqrt_pi * ra_rd * exp(-ra_rd * ra_rd) &
                + (2.0_dp / 3.0_dp * ra_rd * ra_rd - 1.0_dp) * erfc(ra_rd))
        si(1:3) = si(1:3) + t

        ! increment the total stress
        s = s + si

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
