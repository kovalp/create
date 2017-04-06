! copyright info:
!
!                             @Copyright 1998
!                          Fireball2000 Committee
!
! ASU - Otto F. Sankey
!       Kevin Schmidt
!       Jian Jun Dong
!       John Tomfohr
!       Gary B. Adams
!
! Motorola - Alex A. Demkov
!
! University of Regensburg - Juergen Fritsch
!
! University of Utah - James P. Lewis
!                      Richard B. Evans
!                      Kurt R. Glaesemann
!
! Universidad de Madrid - Jose Ortega
!
! Texas A&M - Traian Dumitrica
!
!                  made possible by support from Motorola
!                         fireball@fermi.la.asu.edu

!
! fireball-qmd is a free (GPLv3) open project.

!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.


! dvxc.f
! Program Description
! ===========================================================================
!
! This function computes vxc(n1+n2) - vxc(n1) with the densities ni of atom
! in_i at the distance r_i from their centers.
!
! On input:  in1,in2: atomic indices
!            r, z: geometry information for the charge gradient
!            ix: switch for the derivatives
!
! On output:  vxc: vxc[n1(r1) + n2(r2)] - vxc[n1(r1)]

! We calculate vxc(n1+n2) - vxc(n1). The catch comes in when we compute
! derivatives. We compute neutral, neutral for ideriv1. For other ideriv's we
! have the following KEY:
!
! (xy) means charge on (1,2). Case 1 (KEY=1),
! neutral neutral corresponds to (00) etc.
! KEY = 1,2,3,4,5 for ideriv=1,2,3,4,5
!
!                             dq Atom 2 axis
!                                   |
!                                   + (0+) KEY=5
!                                   |
!                                   |
!                                   |
!                       KEY=2       |  KEY=1
!                      (-0)         |(00)        (+0) KEY=3
!                     --+-----------O-----------+--
!                                   |           dq Atom 1 axis
!                                   |
!                                   |
!                                   |
!                                   |
!                                   +(0-) KEY=4
!                                   |
!
! ===========================================================================
! Original code from Juergen Fritsch

! Code rewritten by:
! James P. Lewis and Richard B. Evans
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        real*8 function dvxc (in1, in2, r, z, r1, iexc, fraction, ix)
        implicit none

        include '../parameters.inc'
        include '../wavefunctions.inc'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer iexc
        integer in1, in2
        integer ix

        real*8 fraction
        real*8 r
        real*8 r1
        real*8 z

! Local Parameters and Data Declaration
! ===========================================================================
        real*8 abohr
        parameter (abohr = 0.529177249d0)

! Local Variable Declaration and Description
! ===========================================================================
        integer in3

        real*8 dens
        real*8 densp
        real*8 denspp
        real*8 densz
        real*8 denszz
        real*8 denspz
        real*8 dexc1c
        real*8 dexc2c
        real*8 dnuxc
        real*8 dnuxcs
        real*8 dnuxc2c
        real*8 dnuxcs2c
        real*8 dvxc1c
        real*8 dvxc2c
        real*8 hartree
        real*8 rin
        real*8 dexc

! Procedure
! ===========================================================================
! By default in3 = in2
        in3 = in2

! Two-center piece: vxc[n1 + n2(r,z)]
! ***************************************************************************
! Interpolate the density and gradients of the density at the given
! point (r, z).
        call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz,
     1                      nnrho, nrho_points, nnz, nz_points, rho2c,
     2                      dens)

! Only interpolate the derivatives if doing GGA exchange-correlation.
        if (iexc .eq. 4 .or. iexc .eq. 5 .or. iexc .eq. 6
     1      .or. iexc .eq. 9 .or. iexc .eq. 10) then
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz,
     1                       nnrho, nrho_points, nnz, nz_points, rhop2c,
     2                       densp)
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz,
     1                       nnrho, nrho_points, nnz, nz_points,
     2                       rhopp2c, denspp)
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz,
     1                       nnrho, nrho_points, nnz, nz_points,
     2                       rhopz2c, denspz)
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz,
     1                       nnrho, nrho_points, nnz, nz_points,
     2                       rhoz2c, densz)
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz,
     1                       nnrho, nrho_points, nnz, nz_points,
     2                       rhozz2c, denszz)
        end if

! Convert to atomic units
        rin = r/abohr
        dens = dens*abohr**3
        densp = densp*abohr**4
        densz = densz*abohr**4
        denspp = denspp*abohr**5
        denspz = denspz*abohr**5
        denszz = denszz*abohr**5

! Here energy and potential due to exchange and correlation are calculated.
        call get_potxc2c (iexc, fraction, rin, dens, densp, denspp,
     1                    densz, denszz, denspz, dexc2c, dvxc2c,
     2                    dnuxc2c, dnuxcs2c)

! One-center piece: vxc[n1(r1)]
! Compute the exchange correlation potential for the one-center case
! ***************************************************************************
! Evaluate the density for the one-center - in1
        call density_calc (iexc, ix, 1, in1, in2, in3, r1, drho,
     1                     dens, densp, denspp)

        rin = r1/abohr
        dens = dens*abohr**3
        densp = densp*abohr**4
        denspp = denspp*abohr**5
        call get_potxc1c (iexc, fraction, rin, dens, densp, denspp,
     1                    dexc1c, dvxc1c, dnuxc, dnuxcs, dexc)

! Answers are in Hartrees convert to eV.
        hartree = 14.39975d0/abohr
        dvxc = hartree*(dvxc2c - dvxc1c)

! Format Statements
! ===========================================================================

        return
        end
