! copyright info:
!
!                             @Copyright 2002
!                            Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio University - Dave Drabold
! University of Regensburg - Juergen Fritsch

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! get_potxc1c.f90
! Program Description
! ===========================================================================
!  This program will access the requested exchange and correlation
!  functionals from the following list:
!          1  LDA   Wigner
!          2  LDA   Hedin/Lundqvist
!          3  LDA   Ceperley/Alder Perdew/Zunger (1980)
!          4  GGA   Perdew/Wang (1991)
!          5  GGA   Becke (1988) X, Perdew (1986) C
!          6  GGA   Perdew/Burke/Ernzerhof (1996)
!          7  LDA   Zhao/Parr
!          8  LDA   Ceperley/Alder Perdew/Wang (1991)
!          9  GGA   Becke (1988) X, Lee/Yang/Parr (1988) C
!         10  GGA   Perdew/Wang (1991) X, Lee/Yang/Parr (1988) C
!         11  LSDA  Vosko/Wilk/Nusair (1980)
!         12  GGA   Becke (1988) X, Lee/Yang/Parr (1988) C
!                   with exact exchange
! The numerical value above is assigned to the variable iexc output
! exchange-correlation potential. This program has been modified to account
! for the fact that the density is a sum of two densities at two different
! centers.  Also, the potential is evaluated at one point in space and the
! integrals (matrix elements) are evaluated elsewhere.
!
! input
!    iexc        xc scheme
!    r           radial coordinate
!    rho         sum of atomic density in au
!    rhop        sum of atomic density gradient (with respect to r) in au
!    rhopp       sum of second gradient (with respect to r) in au
!
! output
!    vpxc        xc potential
!    newexc      xc energy
!    dnuxc
!    dnuxcs
! ===========================================================================
! Code written by:
! James P. Lewis
! Eduardo Mendez
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
!
! Program Declaration
! ===========================================================================
        subroutine get_potxc1c (iexc, fraction, r, rho, rhop, rhopp, newexc, &
     &                          vpxc, dnuxc, dnuxcs, dexc)
        use precision
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iexc

        real*8, intent(in) :: fraction
        real*8, intent(inout) :: r
        real*8, intent(inout) :: rho
        real*8, intent(in) :: rhop
        real*8, intent(in) :: rhopp

! Output
        real*8, intent(out) :: newexc
        real*8, intent(out) :: dexc
        real*8, intent(out) :: vpxc
        real*8, intent(out) :: dnuxc
        real*8, intent(out) :: dnuxcs

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer ix

        real*8 aln
        real*8 dec
        real*8 dex
        real*8 ecp
        real*8 ex
        real*8 exc
        real*8 fx
        real*8 fxc
        real*8 rh
        real*8 rs
        real*8 x
        real*8 zeta

        real*8, dimension (2) :: cpot
        real*8, dimension (2) :: d
        real*8, dimension (2) :: dp
        real*8, dimension (2) :: dpp
        real*8, dimension (2) :: xpot

! Procedure
! ===========================================================================
! Initialize to zero.
        dnuxc = 0.0d0
        dnuxcs = 0.0d0

! If r is really small, then set to manageably small number.
        if (r .lt. 1.0d-4) r = 1.0d-4

! Rho must be positive, but not too small
        if (rho .lt. 1.0d-8) then
         rho = 0.0d0
         dnuxc = 0.0d0
         newexc = 0.0d0
         vpxc = 0.0d0
         return
        else if (rho .lt. 1.0d-5) then
         rho = 1.0d-5
        end if

! Determine exchange-correlation potentials
! exchange (X) only
        if (iexc .eq. 11) then
         zeta = 0.0d0
         d(1) = rho*0.5*(1 + zeta)
         d(2) = rho*0.5*(1 - zeta)
         call lsdavwn (d, dex, dec, xpot, cpot, dnuxc, dnuxcs)
         newexc = dex + dec
         vpxc = xpot(2) + cpot(2)                ! Holds for zeta = 0.0d0 only

! correlation (C) Wigner
        else if (iexc .eq. 1) then
         rh = rho
         call wigner (rh, ex, fx, exc, fxc)
         vpxc = fxc
         dex = ex
         dec = exc - ex

! C Hedin-Lundqvist
        else if (iexc .eq. 2) then
         rh = rho
         if (rh .ne. 0.0d0) then
          rs = 0.62035049d0*rh**(-1.0d0/3.0d0)
          x = rs/21.0d0
          aln = log(1.0d0 + 1.0d0/x)
          ecp = aln + (x**3*aln - x*x) + x/2 - 1.0d0/3.0d0
          dex = -0.458175d0/rs - 0.0225d0*ecp
          vpxc = -0.6109d0/rs - 0.0225d0*aln
         else
          dex = 0.0d0
          vpxc = 0.0d0
         end if

! XC Ceperley - Alder
! as parameterized by Perdew and Zunger, Phys. Rev. B23, 5048 (1981)
        else if (iexc .eq. 3) then
         rh = rho
         call ceperley_alder (rh, ex, fx, exc, fxc, dexc, dnuxc)
         vpxc = fxc
         dex = ex
         dec = exc - ex

! X Perdew, C Perdew, generalized gradient approximation 1992
        else if (iexc .eq. 4) then
         rh = rho
         d = 0.5d0*rho
         dp = 0.5d0*rhop
         dpp = 0.5d0*rhopp
         call ggaxrad1c (3, r, d, dp, dpp, xpot, dex)
         call ggacrad1c (2, r, d, dp, dpp, cpot, dec)
         vpxc = xpot(1) + cpot(1)

! X Becke, C Perdew, generalized gradient approximation
        else if (iexc .eq. 5) then
         d = 0.5d0*rho
         dp = 0.5d0*rhop
         dpp = 0.5d0*rhopp
         call ggaxrad1c (2, r, d, dp, dpp, xpot, dex)
         call ggacrad1c (3, r, d, dp, dpp, cpot, dec)
         vpxc = xpot(1) + cpot(1)

! XC Wigner-scaled LDA of PRA 46, R5320 (1992)
        else if (iexc .eq. 7) then
         rh = rho
         call wigscaled(rh, ex, fx, exc, fxc)
         vpxc = fxc
         dex = ex
         dec = exc - ex

! XC Ceperley-Alder in Perdew-Wang parametrization of 1991
        else if (iexc .eq. 8) then
         d = 0.5d0*rho
         dp = 0.0d0
         dpp = 0.0d0
         call ggaxrad1c (1, r, d, dp, dpp, xpot, dex)
         call ggacrad1c (1, r, d, dp, dpp, cpot, dec)
         vpxc = xpot(1) + cpot(1)

! C Lee-Yang-Parr
        else if (iexc .eq. 9 .or. iexc .eq. 10 .or. iexc .eq. 12) then

! X Becke gga by default
         ix = 2

! X Perdew-Wang gga
         if (iexc .eq. 10) ix = 3
         d = 0.5d0*rho
         dp = 0.5d0*rhop
         dpp = 0.5d0*rhopp
         call ggaxrad1c (ix, r, d, dp, dpp, xpot, dex)
         call ggacrad1c (4, r, d, dp, dpp, cpot, dec)
         if (iexc .ne. 12) then
          vpxc = xpot(1) + cpot(1)
         else
          vpxc = (1.0d0 - fraction)*xpot(1) + cpot(1)
          dex = (1.0d0 - fraction)*dex
         end if

! XC burke-perdew-ernzerhof gga 1996
        else if(iexc .eq. 6) then
         d = 0.5d0*rho
         dp = 0.5d0*rhop
         dpp = 0.5d0*rhopp
         call ggaxrad1c (5, r, d, dp, dpp, xpot, dex)
         call ggacrad1c (5, r, d, dp, dpp, cpot, dec)
         vpxc = xpot(1) + cpot(1)

! If the improper iexc option was entered then the program will stop.
        else
         write (*,*) ' In get_potxc1c.f90 - '
         write (*,*) ' stop: xc option not implemented', iexc
         stop
        end if

! Calculate the exchange energy by combining the exchange and correlation
! energies and subtracting the exchange/correlation potential energy
! This comment seems to be old fashioned
        newexc = dec + dex

! Format Statements
! ===========================================================================
        return
        end
