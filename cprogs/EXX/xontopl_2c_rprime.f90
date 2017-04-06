! copyright info:
!
!                             @Copyright 1999
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
!                      Kurt R. Glaesemann
!
! Universidad de Madrid - Jose Ortega
!
!                  made possible by support from Motorola
!                         fireball@fermi.la.asu.edu
 
!
! RESTRICTED RIGHTS LEGEND
!
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.
 
! xontopl_2c_rprime.f
! Program Description
! ===========================================================================
!       This routine calculates the integral over rprime, which is
! specifically used in the ontop (left) interactions for the exact exchange.
!
! ===========================================================================
! Code rewritten by:
! James P. Lewis
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ===========================================================================
 
! Program Declaration
! ===========================================================================
        subroutine xontopl_2c_rprime (nspec_max, nssh, nalpha, itype,      &
    &                                 rcutoff, nrho, lmax)
        use precision
        use x_exact
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: itype
        integer, intent (in) :: lmax
        integer, intent (in) :: nalpha
        integer, intent (in) :: nrho
        integer, intent (in) :: nspec_max
 
        integer, intent (in) :: nssh (nspec_max)
 
        real*8, intent (in) :: rcutoff
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer irho
        integer irhop
        integer issh
        integer lqn
 
        real*8 psi3
        real*8 psi4
        real*8 r
        real*8 rhomax
        real*8 rhomin
        real*8 rp
        real*8 sumrp
 
        real*8, dimension (:), allocatable :: rhopmult
 
        real*8, external :: psiofr
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! Set integration limits
        rhomin = 0.0d0
        rhomax = rcutoff
 
! Strictly define what the density of the mesh should be. 
        nnrhop = 2*nrho - 1
        drhop = rcutoff/real(nnrhop)

! Set up Simpson's rule factors. First for the z integration and then for
! the rho integration.
        allocate (rpoint(nnrhop))
        allocate (rhopmult(nnrhop))
        rpoint(1) = rhomin
        rpoint(nnrhop) = rhomax
        rhopmult(1) = drhop/3.0d0
        rhopmult(nnrhop) = drhop/3.0d0
        do irho = 2, nnrhop - 1, 2
         rpoint(irho) = rhomin + real(irho - 1)*drhop
         rhopmult(irho) = 4.0d0*drhop/3.0d0
        end do
        do irho = 3, nnrhop - 2, 2
         rpoint(irho) = rhomin + real(irho - 1)*drhop
         rhopmult(irho) = 2.0d0*drhop/3.0d0
        end do
 
! Loop over all of the shells of itype.
        allocate (rprime(nnrhop, 0:2*lmax, nssh(itype)))
        do issh = 1, nssh(itype)
 
! Loop over all the possible quantum numbers.
         do lqn = 0, 2*lmax
 
! Perform the radial integration over r' for each given r.
          do irho = 1, 2*nrho - 1
           r = rpoint(irho)
           if (r .lt. 1.0d-4) r = 1.0d-4
 
           sumrp = 0.0d0
           do irhop = 1, 2*nrho - 1
            rp = rpoint(irhop)
            if (rp .lt. 1.0d-4) rp = 1.0d-4

            if (rp .lt. rcutoff) then
             psi3 = psiofr (itype, nalpha, rp)
             psi4 = psiofr (itype, issh, rp)
 
! Limits from 0 to r.
             if (rp .le. r) then
              sumrp =                                                        &
      &        sumrp + rhopmult(irhop)*psi3*psi4*rp**(lqn + 2)/r**(lqn + 1)
 
! Limits from r to rcutoff
             else
              sumrp = sumrp + rhopmult(irhop)*psi3*psi4*r**lqn/rp**(lqn - 1)
             end if
            end if
           end do
 
! The actual integral as a function of r.
           rprime(irho, lqn, issh) = sumrp
          end do
         end do
        end do
 
! Deallocate Arrays
! ===========================================================================
        deallocate (rhopmult)
 
! Format Statements
! ===========================================================================
 
        return
        end
