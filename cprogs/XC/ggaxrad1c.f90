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
!                      Kurt R. Glaesemann
!
! Universidad de Madrid - Jose Ortega
!
!                  made possible by support from Motorola
!
!
! RESTRICTED RIGHTS LEGEND
!
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! ggaxrad1c.f90
! Program Description
! ===========================================================================
!
!      This routine calculates the exchange potential and energy density.
! Spherical symmetry is used. LSDA - GGA
!
! input
!    mode = 1    LSDA
!    mode = 2    GGA-X Becke
!    mode = 3    GGA-X Perdew
!    mode = 5    GGA-X Burke-Perdew-Ernzerhof
!
! convention
!    yy(1) = spin up, yy(2) = spin down
!
! Martin Fuchs, FHI der MPG, Berlin, 07-1992
!
! ===========================================================================
! Code rewritten to FORTRAN 90 by:
! James P. Lewis
! Campus Box 7260
! Department of Biochemistry and Biophysics
! University of North Carolina
! Chapel Hill, NC 27599-7260
! FAX 919-966-2852
! Office telephone 919-966-4644
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine ggaxrad1c (mode, rin, rho, rhop, rhopp, xpot, xen)
        use precision
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: mode

        real(kind=long), intent (in) :: rin

        real(kind=long), intent (in), dimension (2) :: rho
        real(kind=long), intent (in), dimension (2) :: rhop
        real(kind=long), intent (in), dimension (2) :: rhopp

! Output
        real(kind=long), intent (out) :: xen

        real(kind=long), intent (out), dimension (2) :: xpot

! Local Parameters and Data Declaration
! ===========================================================================
        real(kind=long), parameter :: eps = 1.0d-15

! Local Variable Declaration and Description
! ===========================================================================
        integer ispin

        real(kind=long) density
        real(kind=long) densityp
        real(kind=long) densitypp
        real(kind=long) ex
        real(kind=long) fermik
        real(kind=long) pi
        real(kind=long) r
        real(kind=long) s
        real(kind=long) u
        real(kind=long) v
        real(kind=long) vx

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Initialize
        pi = 3.141592653589793238462643D0

! If r is really small, then set to manageably small number.
        r = rin
        if (rin .lt. 1.0d-4) r = 1.0d-4

! exchange GGA, loop for up & down spin
        xen = 0.0d0
        do ispin = 1, 2
         if (rho(ispin) .le. eps) then
          xpot(ispin) = 0.0d0
         else
          density = 2.0d0*rho(ispin)
          if (mode .eq. 1) then
           call xlda (density, vx, ex)
          else if (mode .eq. 2 .or. mode .eq. 3 .or. mode .eq. 5) then
           densityp = 2.0d0*rhop(ispin)
           densitypp = 2.0d0*rhopp(ispin)
           fermik = 2.0d0*(3.0d0*pi*pi*density)**(1.0d0/3.0d0)

! s = abs(grad d)/(2kf*d)
! u = (grad d)*grad(abs(grad d))/(d**2 * (2*kf)**3)
!      >>  grad(abs(grad d) has mixed derivatives ! <<
! v = (laplacian d)/(d*(2*kf)**2)
           s = abs(densityp)/(fermik*density)
           u = abs(densityp)*densitypp/(density*density*fermik**3)
           v = (densitypp + 2.0d0*densityp/r)/(density*fermik*fermik)

           select case (mode)
            case (2)
             call xbecke (density, s, u, v, ex, vx)
            case (3)
             call exch (density, s, u, v, ex, vx)
            case (5)
             call exchpbe (density, s, u, v, 1, 1, ex, vx)
           end select
          else
           stop 'ggaxrad1c : mode improper'
          end if
          xpot(ispin) = vx
          xen = xen + rho(ispin)*ex
         end if
        end do

! energy
        xen = xen/max(rho(1) + rho(2), eps)

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================

        return
        end
