! copyright info:
!
!                             @Copyright 2004
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Universidad de Madrid - Pavel Jelinek

! Other contributors, past and present:
! Auburn University - Jianjun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Brigham Young University - Hao Wang
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio University - Dave Drabold
! University of Regensburg - Juergen Fritsch

!
! fireball-qmd is a free (GPLv3) open project.

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

 
! threecenter_integral.f90
! Program Description
! ===========================================================================
!       This routine performs the actual three-center integration for all
! types of three center interactions.  A specific configuration of dbcx
! (distance between the bond charge) and rna (the position of the neutral
! atom) are input, and the integral is performed for this configuration.  

! ...........................................................................
!     interaction = 1: bcna  neutral atom
!                   2: xc3c  exchange correlation
!                   3: xc3c  average densities (SNXC and OLSXC)
!
!     for the exchange correlation case,
!
!     gmat(ix,*) stores the derivatives w.r.t charges:
!
!     gmat(1,*) : in1,  0   gmat(2,*) : in1, -q1  gmat(3,*) : in1, +q1
!     gmat(4,*) : in2, -q2  gmat(5,*) : in2, +q2
!     gmat(6,*) : in3, -q3  gmat(7,*) : in3, +q3
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
        subroutine threecenter_integral (in1, in2, in3, index_max, iexc,     &
     &                                   interaction, ispmin, ispmax,        & 
     &                                   ispherical, dbcx, rna, rcutoff1,    &
     &                                   rcutoff2, nr, ntheta, nphi, n1, l1, &
     &                                   m1, n2, l2, m2, gmat)
        use dimensions
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: in1
        integer, intent (in) :: in2
        integer, intent (in) :: in3
        integer, intent (in) :: iexc         ! type of exchange-correlation
        integer, intent (in) :: index_max    ! maximum number of matrix elements
        integer, intent (in) :: interaction  ! number for the interaction
        integer, intent (in) :: ispmax
        integer, intent (in) :: ispmin
        integer, intent (in) :: nphi
        integer, intent (in) :: nr
        integer, intent (in) :: ntheta

        integer, intent (in), dimension (inter_max) :: n1  ! left atom shell
        integer, intent (in), dimension (inter_max) :: n2  ! right atom shell
        integer, intent (in), dimension (inter_max) :: l1  ! angular momentum
        integer, intent (in), dimension (inter_max) :: l2
        integer, intent (in), dimension (inter_max) :: m1  ! m-value in shell
        integer, intent (in), dimension (inter_max) :: m2

        real*8, intent (in) :: dbcx
        real*8, intent (in) :: rcutoff1  ! largest radius of i-th atom (in Ang.)
        real*8, intent (in) :: rcutoff2  ! 1, 2, 3 = left, right, neutral atom

        real*8, dimension (3) :: rna

        logical, intent (in) :: ispherical

! Output
        real*8, dimension (0:10, inter_max) :: gmat

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer index
        integer ip
        integer ir
        integer isorp
        integer it
        integer ix
        integer lleft
        integer mleft
        integer nleft
        integer lright
        integer mright
        integer nright

        real*8 dphi
        real*8 dr
        real*8 dtheta
        real*8 pi
        real*8 r
        real*8 rmax
        real*8 rmin
      
        real*8, dimension (:), allocatable :: cphiy
        real*8, dimension (:), allocatable :: hrho
        real*8, dimension (:), allocatable :: htheta
        real*8, dimension (:), allocatable :: hphi
        real*8, dimension (:), allocatable :: wrho
        real*8, dimension (-3:3) :: phifactor
        real*8, dimension (inter_max) :: psiAmat
        real*8, dimension (inter_max) :: psiBmat
        real*8, dimension (:), allocatable :: sphiy
        real*8, dimension (0:10, inter_max) :: thetamat
        real*8, dimension (0:3,-3:3) :: thfactMU  ! up to f-orbitals
        real*8, dimension (0:3,-3:3) :: thfactNU
        real*8, dimension (:), allocatable :: wtheta
        real*8, dimension (:), allocatable :: wphi
        real*8, dimension (inter_max) :: znormMU
        real*8, dimension (inter_max) :: znormNU

      real*8 h
      integer tempR, tempT, tempP
      integer tmpImid,imid
      real*8 xmin
      real*8 xxp
      real*8 psipsi
      real*8 psiofr,vpot,vnnaofr,dvxc3c
      external dvxc3c, vnnaofr,psiofr
      real*8 avgVmat(0:10,inter_max), fsimp(5000)
      real*8 w1, w2
      integer ix1
      real*8  theta,phi,dc1,ds1,dc,ds,zr,r2,dc2,ds2,xr,yr,r3,  &
     &        cphi,sphi,simp,simp2,simpson,  &
     &        averagephi,averagetheta,  &
     &        prod,prod2,dsth,stuffmunu

! Allocate arrays
! ===========================================================================
        allocate (cphiy (nphi))
        allocate (hrho (nr))
        allocate (hphi (nphi))
        allocate (htheta (ntheta))
        allocate (sphiy (nphi))
        allocate (wrho (nr))
        allocate (wphi (nphi))
        allocate (wtheta (ntheta))

! Procedure
! ===========================================================================
! The size of the matrix is determined by nssh(in1) and nssh(in2).  We do 
! only those matrix elements that are not zero.  The number of the non-zero 
! matrix elements is INMAX
      
! Initialize the integration array:
        gmat (ispmin:ispmax, 1:index_max) = 0.0d0

! Initialize integration factors and other variables. 
        pi = 4.0d0*atan(1.0d0)
        rmin = 0.0d0
        rmax = max (rcutoff1, rcutoff2)

        dr = (rmax - rmin)/float(nr - 1)
        dtheta = pi/float(ntheta - 1)
        dphi = 2.0d0*pi/float(nphi - 1)

        hrho(1) = rmin
        wrho(1) = dr*41.0d0/140.0d0 
        do ir = 2, nr - 1
         hrho(ir) = rmin + (ir - 1.0d0)*dr
         if (mod(ir,6) .eq. 2) wrho(ir) = dr*216.0d0/140.0d0 
         if (mod(ir,6) .eq. 3) wrho(ir) = dr*27.0d0/140.0d0
         if (mod(ir,6) .eq. 4) wrho(ir) = dr*272.0d0/140.0d0
         if (mod(ir,6) .eq. 5) wrho(ir) = dr*27.0d0/140.0d0
         if (mod(ir,6) .eq. 0) wrho(ir) = dr*216.0d0/140.0d0
         if (mod(ir,6) .eq. 1) wrho(ir) = dr*82.0d0/140.0d0
        end do
        wrho(nr) = dr*41.0d0/140.0d0
        hrho(nr) = rmax
    
        htheta(1) = 0.0d0
        wtheta(1) = dtheta*41.0d0/140.0d0
        do it = 2, ntheta - 1
         htheta(it) = (it - 1.0d0)*dtheta
         if (mod(it,6) .eq. 2) wtheta(it) = dtheta*216.0d0/140.0d0
         if (mod(it,6) .eq. 3) wtheta(it) = dtheta*27.0d0/140.0d0
         if (mod(it,6) .eq. 4) wtheta(it) = dtheta*272.0d0/140.0d0
         if (mod(it,6) .eq. 5) wtheta(it) = dtheta*27.0d0/140.0d0
         if (mod(it,6) .eq. 0) wtheta(it) = dtheta*216.0d0/140.0d0
         if (mod(it,6) .eq. 1) wtheta(it) = dtheta*82.0d0/140.0d0
        end do
        wtheta(ntheta) = dtheta*41.0d0/140.0d0
        htheta(ntheta) = pi

        hphi(1) = 0.0d0
        wphi(1) = dphi*41.0d0/140.0d0
        do ip = 2, nphi - 1
         hphi(ip) = (ip - 1.0d0)*dphi
         if (mod(ip,6) .eq. 2) wphi(ip) = dphi*216.0d0/140.0d0
         if (mod(ip,6) .eq. 3) wphi(ip) = dphi*27.0d0/140.0d0
         if (mod(ip,6) .eq. 4) wphi(ip) = dphi*272.0d0/140.0d0
         if (mod(ip,6) .eq. 5) wphi(ip) = dphi*27.0d0/140.0d0
         if (mod(ip,6) .eq. 0) wphi(ip) = dphi*216.0d0/140.0d0
         if (mod(ip,6) .eq. 1) wphi(ip) = dphi*82.0d0/140.0d0
        end do
        wphi(nphi) = dphi*41.0d0/140.0d0
        hphi(nphi) = 2.0d0*pi

! ===================================================================
! Here is the correct list:
!   m:     -2       -1         0        1         2
!          xy       yz      3z^2-r^2   xz       x^2-y^2
! sq15 * (  1        1        1/sq12    1       1/2  )
! ====================================================================
        do index = 1, index_max
         lleft = l1(index)
         lright = l2(index)
         if (lleft .eq. 0) znormMU(index) = 1.0d0
         if (lright .eq. 0) znormNU(index) = 1.0d0
         if (lleft .eq. 1) znormMU(index) = sqrt(3.0d0)
         if (lright .eq. 1) znormNU(index) = sqrt(3.0d0)
         if (lleft .eq. 2) znormMU(index) = sqrt(15.0d0)
         if (lright .eq. 2) znormNU(index) = sqrt(15.0d0)

! First we calculate the m values.
         mleft = m1(index)
         mright = m2(index)
         if (lleft .eq. 2) then
           if (mleft .eq. 0) znormMU(index) = znormMU(index)/sqrt(12.0d0)
           if (mleft .eq. 2) znormMU(index) = znormMU(index)/2.0d0
         end if
         if (lright .eq. 2) then
           if (mright .eq. 0) znormNU(index) = znormNU(index)/sqrt(12.0d0)
           if (mright .eq. 2) znormNU(index) = znormNU(index)/2.0d0
         end if
        end do
 
! Note that isorp = 0, 1, 2, 3 for neutral atom, s part, p part, d part.
! We add 1 because of 0 being the neutral atom.
!
! Set up some constants.
! The phi integration does not depend on the r integration.
! Therefore, it is done outside the r loop.
        do ip = 1, nphi
         phi = hphi(ip)
         cphiy(ip) = cos(phi)
         sphiy(ip) = sin(phi)
        end do
     
! ----------------------------------------------------------------------------
! Integration over r
! ----------------------------------------------------------------------------
        do ir = 1, nr
         r = hrho(ir)

! Zero out the array for theta integral.
         thetamat(ispmin:ispmax, index_max) = 0.0d0

!        nli, nni tells us which shell is used for each atom.
!
         do index=1,index_max
           nleft=n1(index)
           psiAmat(index)=psiofr(in1,nleft,r)
           if(ispherical) psiAmat(index) = sqrt(psiAmat(index)**2.0d0)
         end do

 
!
! ========================================================
!              Do integral over theta:
! ========================================================
!
         do 208 it=1, ntheta 
!           
           theta=htheta(it)
!          Theta stuff for atom 1
           dc1=cos(theta)
           ds1=sin(theta)
           dc=dc1
           ds=ds1
!
           r2 = r**2 + dbcx**2 - 2*r*dbcx*dc
           if (r2 .le. 0.0) then
            r2 = 0.0d0
           else
            r2 = sqrt(r2)
           end if
!
           if(r2.gt.rcutoff2)go to 208       ! outside integration range
!                                       ! other theta values might work
!
           do index=1,index_max
             nright = n2(index)
             psiBmat(index)=psiofr(in2,nright,r2)
             if(ispherical) psiBmat(index) = sqrt(psiBmat(index)**2.0d0)
           enddo
!
           zr=r*dc
!
!          Theta stuff for atom 2.
!          Be careful for r2 very small.
!          Find cos(theta2), sin(theta2).
           if(r2.gt.0.00001)then
             dc2=(zr-dbcx)/r2
             ds2=ds1*r/r2
           else
             dc2=1.0d0
             ds2=0.0d0
           end if
!
!
! -------------------------------------------------------
!          Theta factors for A and B.
! -------------------------------------------------------
!          Use the (l,m) notation. For example 3z^2-1 becomes
!          (2,0), px becomes (1,1), and so on.
!
!          ATOM A .........................................
!
!          S
           thfactMU(0,0)=1.0d0
!
!          P
!          Note: We order the orbitals here x,y,z (or pi,pi',sig)
!
!          watch it: theta factors for d-orbitals should
!                    also contain the m-dependency of the
!                    prefactors znormMU and znormNU.
!                    This is not the case so far.
!
           thfactMU(1,1)=ds1
           thfactMU(1,-1)=ds1
           thfactMU(1,0)=dc1
!
!          D Order of d-orbitals is 3z^2-1, x^2-y^2, xz, xy, yz
!
           thfactMU(2,0)=3.0d0*dc1*dc1-1.0d0
           thfactMU(2,2)=ds1*ds1
           thfactMU(2,1)=ds1*dc1
           thfactMU(2,-2)=ds1*ds1
           thfactMU(2,-1)=ds1*dc1
!
!          ATOM B .............................................
!
!          S
           thfactNU(0,0)=1.0d0
!
!          P
!
           thfactNU(1,1)=ds2
           thfactNU(1,-1)=ds2
           thfactNU(1,0)=dc2
!
!          D Order of d-orbitals is 3z^2-1, x^2-y^2, xz, xy, yz
!
           thfactNU(2,0)=3.0d0*dc2*dc2-1.0d0
           thfactNU(2,2)=ds2*ds2
           thfactNU(2,1)=ds2*dc2
           thfactNU(2,-2)=ds2*ds2
           thfactNU(2,-1)=ds2*dc2
!
! --------------------------------------------------------------
!          Done with theta factors.
! --------------------------------------------------------------
 
           do index=1,index_max
             do ix=ispmin,ispmax
               avgVmat(ix,index)=0.0d0
             end do
           end do
!
! ========================================================
!          Do integral over phi:
! ========================================================
 
!          average over phi (divide by 2*pi)
           averagephi=0.5d0/pi
 
!          The phi factors depend only on m.
 
           phifactor(0)=1.0d0
 
           do 244 ip=1, nphi
!

             phi=hphi(ip)
             cphi=cphiy(ip)
             sphi=sphiy(ip)
!
!            Do the phi integral, with the phifactors.
!            Note: We order the p-orbitals
!            here x,y,z (or pi,pi',sig), NOT z,x,y.
!            Note that px, and xz now are +1. And so on!
!
             phifactor(1)=cphi
             phifactor(-1)=sphi
! d's
             phifactor(2)=cphi*cphi-sphi*sphi
             phifactor(-2)=cphi*sphi
 
             xr=r*ds*cphi
             yr=r*ds*sphi
!
             r3=sqrt((xr-rna(1))**2+(yr-rna(2))**2+(zr-rna(3))**2)
!
! ---------------------------------------------------
!
             do  iX=ispmin,ispmax
!
               IF(interaction .EQ. 1) vpot=vnnaofr(in3,iX,r3)
!
               IF(interaction .EQ. 2) then
                IX1=IX+1
                vpot = dvxc3c (iexc, r, r2, r3, in1, in2, in3, IX1)
               END IF
! xc3c_SN
               IF(interaction .EQ. 3) then
                  psipsi = psiofr(in3,ix,r3)
                  vpot=(psipsi**2)/(4.0d0*pi)
               ENDIF
!
!              note: dc,ds defined at the beginning of the theta loop
!
               prod=vpot*wphi(ip)*averagephi
               do index=1,index_max
                 mleft=m1(index)
                 mright=m2(index)
                 avgVmat(ix,index)=avgVmat(ix,index) +   &
     &                           prod*phifactor(mleft)*phifactor(mright)
               end do
             end do
!
! ===============================================
!          The end of the phi integral.
 244       continue
! ===============================================
!
           dsth=ds1
           averagetheta=0.5d0
           prod=wtheta(it)*averagetheta*dsth
 
           do 1407 index=1,index_max
             lleft=l1(index)
             mleft=m1(index)
             lright=l2(index)
             mright=m2(index)
             stuffmunu=prod*thfactMU(lleft,mleft)*      &
     &                      thfactNU(lright,mright)*psiBmat(index)
             do 1408 ix=ispmin,ispmax
               thetamat(ix,index)=thetamat(ix,index)+           &
     &                          avgVmat(ix,index)*stuffmunu
 1408        continue
 1407      continue
!
!
!
 208    continue
!       ========================================
!       The end of the integral over theta (loop 208).
!       =====================================================
 
!       now finish off the r integral!
 
        prod=wrho(ir)*r*r
 
        do index=1,index_max
          prod2=prod*psiAmat(index)
          do ix=ispmin,ispmax
            gmat(ix,index)=gmat(ix,index)+prod2*thetamat(ix,index)
          end do
        end do
 
! ----------------------------------------------------------------------------
! Integration over r - end
! ----------------------------------------------------------------------------
        end do
 
! Finally, the normalization factors for the different orbitals.  For instance,
! a p-orbital is sqrt(3)*x/r*psi(r).  The sqrt(3) factor (and ones like it) 
! are now included.
        do index = 1, index_max
         prod = znormMU(index)*znormNU(index)
         do ix = ispmin, ispmax
          gmat(ix,index) = gmat(ix,index)*prod
         end do
        end do
 
! ================================================================
!     SUMMARY
! We have computed gmat(isorp,mu,nu). mu, and nu are 1 to 9 for
! sp^3d^5.
! We are in molecular coordinates with z along sigma, x along pi, and
! y along pi'. Many matrix elements og gmat are zero. We have computed them
! anyway, and indeed find they are zero. Just to avoid at a later time,
! any trouble with roundoffs, I will now set those that are supposed to be
! zero, identically to zero.
! Here is the matrix, and the zero's are indicated.
!
!          \ s   x   y   z    3z^2-r^2    x^2-y^2    xz    xy      yz
!
!  s                 0                                     0        0
!
!  x                 0                                     0        0
!
!  y         0   0       0        0          0        0
!
!  z                 0                                     0        0
!
! 3z^2-r^2           0                                     0        0
!
! x^2-y^2            0                                     0        0
!
!  xz                0                                     0        0
!
!  xy        0   0       0        0          0        0
!
!  yz        0   0       0        0          0        0
!

! Deallocate Arrays 
! ===========================================================================
        deallocate (hrho, hphi, htheta)
        deallocate (wrho, wphi, wtheta)

! Format Statements
! ===========================================================================
        return
        end

