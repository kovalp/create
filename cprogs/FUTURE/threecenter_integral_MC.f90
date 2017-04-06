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
!            Jun Wang
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
 
! fireball-qmd is a free (GPLv3) open project.

!
! This program is free software: you can redistribute it and/or modify
! by the Government is subject to restrictions as set forth in
! subdivsion { (b) (3) (ii) } of the Rights in Technical Data and
! Computer Software clause at 52.227-7013.
 
 
!     ============================================================
       subroutine threecenter_integral(dbcx,rna,rc1,rc2,nrr,nt,nphi,    &
     &                   in1,in2,in3,gmat,index_max,n1,l1,m1,n2,l2,m2,  &
     &                   iexc,interaction,ispmin,ispmax)
!     ============================================================
!
!     Calculates the matrix elements of the neutral atom for the
!     specific configuration of dbcx (dist. between bc) and rna=position of
!     neutral atom.
!
!     ------------------------------------------------------------
!
      implicit none
!
      include '../parameters.inc'
!
!     ------------------------------------------------------------
!
!     Passed variables:
!
      integer, intent(in) :: in1,in2,in3    ! types of the atoms
      integer, intent(in) :: iexc           ! exchange-correlation approximation
      integer, intent(in) :: nrr,nt,nphi    ! number of integration points
!
      integer, intent(in) :: index_max      ! maximal number of matrix elements
      integer, intent(in) :: interaction    ! number for the interaction
      integer, intent(in) :: ispmin,ispmax  ! range for isorp loop
!
      integer, intent(in) :: n1 (inter_max) ! shell number of left  atom
      integer, intent(in) :: n2 (inter_max) !                 right atom
      integer, intent(in) :: l1 (inter_max) ! angular momentum of that shell
      integer, intent(in) :: l2 (inter_max) !
      integer, intent(in) :: m1 (inter_max) ! m-value in that shell
      integer, intent(in) :: m2 (inter_max) !
!
      real*8, intent(in) :: rc1           ! i=1,2,3  for left,right, neutral atm
      real*8, intent(in) :: rc2
 
      real*8, intent(in) :: dbcx          ! bond charge distance  (A)
      real*8, intent(in) :: rna(3)        ! neutral atom location (A)
!
      real*8, intent(out) :: gmat(0:10,inter_max)   ! result from the integrator
!
!     External function calls
      real*8 psiofr,vnnaofr,dvxc3c
      external dvxc3c, vnnaofr, psiofr
 
!
!     ............................................................
!
!     interaction = 1: bcna  neutral atom
!                   2: xc3c  exchange correlation
!
!     for the exchange correaltion case,
!
!     gmat(ix,*) stores the derivatives w.r.t charges:
!
!     gmat(1,*) : in1,0    gmat(2,*) : in1,-q1  gmat(3,*) : in1,+q1
!     gmat(4,*) : in2,-q2  gmat(5,*) : in2,+q2
!     gmat(6,*) : in3,-q3  gmat(7,*) : in3,+q3
!
!     ............................................................
!
!     ------------------------------------------------------------
!
!     Internal variables:
!
!
      real*8 vpot(0:10)
 
      real*8 znormMU(inter_max),znormNU(inter_max),   &
    &        psiAmat(inter_max),psiBmat(inter_max)
!
!     dimensions up to f-orbitals
      real*8 thfactMU(0:3,-3:3),thfactNU(0:3,-3:3),phifactor(-3:3)
!
!
      integer nn1,nl1,nm1,nn2,nl2,nm2,ir,ix,ix1,inm,isorpX,npoints
!
      real*8  sq3,sq15,pi,r,theta,phi,rmax,dc1,ds1,zr,r2,dc2,ds2,xr,yr,r3,  &
     &        cphi,sphi,stuffmunu,rmin,dr,dphi,dtheta,dvolume
      real*8 sq12,sq4
 
      real*8, allocatable, dimension(:,:), save :: data
!
! =====================================================================
!     The size of the matrix is determined by Nsh(in1) and Nsh(in2)
!     We do only those matrix elements that are not zero
!     The number of the non-zero matrix elements is INMAX
!
!     Normalization factors
      sq3=sqrt(3.0d0)
      sq15=sqrt(15.0d0)
      sq12=sqrt(12.d0)
      sq4 =sqrt(4.d0)
      pi=3.141592653589793238462643D0
 
      npoints = nrr*nt*nphi
 
!
!     Open MC integration file - if not already done
!
      if(.not. allocated(data))then
        allocate (data(1:3,1:npoints))
        open (unit = 99, file = 'Halton.dat', status = 'old')
        do inm=1,npoints
          read(99,*,END=2000,ERR=1000) data(1,inm),data(2,inm),data(3,inm)
        end do
        close(99)
      end if
 
! ===================================================================
! Here is the correct list:
!   m:     -2       -1         0        1         2
!          xy       yz      3z^2-r^2   xz       x^2-y^2
! sq15 * (  1        1        1/sq12    1       1/sq4  )
!
! ====================================================================
      do 313 inm=1,index_max
        nl1=l1(inm)
        nl2=l2(inm)
        if(nl1.eq.0)znormMU(inm)=1.0d0
        if(nl2.eq.0)znormNU(inm)=1.0d0
        if(nl1.eq.1)znormMU(inm)=sq3
        if(nl2.eq.1)znormNU(inm)=sq3
        if(nl1.eq.2)znormMU(inm)=sq15
        if(nl2.eq.2)znormNU(inm)=sq15
! First we calculate the m values.
        nm1=m1(inm)
        nm2=m2(inm)
! working on d for mu
        if(nl1.eq.2)then
          if(nm1.eq.0)znormMU(inm)=znormMU(inm)/sq12
          if(nm1.eq.2)znormMU(inm)=znormMU(inm)/sq4
        end if
! working on d for nu
        if(nl2.eq.2)then
          if(nm2.eq.0)znormNU(inm)=znormNU(inm)/sq12
          if(nm2.eq.2)znormNU(inm)=znormNU(inm)/sq4
        end if
 313  continue
 
!
!     isorp=0,1,2,3, for neutral atom, s part, p part, d part.
!     We add 1 because of 0 being the neutral atom.
!
      rmin = 0.0d0
      rmax = dmax1(rc1,rc2)
!
 
      do inm=1,index_max
        do isorpX=0,10
          gmat(isorpX,inm)=0.0d0
        end do
      end do
!
! ========================================================
!                 Do integral over r, theta, phi:
! ========================================================
      dr=rmax-rmin
      dtheta=pi
      dphi=2.0d0*pi
 
      do 310 ir=1,npoints    ! npoints is the total number of points in space
!       r=0 to rmax
!       theta=0 to pi
!       phi=0 to 2pi
        r    =data(1,ir)*rmax
        theta=data(2,ir)*pi
        phi  =data(3,ir)*pi*2
 
!       Theta stuff for atom 1
        dc1=cos(theta)
        ds1=sin(theta)
!
        r2=sqrt(r*r+dbcx*dbcx-2*r*dbcx*dc1)
!
        if(r2.gt.rc2)go to 310       ! outside integration range
 
        do inm=1,index_max
          nn1=n1(inm)
          psiAmat(inm)=psiofr(in1,nn1,r)
        end do
        do inm=1,index_max
          nn2 = n2(inm)
          psiBmat(inm)=psiofr(in2,nn2,r2)
        enddo
!
        zr=r*dc1
!
!       Theta stuff for atom 2.
!       Be careful for r2 very small.
!
        if(r2.gt.0.00001)then
          dc2=(zr-dbcx)/r2
          ds2=ds1*r/r2
        else
          dc2=1.0d0
          ds2=0.0d0
        end if
 
        cphi=cos(phi)
        sphi=sin(phi)
 
!
!       Theta factors for A and B.
!
!       Use the (l,m) notation. For example 3z^2-1 becomes
!       (2,0), px becomes (1,1), and so on.
!
!       ATOM A .........................................
!
!       S
        thfactMU(0,0)=1.0d0
!
!       P
!       Note: We order the orbitals here x,y,z (or pi,pi',sig)
!
!       watch it: theta factors for d-orbitals should
!                 also contain the m-dependency of the
!                 prefactors znormMU and znormNU.
!                 This is not the case so far.
!
        thfactMU(1,1)=ds1
        thfactMU(1,-1)=ds1
        thfactMU(1,0)=dc1
!
!       D Order of d-orbitals is 3z^2-1, x^2-y^2, xz, xy, yz
!
        thfactMU(2,0)=3.0d0*dc1*dc1-1.0d0
        thfactMU(2,2)=ds1*ds1
        thfactMU(2,1)=ds1*dc1
        thfactMU(2,-2)=ds1*ds1
        thfactMU(2,-1)=ds1*dc1
!
!       ATOM B .............................................
!
!       S
        thfactNU(0,0)=1.0d0
!
!       P
!
        thfactNU(1,1)=ds2
        thfactNU(1,-1)=ds2
        thfactNU(1,0)=dc2
!
!       D Order of d-orbitals is 3z^2-1, x^2-y^2, xz, xy, yz
!
        thfactNU(2,0)=3.0d0*dc2*dc2-1.0d0
        thfactNU(2,2)=ds2*ds2
        thfactNU(2,1)=ds2*dc2
        thfactNU(2,-2)=ds2*ds2
        thfactNU(2,-1)=ds2*dc2
!
!       Done with theta factors.
!
!       Do the phi integral, with the phifactors.
!       Note: We order the p-orbitals
!       here x,y,z (or pi,pi',sig), NOT z,x,y.
!       Note that px, and xz now are +1. And so on!
!
!       The phi factors depend only on m.
 
        phifactor(0)=1.0d0
 
        phifactor(1)=cphi
        phifactor(-1)=sphi
! d's
        phifactor(2)=cphi*cphi-sphi*sphi
        phifactor(-2)=cphi*sphi
 
        xr=r*ds1*cphi
        yr=r*ds1*sphi
 
        r3=sqrt((xr-rna(1))**2+(yr-rna(2))**2+(zr-rna(3))**2)
!
! ---------------------------------------------------
!
        do  iX=ispmin,ispmax
!
          IF(interaction .EQ. 1) vpot(iX)=vnnaofr(in3,iX,r3)
 
          IF(interaction .EQ. 2) then
           IX1=IX+1
           vpot(iX)=dvxc3c(iexc, r, r2, r3, in1, in2, in3, IX1)
          END IF
!
        end do
 
        do inm=1,index_max
          nl1=l1(inm)
          nm1=m1(inm)
          nl2=l2(inm)
          nm2=m2(inm)
!         average over theta (divide by 2)
!         average over phi   (divide by 2*pi)
          stuffmunu=0.5d0*(0.5d0/pi)*r*r*ds1*              &
     &              thfactMU(nl1,nm1)*thfactNU(nl2,nm2)*   &
     &              phifactor(nm1)*phifactor(nm2)*         &
     &              psiBmat(inm)*psiAmat(inm)
          do ix=ispmin,ispmax
            gmat(ix,inm)=gmat(ix,inm)+vpot(ix)*stuffmunu
          end do
        end do
 
 310  continue
!
!     Finally, the normalization factors for the different orbitals.
!     For instance a p orbital is sqrt(3) * x/r * r(r). The sqrt(3)
!     factor (and ones like it) are now included.
!
      dvolume=dr*dtheta*dphi/npoints
 
      do inm=1,index_max
        do ix=ispmin,ispmax
           gmat(ix,inm)=gmat(ix,inm)*znormMU(inm)*znormNU(inm)*dvolume
        end do
      end do
 
!
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
!
       return
 1000  write(*,*) 'An error accured while reading Halton.dat'
       stop
 2000  write(*,*) 'Halton.dat file does not have enough points'
       stop
       end
 
