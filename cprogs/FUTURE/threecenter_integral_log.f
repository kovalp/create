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
 
! RESTRICTED RIGHTS LEGEND
!
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in
! subdivsion { (b) (3) (ii) } of the Rights in Technical Data and
! Computer Software clause at 52.227-7013.
 
 
!     ============================================================
       subroutine threecenter_integral(dbcx,rna,rc1,rc2,nrr,nt,nphi,
     1                   in1,in2,in3,gmat,index_max,n1,l1,m1,n2,l2,m2,
     2                   iexc,interaction,ispmin,ispmax)
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
      integer in1,in2,in3      ! types of the atoms
      integer iexc             ! exchange-correlation approximation
      integer nrr,nt,nphi      ! number of integration points
!
      integer index_max        ! maximal number of matrix elements
      integer interaction      ! number for the interaction
      integer ispmin,ispmax    ! range for isorp loop
!
      integer n1 (inter_max)   ! shell number of left  atom
      integer n2 (inter_max)   !                 right atom
      integer l1 (inter_max)   ! angular momentum of that shell
      integer l2 (inter_max)   !
      integer m1 (inter_max)   ! m-value in that shell
      integer m2 (inter_max)   !
!
      real*8 rc1           ! i=1,2,3  for left,right, neutral atm
      real*8 rc2
 
      real*8 dbcx          ! bond charge distance  (A)
      real*8 rna(3)        ! neutral atom location (A)
!
      real*8 gmat(0:10,inter_max)   ! result from the integrator
!
      real*8 psiofr,vpot,vnnaofr,dvxc3c
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
      real*8 thetamat(0:10,inter_max)
!
      real*8 znormMU(inter_max),znormNU(inter_max),
     &       psiAmat(inter_max),psiBmat(inter_max)
!
!
      real*8 thfactMU(0:3,-3:3),    !  dimensions up to f-orbitals
     &       thfactNU(0:3,-3:3),
     &       phifactor(-3:3)
!
      real*8 avgVmat(0:10,inter_max), fsimp(5000)
!
      integer   numbphi
      parameter (numbphi=5000)
      real*8    phiy(numbphi),
     &          cphiy(numbphi),sphiy(numbphi)
!
!
      integer nn1,nl1,nm1,nn2,nl2,nm2,i,ir,ix,ix1,it,ip,inm,
     &        nmax,isorpX
!
      real*8  sq3,sq15,pi,dr,dtheta,dphi,r,theta,phi,
     &        rmin,rmax,dc1,ds1,dc,ds,zr,r2,dc2,ds2,xr,yr,r3,
     &        cphi,sphi,simp,simp2,simpson,
     &        averagephi,averagetheta,
     &        prod,prod2,dsth,stuffmunu,thrd,ntinv,nphiinv,rloggy
      real*8 sq12,sq4
!
!
!     External Variables
!
      real*8 funky
      external funky
!     funky is defined at end of file.  This way a log, exp, or any other
!     function can be used to define the r grid
 
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
      thrd=1.0d0/3.0d0
      ntinv=1.0d0/dfloat(nt-1)
      nphiinv=1.0d0/dfloat(nphi-1)
 
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
 
      if(nphi.gt.numbphi)then
        write(*,*)' nphi=',nphi
        write(*,*)' numbphi=',numbphi
        write(*,*)' In threecenter_integral-- redimension numbphi'
        stop 'error in threecenter_integral'
      end if
!
!     =======================================================
!
!     Simpson rule factors 1,4,2,4,2,4,2,.....,2,4,1
!     (all divided by 3)
!     ***careful --- remember fsimp(ilast)=1.d0. But since
!     we dont know ilast, we must figure it out when needed.
!     We dont ALWAYS go to nmax (nrr,nt,nphi don't have to be the same)!
!
      nmax=max(nrr,nt,nphi)
!
      fsimp(1)=thrd
      do 5000 i=2,nmax-1
       if(mod(i,2).ne.0)then
        fsimp(i)=2.d0*thrd
       else
        fsimp(i)=4.d0*thrd
       end if
 5000 continue
      fsimp(nmax)=thrd
!
! ========================================================
!
      rmin=0.0d0
!     rmax=rc1
      rmax = dmax1(rc1,rc2)
!
      dtheta=pi*ntinv
      dphi=2.d0*pi*nphiinv
!
!     Set up some constants.
!     The phi integration does not depend on the r integration.
!     Therefore, it is done outside the r loop.
!
      do 7244 ip=1,nphi
        phiy(ip)=2.d0*pi*dfloat(ip-1)*nphiinv
        phi=phiy(ip)
        cphiy(ip)=cos(phi)
        sphiy(ip)=sin(phi)
 7244 continue
! ========================================================
 
      do 1451 inm=1,index_max
        do 1451 isorpX=ispmin,ispmax
          gmat(isorpX,inm)=0.0d0
 1451 continue
!
! ========================================================
!                 Do integral over r:
! ========================================================
!
 
      rloggy=(rmax-rmin)/funky(dfloat(nrr))
      do 310 ir=1,nrr
         r=rmin+funky(dfloat(ir))*rloggy
         if (ir .eq. 1) then
           dr=(funky(dfloat(ir)+0.5d0)-funky(dfloat(ir)      ))*rloggy
         else if (ir .eq. nrr) then
           dr=(funky(dfloat(ir)      )-funky(dfloat(ir)-0.5d0))*rloggy
         else
           dr=(funky(dfloat(ir)+0.5d0)-funky(dfloat(ir)-0.5d0))*rloggy
         end if
!
!        dorbital addition OFS
!        Fireball2000: Jose y Jandro.
!        Zero out the array for theta integral.
!
         do inm=1,index_max
           do ix=ispmin,ispmax
             thetamat(ix,inm)=0.0d0
           end do
         end do
!
!        nli, nni tells us which shell is used for each atom.
!
         do inm=1,index_max
           nn1=n1(inm)
           psiAmat(inm)=psiofr(in1,nn1,r)
         end do
 
!
! ========================================================
!              Do integral over theta:
! ========================================================
!
         do 208 it=1,nt
!
           theta=pi*dfloat(it-1)*ntinv
!          Theta stuff for atom 1
           dc1=cos(theta)
           ds1=sin(theta)
           dc=dc1
           ds=ds1
!
           r2=sqrt(r*r+dbcx*dbcx-2*r*dbcx*dc)
!
           if(r2.gt.rc2)go to 208       ! outside integration range
!                                       ! other theta values might work
!
           do inm=1,index_max
             nn2 = n2(inm)
             psiBmat(inm)=psiofr(in2,nn2,r2)
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
 
           simp=fsimp(it)
           if(it.eq.nt)simp=thrd
           do inm=1,index_max
             do ix=ispmin,ispmax
               avgVmat(ix,inm)=0.0d0
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
 
           do 244 ip=1,nphi
!
             simp2=fsimp(ip)
             if(ip .eq. nphi)simp2=thrd
             simpson=simp2*dphi
             phi=phiy(ip)
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
 
               IF(interaction .EQ. 2) then
                IX1=IX+1
                vpot = dvxc3c (iexc, r, r2, r3, in1, in2, in3, IX1)
               END IF
!
!              note: dc,ds defined at the beginning of the theta loop
!
               prod=vpot*simpson*averagephi
               do inm=1,index_max
                 nm1=m1(inm)
                 nm2=m2(inm)
                 avgVmat(ix,inm)=avgVmat(ix,inm) +
     1                           prod*phifactor(nm1)*phifactor(nm2)
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
           simpson=simp*dtheta
           prod=simpson*averagetheta*dsth
 
           do 1407 inm=1,index_max
             nl1=l1(inm)
             nm1=m1(inm)
             nl2=l2(inm)
             nm2=m2(inm)
             stuffmunu=prod*thfactMU(nl1,nm1)*
     &                      thfactNU(nl2,nm2)*psiBmat(inm)
             do 1408 ix=ispmin,ispmax
               thetamat(ix,inm)=thetamat(ix,inm)+
     1                          avgVmat(ix,inm)*stuffmunu
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
 
!       We don't use simpson's rule, since it's not an equally space grid
        prod=dr*r*r
 
        do inm=1,index_max
          prod2=prod*psiAmat(inm)
          do ix=ispmin,ispmax
            gmat(ix,inm)=gmat(ix,inm)+prod2*thetamat(ix,inm)
          end do
        end do
!
!
!
 310  continue
!
! ========================================================
!     The end of the integral over r (loop 310).
! ========================================================
!
!     Finally, the normalization factors for the different orbitals.
!     For instance a p orbital is sqrt(3) * x/r * r(r). The sqrt(3)
!     factor (and ones like it) are now included.
!
      do inm=1,index_max
        prod=znormMU(inm)*znormNU(inm)
        do ix=ispmin,ispmax
           gmat(ix,inm)=gmat(ix,inm)*prod
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
       end
 
       real*8 function funky(x)
       real*8 x
!      funky(1) must equal zero
       funky=exp(x)-exp(1.0D0)
       end
