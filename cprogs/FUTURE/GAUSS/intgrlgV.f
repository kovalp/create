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

! RESTRICTED RIGHTS LEGEND
!
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in
! subdivsion { (b) (3) (ii) } of the Rights in Technical Data and
! Computer Software clause at 52.227-7013.

!-----------------------------------------------------------------------
        subroutine intgrlgV(imqn,inorb,nz,nrho,d,rc,sum,
     &                      in3,itype,alpha,isorp,rc3)
!-----------------------------------------------------------------------
 
! Subroutine intgrl performs the integral over rho and z
! (the phi part having been done by hand and giving
! the "ifaktor"s below).
 
        implicit none
 
        integer iat,imqn,in3,inorb
        integer irho,isorp,itype,iz,nnrho,nnz,nrho,nz
        integer ifaktor1, ifaktor2,ifaktor3
 
        real*8 alpha,d,drho,dz,fact,faktor,pi,rc,rc3,rho
        real*8 rhomin,sum,vnnaofr,xmult,xxmult,ymult,yymult
        real*8 zmax,zmin
 
!       common /stuff8/iorb(2),idir(2)
        real*8 psi(2),z(2),r(2)
        dimension ifaktor1(20),ifaktor2(20),ifaktor3(20)
        dimension inorb(2)
 
        pi=3.141592653589793238462643D0
 
! These are the faktor(i) in the "MAGIC FORMULA":
!       faktor(i) = ifaktor1(i)*sqrt(ifaktor3(i))/ifaktor2(i)
 
        data ifaktor1/1,1,3,3,1,1,3,5,15,15,
     &                1,1,3,1,3,15,7,21,105,35/
        data ifaktor2/2,2,2,4,4,4,4,8,4,16,
     &                4,4,16,8,16,16,8,32,16,32/
        data ifaktor3/1,3,1,1,5,15,5,1,1,1,
     &                7,21,14,35,70,7,1,1,1,1/
 
! The factors of (4/3)drho (2/3)drho etc
! below are from simpson's rule.
! (see any basic numerical analysis text)
 
        sum=0.0d0
        if(abs(alpha).lt.1.e-10)return
 
        if (isorp.eq.0) then
         zmin=d-rc3
         zmax=5.0d0/sqrt(alpha)
         rc=rc3
        else
         zmin=-5.0d0/sqrt(alpha)
         zmax=5.0d0/sqrt(alpha)
         rc=5.0d0/sqrt(alpha)
        end if
 
!        zmin=d-rc3
!        zmax=rc
!        dz=( rc-(d-rc) )/dfloat(nz*2)
 
        dz=(zmax-zmin)/dfloat(nz*2)
 
        xxmult=2.0d0*dz/3.0d0
        rhomin=0.0d0
 
!       rhomax=rc
        drho=rc/dfloat(nrho*2)
        yymult=2.0d0*drho/3.0d0
        nnz=2*nz+1
        nnrho=2*nrho+1
 
        faktor = dfloat(ifaktor1(itype))/dfloat(ifaktor2(itype))*
     &           sqrt(dfloat(ifaktor3(itype)))
 
! new faktor factors for gV.
 
        do iat=1,2
! s state :
                 if(inorb(iat).eq.0)
     &              faktor=sqrt(4.0d0*pi)*faktor
! p states:
                 if(inorb(iat).eq.1)
     &              faktor=sqrt(4.0d0*pi/3.0d0)*faktor
! d states:
                 if(inorb(iat).eq.2) then
                    faktor=sqrt(4.0d0*pi/15.0d0)*faktor
                 end if
! f states:
                 if(inorb(iat).eq.3) then
                    if(imqn.eq.0)
     &                 faktor=sqrt(16.0d0*pi/7.d0)*faktor
                    if(imqn.eq.1)
     &                 faktor=sqrt(32.0d0*pi/21.0d0)*faktor
                    if(imqn.eq.2)
     &                 faktor=sqrt(4.0d0*pi/105.0d0)*faktor
                    if(imqn.eq.3)
     &                 faktor=sqrt(32.0d0*pi/35.0d0)*faktor
                 end if
        end do
! done with faktor
 
! Integration is over z (z-axis points from atom 1 to atom 2)
! and rho (rho is radial distance from z-axis).
 
        do 10 iz=1,nnz
           z(1)=zmin+dfloat(iz-1)*dz
           z(2)=z(1)-d
           xmult=xxmult
           if(mod(iz,2).eq.0)xmult=2*xxmult
           if(iz.eq.1.or.iz.eq.nnz)xmult=0.5d0*xxmult
 
           do 11 irho=1,nnrho
              rho= rhomin+dfloat(irho-1)*drho
              r(1) = sqrt(z(1)**2+rho**2)
              r(2) = sqrt(z(2)**2+rho**2)
 
! Skip if r outside of natm potential range.
! ttest
           if(isorp.eq.0) then
                 if(r(2).ge.rc) go to 11
           else
                 if(r(1).ge.rc) go to 11
           end if
 
              ymult=yymult
              if(mod(irho,2).eq.0)ymult=2*yymult
              if(irho.eq.1.or.irho.eq.nnrho)ymult=0.5d0*yymult
              fact=xmult*ymult*faktor
 
! Calculate psi for each of the two atoms.
! Again, inorb=0,1,2,3 means s,p,d,f-state.
! and    imqn=0,1,2,3 means sigma,pi,delta,phi.
 
! first get radial part:
 
              psi(1)=(alpha**1.5d0)*exp(-alpha*r(1)**2)
              if(isorp.eq.0)psi(2)=vnnaofr(in3,isorp,r(2))
              if(isorp.ne.0)then
               if(r(2).gt.rc3)then
                 psi(2)=1.0d0/r(2)
               else
                 psi(2)=vnnaofr(in3,isorp,r(2))
               end if
              end if
 
              do iat=1,2
! p states:
              if(inorb(iat).eq.1) then
                    if(imqn.eq.0)psi(iat)=psi(iat)*z(1)
                    if(imqn.eq.1)psi(iat)=psi(iat)*rho
                  end if
 
! d states:
              if(inorb(iat).eq.2) then
                    if(imqn.eq.0)psi(iat)=psi(iat)*
     &                                   (2*z(1)**2-rho**2)
                    if(imqn.eq.1)psi(iat)=psi(iat)*rho*z(1)
                    if(imqn.eq.2)psi(iat)=psi(iat)*rho**2
              end if
 
! f states:
              if(inorb(iat).eq.3) then
                    if(imqn.eq.0)psi(iat)=psi(iat)*z(1)*
     &                              (2*z(1)**2-3*rho**2)
                    if(imqn.eq.1)psi(iat)=psi(iat)*rho*
     &                                      (4*z(1)**2-rho**2)
                    if(imqn.eq.2)psi(iat)=psi(iat)*rho**2*z(1)
                    if(imqn.eq.3)psi(iat)=psi(iat)*rho**3
              end if
             end do
              sum=sum+fact*psi(1)*psi(2)*rho
11         continue
10      continue
 
        return
        end
