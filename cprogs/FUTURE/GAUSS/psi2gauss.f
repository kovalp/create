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

        subroutine psi2gauss(ispec,nssh,lssh,signature,what)
 
! p2g = psi to gauss.
! The program approximates the radial part of an atomic
! wave function as a sum of gaussians.  For s wave funcions
! the gaussians are exp(-alpha*r**2), for p they are
! r*exp(-alpha*r**2), for d they are r**2*exp(-alpha*r**2),
! etc...
!===================================================================
! lmax is the maximum l (quantum number) value and is 1
! because we're only looking at s and p right now.
!
! ng is the number of gaussians (4 should be good).
!
! nmax is the maximum number of points for the wave function.
!
! nmag is a number that determines the range of alpha values.
 
        implicit none
 
        real*8 al,alidx,alpha,c,coef,dr,error,g,pi,psi2,r
        real*8 ra,rb,rca,rmid,rmin,rnorm,rnorm2,rpsi,s
        real*8 sinv,sum2
 
        integer i,index,ispec,issh,iunit,j,k,lq,lssh,m
 
        integer magic,nbig,ng,nmag,nmax,nng,npnts,nssh
 
!       integer lmax
!       parameter (lmax=3)
 
        include '../parameters.inc'
        include '../wavefunctions.inc'
 
        parameter(ng=20,pi=3.1415926536d0,nmax=5000,nmag=19)
 
! temporary. nng=number of gaussians in the expansions.
        parameter(nng=2)
 
        dimension nssh(nspec_max),lssh(nspec_max,nsh_max)
        dimension ra(ng), S(ng,ng), g(ng,nmax)
        dimension Sinv(ng,ng), rb(nmag,ng), psi2(nmag,nmax)
        dimension rnorm2(nmag), error(nmag)
        dimension al(ng), M(ng,ng), C(ng,ng)
        character froot*80, suffix*10
        character*70 bar
        character    what(1:nspec_max)*70,signature*30
 
        write(*,*)' Fitting the wavefunctions of species: ',ispec
        write(*,*)' Using ',nng,' gaussians.'
 
!===================================================================
 
! start writing to file
 
        bar='====================================================='
 
 
        froot='coutput/p2g'
        suffix='dat'
        iunit=36
        index=ispec
        call iofile(froot,suffix,index,iunit)
 
        write(36,*)bar
        write(36,*)' created by: ',signature
        write(36,*)what(ispec)
        write(36,*)bar
        write(36,*)' number of shells='
        write(36,*)nssh(ispec)
        write(36,*)' number of gaussians='
        write(36,*)nng
        write(36,*)
 
! BIGGEST LOOP.  loop over shells.
        do issh=1,nssh(ispec)
 
        lq=lssh(ispec,issh)
        rca=rrc(issh,ispec)
        npnts=npoints(issh,ispec)
        dr=(rca/(npnts-1))
 
! BIGGER LOOP.  First time through it samples a course range
! of groups of alphas and finds the best group.  Second time
! it looks for a better group of alphas in the vicinity of
! the best one it found the first time around.
! mid is the middle of the range of alidx values sampled.
! After the first run through the BIGGER loop, rmid is set
! equal to the alidx of the best group...
! Initially:
        rmid=2.0d0
! START BIGGER LOOP
        do nbig=1,2
! START BIG LOOP
        do magic=1,nmag
! The 10.**nbig is just to modify the coarseness of the
! range of alidx values.
        alidx=rmid+dfloat(magic-10)/dfloat(10**nbig)
! >>> Get the alpha values <<<
        do i=1,nng
         al(i)=((alidx*i)/rca)**2
        end do
! Get the gaussians.
        if (lq.eq.0) then
         do i=1,nng
          do k=1,nmax
           r=(k-1)*dr
           g(i,k)=exp(-al(i)*r**2)
          end do
         end do
        else if (lq.eq.1) then
         do i=1,nng
          do k=1,nmax
           r=(k-1)*dr
           g(i,k)=r*exp(-al(i)*r**2)
          end do
         end do
        else if (lq.eq.2) then
         do i=1,nng
          do k=1,nmax
           r=(k-1)*dr
           g(i,k)=r**2*exp(-al(i)*r**2)
          end do
         end do
        else if (lq.eq.3) then
         do i=1,nng
          do k=1,nmax
           r=(k-1)*dr
           g(i,k)=r**3*exp(-al(i)*r**2)
          end do
         end do
        end if
! Initialize ra and rb.
! ra(i) will later be given the value of the integral
! of psi times the ith gaussian over all space.
! rb(*,i) will be the coefficient of the ith gaussian
! in the expansion of psi.
        do i=1,nng
         ra(i)=0.0d0
         rb(magic,i)=0.0d0
        end do
! Initialize Sij.
! Sij=the scalar product of the ith with the jth gaussian.
! It can be worked out analytically but for now I just
! do it the dumb way...
        do i=1,nng
         do j=i,nng
          S(i,j)=0.0d0
          do k=1,nmax-1
           r=(k-1)*dr
           S(i,j)=S(i,j)+dr*2*pi*(r**2*g(i,k)*g(j,k)
     1            +(r+dr)**2*g(i,k+1)*g(j,k+1))
           end do
           S(j,i)=S(i,j)
         end do
        end do
! Get the inverse of S.
        do i=1,nng
         do j=1,nng
          Sinv(i,j)=S(i,j)
         end do
        end do
        call invert(Sinv,nng,ng,M,C)
 
! integrate psi*g over all space.
! Use trapezoid rule.
 
        do i=1,nng
         do j=1,npnts-1
          r=(j-1)*dr
          ra(i)=ra(i)+dr*2*pi*(r**2*g(i,j)*psi(j,issh,ispec)
     1          +(r+dr)**2*g(i,j+1)*psi(j+1,issh,ispec))
         end do
        end do
 
! Write out the coefficients of the gaussians.
        do i=1,nng
         do j=1,nng
          rb(magic,i)=rb(magic,i)+Sinv(i,j)*ra(j)
         end do
        end do
 
! Write psi2=the sum of gaussians that approximates psi.
        do i=1,nmax
         r=(i-1)*dr
         psi2(magic,i)=0.0d0
         do j=1,nng
          psi2(magic,i)=psi2(magic,i)+rb(magic,j)*g(j,i)
         end do
        end do
 
! normalize psi2 and psi such that the
! integral from zero to infinity of (r*psi)**2
! is equal to one.  (Trapezoid rule again).
! SHOULD DO THIS ANALYTICALLY!!!!!!
        sum2=0.0d0
        do i=1,nmax-1
         r=(i-1)*dr
         sum2=sum2+dr*0.5d0*((r*psi2(magic,i))**2+
     1       ((r+dr)*psi2(magic,i+1))**2)
        end do
 
       rnorm2(magic)=1.0d0/sqrt(sum2)
! Estimate the error.
       error(magic)=0.0d0
       do i=1,nmax
        r=(i-1)*dr
! If r is greater than rc then psi is zero...
        if (i.gt.npnts) then
           rpsi=0.0d0
        else
           rpsi=r*psi(i,issh,ispec)
        end if
        error(magic)=error(magic)+
     1     (rpsi*abs(rpsi)-
     2      r**2*psi2(magic,i)*abs(psi2(magic,i)))**2
       end do
       error(magic)=sqrt(error(magic))
! END OF BIG LOOP:
       end do
 
! which has the smallest error?
       rmin=error(1)
       magic=1
       do i=2,nmag
        if (error(i).le.rmin) then
          rmin=error(i)
          magic=i
        end if
       end do
 
! remember: rmid is the middle alidx value sampled.
! At the end it becomes the best alidx value.
       rmid=rmid+dfloat(magic-10)/dfloat(10**nbig)
       rnorm=rnorm2(magic)
! END OF BIGGER LOOP
       end do
 
! >>>Write alpha values and the coefficients using the
! proper normalization (the expansion is a sum of b*exp(al*r**2)
! or b*r*exp(al*r**2) for a p wave function.
 
       write(36,*)' '
       write(36,*)' l='
       write(36,*)lq
       write(36,*)' alpha, coef ='
       do i=1,nng
        alpha=(i*rmid/rca)**2
        coef=rb(magic,i)*rnorm
        write(36,*)alpha, coef
       end do
 
! END OF BIGGEST LOOP
       end do
 
       close(unit=36)
 
      end
