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

! >>> subroutine gV.f <<<
! Subsection: create.
! Computes the integral of
!  (alpha)**(1.5)*f1*f2*exp(-alpha*r**2)*V(rvec-d*zhat)
! where V is one of Vna, Vss, Vpp, or Vdd, and
! f1 and f2 are each one of the following:
!
!       s  =  1
!       p1 =  x
!       p2 =  y
!       p3 =  z
!       d1 =  xy
!       d2 =  (x**2-y**2)/2
!       d3 =  xz
!       d4 =  yz
!       d5 =  (2*z**2-x**2-y**2)/sqrt(12)
!
! Only ten of the possible combinations of f1 and f2
! give a non-zero and unique result. These are computed
! then tabulated in the following order:
!
!       1. ss
!       2. sp_sigma     (s  and p3)
!       3. pp_sigma     (p3 and p3)
!       4. pp_pi        (p1 and p1 or p2 and p2)
!       5. sd_sigma     (s  and d5)
!       6. pd_sigma     (p3 and d5)
!       7. pd_pi        (p1 and d3 or p2 and d4)
!       8. dd_sigma     (d5 and d5)
!       9. dd_pi        (d3 and d3 or d4 and d4)
!      10. dd_delta     (d1 and d1 or d2 and d2)
!
! The integrals are computed on a 2-d grid of d and sqrt(alpha) values.
!
!-----------------------------------------------------------------------
        subroutine gV(nz,nrho,ndd,in3,rca,isorp,
     1              nspec,signature,what)
!-----------------------------------------------------------------------
! questions and comments should be sent to tomfohr@asu.edu
!-----------------------------------------------------------------------
! rca  - an array whith the Rc's for all atoms;
! nz   - z-direction grid (quadrature.inc);
! nrho - rho-direction grid (quadrature.inc);
! ndd  - dbc grid (grids.inc);
! in3  - type of atom.
!-----------------------------------------------------------------------
 
        implicit none
 
        integer i,id,idx1,idx2,idxrtal,imqn,in3
        integer istop,itpw,itype,iunit,matel,maxmat,ndd,ndime
        integer nrho,nrtal,nspec,nz,isorp
 
        real*8 drtal,hold,rc3,rca,rcin,rcmax,rcmin,rmax,sum
        real*8 alpha,alphamax,d,dr
 
        include '../parameters.inc'
 
! maxmat is the maximum number of matrix elements.
! matel is how many are computed:
! 10 will go up to d orbitals. 20 if f-orbitals are included.
! nrtal=number of (root) alpha values.
        parameter(maxmat=20,matel=10,nrtal=20)
 
        integer*4    mqn(maxmat)
        integer*4    inorb1(maxmat),inorb2(maxmat),inorb(2)
        dimension    hold(maxmat,1107)
        dimension    rca(nspec_max)
        character*70 title,bar,mesg70
        character*1  ang
        character    froot*80, suffix*10
        character*40 fname
        character    what(1:nspec_max)*70,signature*30
 
        bar='===========================================================
     &======='
 
! 4/98: inorb and mqn are used in the integralo subroutine.
! The inorb tell the type of orbital: 0,1,2,3 = s,p,d,f.
! The mqn (m quantum number) tell which type of p d or f:
! 0,1,2,3 = sigma,pi,delta,phi.
! This is explained below after the list of
! the matrix elements (sp_sigma, pp_pi, etc.)
 
        data mqn      /0,0,0,1,0,0,1,0,1,2,
     &                 0,0,1,0,1,2,0,1,2,3/
        data inorb1   /0,0,1,1,0,1,1,2,2,2,
     &                 0,1,1,2,2,2,3,3,3,3/
        data inorb2   /0,1,1,1,2,2,2,2,2,2,
     &                 3,3,3,3,3,3,3,3,3,3/
 
! get rcmax and rcmin
        rcmax=rca(1)
        rcmin=rca(1)
        do i=2,nspec
           if (rca(i).gt.rcmax) rcmax=rca(i)
           if (rca(i).lt.rcmin) rcmin=rca(i)
        end do
 
! set range of alpha values:  zero to alphamax.
! for now !!!!!TEMPORARILY!!!!! set alphamax to (20./rcmin)**2
        alphamax=(20.0D0/rcmin)**2
! rtal = sqrt(alpha).  nrtal is the number of sqrt(alpha)'s
! ("number root alpha's")
! and drtal is the spacing.
!       nrtal=5
        drtal=sqrt(alphamax)/dfloat(nrtal-1)
 
!-----------------------------------------------------------------------
 
        write(*,*)' '
        write(*,*)' *------------------------------*'
        write(*,*)' | Welcome to the gV subroutine |'
        write(*,*)' *------------------------------*'
        write(*,*)' '
        write(*,*)' Who and when: ',signature
        write(*,*)' '
 
! ang is the Angstrom symbol (Å)
        ang = char(197)
 
        write(*,*)' '
        write(*,*)' Integration in cylinder coordinates.'
        write(*,8975)nz,nrho
8975    format('  Number of points for nz,nrho grids =',2i5)
 
! Check dimensions.
 
        ndime=1107
        if(ndd.gt.ndime)stop ' Overlap: Bad dimensions.......'
 
        rc3=rca(in3)
        rmax=rc3+rcmax
 
        dr=rmax/dfloat(ndd-1)
        write(*,7749)ndd
7749    format('  Number of points for distance between two atoms:',i5)
        write(*,6459)dr,ang
6459    format('  Point seperation is ',f9.6,' ',a1,'.')
        istop=9999
 
! Nuclear charges of the two species:
 
!-----------------------------------------------------------------------
! Now begin writing to the output file
!-----------------------------------------------------------------------
 
! Create string with filename:
 
        froot='coutput/gV'
        suffix='dat'
        idx1=in3
        idx2=isorp
        iunit=36
        call iofile2c(froot,suffix,idx1,idx2,iunit,fname)
 
! Set up the header strings:
 
        write(title,1492)signature
1492    format('  gV created by ',a45)
        write(mesg70,1493)nz,nrho,rmax,ndd,nrtal
1493    format(' nz= ',i4,' nrho= ',i4,' R= ',
     &         f8.4, ' (A)',' ndd= ',i4,' nrtal= ',i4)
 
        write(36,119)bar
        write(36,119)title
119     format(a70)
        write(36,55)what(in3)
55      format(a70)
        write(36,119)bar
        write(36,119)mesg70
 
! Standard information ----- nearest neighbor distance.
 
        write(36,*)'rmax and ndd ='
        write(36,*) rmax,ndd
        write(36,*)'rtalmax and nrtal ='
        write(36,*) sqrt(alphamax),nrtal
 
!---------------------------------------------------------------
! do the integrations.
!---------------------------------------------------------------
 
! d = distance between the atoms.
 
        d=-dr
        do 1110 id =1,ndd
        if(id.gt.istop)go to 1110
        d=d+dr
! loop over range of sqrt(alpha) values
        do 1974 idxrtal = 1,nrtal
 
         alpha=((idxrtal-1)*drtal)**2
         itpw=0
 
! tfix let it go 1,10 for d-orbitals and 1,20 for f-orbs.
         do 787 itype = 1,matel
 
              inorb(1)=inorb1(itype)
              inorb(2)=inorb2(itype)
              imqn=mqn(itype)
 
              rcin=rmax
              call intgrlgV(imqn,inorb,nz,nrho,d,rcin,sum,
     &                      in3,itype,alpha,isorp,rc3)
              itpw = itpw + 1
              hold(itpw,id)=sum
 
787      continue
 
         write(36,*)'d=',d
         write(36,*)'sqrt(alpha)=',sqrt(alpha)
         do itype=1,matel
          write(36,*)hold(itype,id)
         end do
 
1974     continue
1110    continue
 
        write(*,*)'  '
        write(*,6573)fname
6573    format('  Finished. Wrote output to ',a40)
        write(*,*)'  '
        write(*,119)bar
        write(*,*)' '
        close(unit=36)
        return
        end
