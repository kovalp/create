c copyright info:
c
c                             @Copyright 2005
c                           Fireball Committee
c Brigham Young University - James P. Lewis, Chair
c Arizona State University - Otto F. Sankey
c Lawrence Livermore National Laboratory - Kurt Glaesemann
c Universidad de Madrid - Jose Ortega
c Universidad de Madrid - Pavel Jelinek
c Brigham Young University - Hao Wang

c Other contributors, past and present:
c Auburn University - Jian Jun Dong
c Arizona State University - Gary B. Adams
c Arizona State University - Kevin Schmidt
c Arizona State University - John Tomfohr
c Motorola, Physical Sciences Research Labs - Alex Demkov
c Motorola, Physical Sciences Research Labs - Jun Wang
c Ohio University - Dave Drabold
c University of Regensburg - Juergen Fritsch

c
c RESTRICTED RIGHTS LEGEND
c Use, duplication, or disclosure of this software and its documentation
c by the Government is subject to restrictions as set forth in subdivision
c { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
c clause at 52.227-7013.
c
c Program Description
c ===========================================================================
c       gauss_density.f
c
c       This program calculate the electron density using the wavefunctions
c       of a given atom and return to gauss_create.f fitting to gaussians
c       MHL
c       It was modified to get electron density for each shell in additon to
c       the total electron density.
c ===========================================================================
c Code written by:
c	Hao Wang
c	Department of Physics & Astronomy
c	Arizona State University
c	Tempe, AZ 85287
c	Hao.Wang@asu.edu
c
! Program Description
c ===========================================================================
c Old statement
c        subroutine gauss_density(specfile,rc,npnts,rnorm,dens)
c MHL
        subroutine gauss_density (ish, specfile, rc, npnts, rnorm,      
     1                            dens)
        implicit none
c Input
c ===========================================================================
        integer nsh_max
        integer wfmax_points
c MHL
        integer ish
        parameter (nsh_max=8)
        parameter (wfmax_points=10000)

        character*10 specfile
c Local
c ===========================================================================
        integer nzx
        integer nssh 
        integer lssh (nsh_max)
        integer lqn 
        integer issh
        integer ip
        integer npnts
        integer npnts1
        integer npnts2
        integer npnts3
        integer idenp
        integer nwf

        real*8 dens(wfmax_points)
        real*8 rcutoff(nsh_max)

        character*25 wavefxn

        real*8 rcutoff_max_diff
        real*8 rcutoff1
        real*8 rcutoff2
        real*8 rcutoff3
        real*8 rc
        real*8 rnorm
        real*8 rcin
        real*8 xnocc(nsh_max)
        real*8 xnoccin

        real*8 a
        real*8 b
        real*8 r
        real*8 psir
        real*8 x
        real*8 dx

        real*8 sum,pi
        real*8 dens_r_int
        real*8 xnocc_sum
          
        real*8 Gpsiofr
c Common block
c ===========================================================================
        integer npoints (nsh_max)

        real*8 drr (nsh_max)
        real*8 rr (wfmax_points, nsh_max)
        real*8 rrc (nsh_max)
        real*8 psi (wfmax_points, nsh_max)

c Proceedure
c ===========================================================================   
        pi=4.0*atan(1.0)
        rnorm = 0.d0
c
c Opening the specfile
        open(unit=29,file=specfile,status='unknown')
        read(29,*)
        read(29,*) nzx
        read(29,*)
        read(29,*)
        read(29,*)
        read(29,*) nssh

c Loop over the number of shells:
c Read the L-quantum number (0, 1, 2, and 3 => s, p, d, and f),
c the occupation number, cutoff radius (in bohr), wavefunction
c and neutral atom potential for each shell.
c
        do issh = 1, nssh
        read (29,*)  lssh(issh)
        read (29,*)  xnocc(issh)
        read (29,*)  rcutoff(issh)
        read (29,25) wavefxn
        read(29,*)
25      format(a25)
 
        rcin= rcutoff(issh)
        xnoccin=xnocc(issh)
        lqn=lssh(issh)
c Note: Greadpsi is the Gauss (G) version of readpsi.
        call Greadpsi(issh,lqn,nzx,rcin,xnoccin,wavefxn,
     1               drr,rr,rrc,npoints,psi)
c
c Get the biggest rccutoff and npoints among the involved shell
        if(issh.eq.1)then
        npnts=npoints(issh)
        rc=rrc(issh)

        npnts1=npoints(issh)
        rcutoff1=rrc(issh)
        endif
c
        if(issh.eq.2)then
        npnts2=npoints(issh)
        rcutoff2=rrc(issh)
        endif
c
        if(issh.eq.3)then
        npnts3=npoints(issh)
        rcutoff3=rrc(issh)
        endif
c
        if(issh.gt.1)then
c
          if(issh.eq.2)then
          rc = max(rcutoff1, rcutoff2)
          npnts = max(npnts1, npnts2)
          rcutoff_max_diff = abs(rcutoff1 - rcutoff2)
          endif
c
          if(issh.eq.3)then
          rc = max(rcutoff1, rcutoff2, rcutoff3)
          npnts = max(npnts1, npnts2, npnts3)
          rcutoff_max_diff = max(abs(rcutoff1 - rcutoff2),
     1                           abs(rcutoff2 - rcutoff3),
     2                           abs(rcutoff3 - rcutoff1))
          endif
c
c If the cutoffs for different shells are too drastically different,
c then the grid size should be increased.
c Otherwise, the number of non-zero points may be too few.
c
           if (rcutoff_max_diff .gt. 1.0d0) then
           write (*,*) '  '
           write (*,*) ' ===============   WARNING   ================='
           write (*,*)
           write (*,*) ' rcutoff_max_diff :',rcutoff_max_diff
           write (*,*)
           write (*,*) ' You have at least two shells which have '
           write (*,*) ' rcutoff''s which differ by more than 2.0 '
           write (*,*) ' Angstroms. It is highly advisable that you '
           write (*,*) ' increase the number of mesh points, so as to '
           write (*,*) ' avoid a case where you may end up with many '
           write (*,*) ' zeros, and too few non-zero elements in your '
           write (*,*) ' grid. '
           write (*,*)
           write (*,*) ' ========= PLEASE CHECK INPUT FILES =========='
           write (*,*) '  '
           write (*,*) '  '
           stop
           end if
c
        endif
c
c end of issh loop
        end do
        close(29)

c Calculate electron density --- dens(i)
        A=0.d0
        B=rc
        dx=(B-A)/(npnts-1.d0)
        do ip=1,npnts
        dens(ip)=0.d0
        r=A+(ip-1)*dx
c total electron density
        if (ish.eq.-1) then
        do issh=1,nssh
          psir=Gpsiofr(issh,r,drr,rr,rrc,npoints,psi)
c Note about density. Since we integrate over all rvector (r, theta, phi),
c we need to include the Y00 in psi for the density. Y00**2=1/(4.*pi).
c This applies to p or d states as well. 1/3 sum(m) |Y1m|^2=1/(4*pi) etc.
c Correction. 1/4pi added Oct. 1, 2004.
          dens(ip)=dens(ip)+(psir**2)*xnocc(issh)/(4.*pi)
        end do
        else
        issh=ish+1
c electron density for each shell
c issh=0 for s-shell, 1 for p-shell, etc.
        psir=Gpsiofr(issh,r,drr,rr,rrc,npoints,psi)
c Note about density. Since we integrate over all rvector (r, theta, phi),
c we need to include the Y00 in psi for the density. Y00**2=1/(4.*pi).
c This applies to p or d states as well. 1/3 sum(m) |Y1m|^2=1/(4*pi) etc.
c edens for each shell is not multiplied by occupation number
c Correction. 1/4pi added Oct. 1, 2004.
        dens(ip)=(psir**2)/(4.*pi)
        end if
        end do
c Calculate normalization constant of calculated electron density
c integrate from A to B  (B > A) --- rnorm
        A=0.
        B=rc
        dx=(B-A)/(npnts-1.)
        sum=((dens(1)*A)**2+(dens(npnts)*B)**2)/2.
        do ip=1,npnts-2
        x=A+ip*dx
        sum=sum+(dens(ip+1)*x)**2
        enddo
        rnorm=sum*dx

c
c Check the self-consistency of the calculated electron density
c integrate from A to B  (B > A) --- dens_r_int
c MHL
        if (ish.eq.-1) then
        A=0.
        B=rc
        dx=(B-A)/(npnts-1.)
        sum=(dens(1)*A**2+dens(npnts)*B**2)/2.
        do ip=1,npnts-2
        x=A+ip*dx
c Correction. 1/4pi added Oct. 1, 2004.
        sum=sum+dens(ip+1)*(x**2)*(4.0*pi)
        end do
        dens_r_int=sum*dx

c summation of  occupations for all shell --- xnocc_sum
        sum=0.d0
        do issh=1,nssh
        sum=sum+xnocc(issh)
        end do 
        xnocc_sum=sum

        if(dabs(dens_r_int - xnocc_sum) .gt. 1.e-4)then
          write(*,*)' xnocc_sum ',xnocc_sum,' dens_r_int ',dens_r_int
          write(*,*)'  '
          write(*,*)' Error detected between the normolization '
          write(*,*)' and the summation of real occupations '
          stop 'error in gauss_density.f'
        end if
c ish        
        end if
        return
        end

! ===========================================================================
! Gpsiofr.f
! Program Description
! ===========================================================================
!       This function returns the values Gpsiofr(r) for the corresponding
! shell of the atomtype.  The radial functions are normalized as:
 
!  int ( Gpsiofr**2  r**2  dr ) = 1.0
 
! The wavefunctions for each atom must be read in by calling Greadpsi
! separately for each atom type.
! Note: Greadpsi is the Gauss (G) version of readpsi.
 
! The value of r input to this function must be in angstrom units.
!
! ===========================================================================
! Program Declaration
! ===========================================================================
        function Gpsiofr (issh,r,drr,rr,rrc,npoints,psi)
        implicit none

        real*8 Gpsiofr
 
! Argument Declaration and Description
! ===========================================================================
        integer nsh_max
        integer wfmax_points
        parameter (nsh_max=8)
        parameter (wfmax_points=10000)

        integer issh
        real*8 r
 
! Local Parameters and Data Declaration
! ===========================================================================
        integer norder
        parameter (norder = 5)

        integer norder2
        parameter (norder2 = norder/2)

        logical splines
        parameter (splines = .true.)
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ileft, imid, iright
        integer iprod
        integer isum
        integer max_point
 
        real*8 prod
        real*8 rr1
        real*8 rr2
        real*8 psiofr_tmp

!       Cubic spline specific variables
        real*8 L(0:norder),mu(0:norder),Z(0:norder),alpha(0:norder)
        real*8 a(0:norder),b(0:norder),c(0:norder),d(0:norder)
        real*8 h
        integer iam
        integer i,j
        real*8 xmin
        real*8 xxp

! Common block 
! ===========================================================================  

        integer npoints (nsh_max)

        real*8 drr (nsh_max)
        real*8 rr (wfmax_points, nsh_max)
        real*8 rrc (nsh_max)
        real*8 psi (wfmax_points, nsh_max)

! Procedure
! ===========================================================================
! Special cases
       if (r .gt. rrc(issh)) then
         Gpsiofr = 0.0d0
         return
       else if (r .lt. 0.0d0) then
         Gpsiofr = psi(1,issh)
         return
       end if
 
! note : the points are equally spaced
       imid = int(r/drr(issh)) + 1
 
! Find starting point for the interpolation
       ileft = imid - norder2
 
       max_point = npoints(issh)
       if (ileft .lt. 1) then
        ileft = 1
       else if (ileft + norder .gt. max_point) then
        ileft = max_point - norder
       end if
 
! Find ending point for the interpolation
       iright = ileft + norder
 
       if (.not. splines) then

! Now interpolate with polynomials
        psiofr_tmp = 0.0d0
        do isum = ileft, iright
         prod = 1.0d0
         rr1 = rr(isum,issh)
         do iprod = ileft, isum - 1
          rr2 = rr(iprod,issh)
          prod = prod*(r - rr2)/(rr1 - rr2)
         end do
         do iprod = isum + 1, iright
          rr2 = rr(iprod,issh)
          prod = prod*(r - rr2)/(rr1 - rr2)
         end do
        psiofr_tmp = psiofr_tmp + psi(isum,issh)*prod
        end do
        Gpsiofr = psiofr_tmp

       else
! Now interpolate with "natural" splines with f''(x)=0 at end points

        do i=0,norder
          a(i)=psi(i+ileft,issh)
        end do
! note : We use one value of h, so points must be equally spaced
        h=drr(issh)
        do i=1,norder-1
          alpha(i)=3.0d0*(a(i+1)-2*a(i)+a(i-1))/h
        end do

        L(0)=1
        mu(0)=0
        Z(0)=0
        c(0)=0
        do i=1,norder-1
          L(i)=(4.0d0-mu(i-1))*h
          mu(i)=h/L(i)
          Z(i)=(alpha(i)-h*Z(i-1))/L(i)
        end do
        L(norder)=1
        mu(norder)=0
        Z(norder)=0
        c(norder)=0

!       Do not go off of the end
        if(imid .eq. ileft)imid=iright-1
!       What curve section do we use?
        iam=imid-ileft

!       Don't need 0 to iam-1
        do j=norder-1,iam,-1
          c(j)=z(j)-mu(j)*c(j+1)
          b(j)=(a(j+1)-a(j))/h-h*(c(j+1)+2.0d0*c(j))/3.0d0
          d(j)=(c(j+1)-c(j))/(3.0d0*h)
        end do

        xmin = 0.0d0
        xxp=(r-(xmin+(imid-1)*h))

        Gpsiofr=a(iam)+b(iam)*xxp+c(iam)*xxp**2+d(iam)*xxp**3

       end if
 
       return
       end
! ===========================================================================
 
! Greadpsi.f
! Note: Greadpsi is the Gauss (G) version of readpsi.
! Program Description
! ===========================================================================
!        This subroutine Greadpsi reads the values for the radial wavefunctions
! of the different (s,p,d,f) orbitals off the data file atomicsymbol.wfl
! which comes from rcatms.f. The normalization will also be checked here.
! Note: Psi is already in Angstrom units on the data file.
 
! rc must be input to this subroutine in abohr units.
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine Greadpsi(issh,lqn,nzx,rcutoff,xnoccin,filein,
     1                     drr,rr,rrc,npoints,psi) 

        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer issh
        integer lqn
        integer nzx
 
        real*8 rcutoff
        real*8 xnoccin

        character*25 filein
 
        integer nsh_max
        integer wfmax_points
        parameter(nsh_max=8)
        parameter(wfmax_points=10000)

! Local Variable Declaration and Description
! ===========================================================================
        integer inum
        integer ipoint
        integer iremainder
        integer lqnwf
        integer mesh
        integer nzxwf
 
        real*8 r
        real*8 rc
        real*8 rc_max
        real*8 rcutoffwf
        real*8 sum
        real*8 xnoccwf
        real*8 abohr
 
        real*8 psitemp(wfmax_points)
 
        character*25 fileinwf
        real*8 psirr 

        logical iammaster

! Common block 
! ===========================================================================

        integer npoints (nsh_max)

        real*8 drr (nsh_max)
        real*8 rr (wfmax_points, nsh_max)
        real*8 rrc (nsh_max)
        real*8 psi (wfmax_points, nsh_max)
! ===========================================================================

! Procedure
        iammaster = .false.
! some constant:
        abohr=0.529177

! Open the input file
        open (unit = 30, file = filein, status = 'old')

        if (iammaster) then
         write (*,*) ' '
         write (*,*) '*-----------------------------------------------*'
         write (*,*) '|               Welcome to GREADPSI             |'
         write (*,*) '| Reading the radial wavefunction of your atom  |'
         write (*,*) '*-----------------------------------------------*'
         write (*,*) ' '
        end if ! end master
 
        read (30,90) fileinwf
        read (30,*) nzxwf
        read (30,*) mesh
        read (30,*) rcutoffwf, rc_max, xnoccwf
        read (30,*) lqnwf
        if (iammaster) write (*,91) fileinwf
 
! Perform some checks
        if (nzxwf .ne. nzx) then
         write (*,*) ' nzxwf = ', nzxwf, ' nzx = ', nzx
         write (*,*) ' The Z number in the wavefunction file, '
         write (*,*) ' for this shell, does not match the Z number '
         write (*,*) ' that you put into the create.input file. '
         write (*,*) ' Double check everything and rerun creator.'
         stop 'error in Greadpsi'
! Note: Greadpsi is the Gauss (G) version of readpsi.
        end if
 
        if (rcutoffwf .le. (rcutoff - 1.0d-2) .or.
     1      rcutoffwf .ge. (rcutoff + 1.0d-2)) then
         write (*,*) ' rcutoffwf = ', rcutoffwf, ' rcutoff = ', rcutoff
         write (*,*) ' The cutoff radius in the wavefunction file, for '
         write (*,*) ' this shell, does not match the cutoff radius '
         write (*,*) ' that you put into your create.input file. '
         write (*,*) ' Double check everything and rerun creator.'
         stop 'error in Greadpsi'
        end if
 
        if (xnoccwf .ne. xnoccin) then
         write (*,*) ' xnoccwf = ', xnoccwf, ' xnoccin = ', xnoccin
         write (*,*) ' The occupation number in the wavefunction file, '
         write (*,*) ' for this shell, does not match the occupation '
         write (*,*) ' number that you put into your create.input file.'
         write (*,*) ' Double check everything and rerun creator.'
         stop 'error in Greadpsi'
        end if
 
        if (lqnwf .ne. lqn) then
         write (*,*) ' lqnwf = ', lqnwf, ' lqn = ', lqn
         write (*,*) ' The l quantum number in the wavefunction file, '
         write (*,*) ' for this shell, does not match the l quantum '
         write (*,*) ' number that you put into your create.input file.'
         write (*,*) ' Double check everything and rerun creator.'
         stop 'error in Greadpsi'
        end if
 
        if(mesh .gt. wfmax_points) then
         write (*,*) ' Error error ***** in Greadpsi. '
         write (*,*) ' Dimension of wavefunction = ', wfmax_points
         write (*,*) ' We are asking for mesh = ', mesh
         write (*,*) ' Redimension wfmax_points. '
         stop 'error in Greadpsi'
        end if
 
! Set some things up
        rc = rcutoffwf*abohr
!        rc = rcutoffwf
 
        if (iammaster) then
         write (*,*) '  '
         write (*,200) rcutoffwf, rc
        end if ! end master
        npoints(issh) = mesh
        rrc(issh) = rc
        drr(issh) = rc/dfloat(mesh - 1)

! Read in the points
        inum = idint(dfloat(mesh)/4)
        iremainder = mesh - (inum*4)
        do ipoint = 1, mesh - iremainder, 4
         read (30,100) psitemp(ipoint), psitemp(ipoint+1),
     1                 psitemp(ipoint+2), psitemp(ipoint+3)
        end do
 
        if (iremainder .eq. 1) then
         read (30,100) psitemp(mesh)
        else if (iremainder .eq. 2) then
         read (30,100) psitemp(mesh-1), psitemp(mesh)
         else if (iremainder .eq. 3) then
         read (30,100) psitemp(mesh-2), psitemp(mesh-1), psitemp(mesh)
        end if
        close (30)
! Write psitemp to the wavefunction psi

        do ipoint = 1, mesh
         psi(ipoint,issh) = psitemp(ipoint)
        end do

! Check normalization
        if (iammaster) then
         write (*,*) '  '
         write (*,*) ' Checking normalization [NORM(l) should be 1]'
         write (*,*) '  '
        end if ! end master
        r = - drr(issh)
        do ipoint = 1, mesh
         r = r + drr(issh)
         rr(ipoint,issh) = r
        end do
        sum = 0.0d0
        do ipoint = 1, mesh
         if (ipoint .ne. 1 .or. ipoint .ne. mesh) then
          sum = sum + drr(issh)*rr(ipoint,issh)**2
     1                         *psi(ipoint,issh)**2
         else
          sum = sum + 0.5d0*drr(issh)*rr(ipoint,issh)**2
     1                               *psi(ipoint,issh)**2
         end if
        end do
        if (iammaster) then
         write (*,300) issh, sum
         write (*,*) '  '
         write (*,*) ' *--------------- END GREADPSI -----------------*'
         write (*,*) '  '
        end if ! end master
 
! Code to just make ftnchek happy, because rc_max is not ever used
        if (rc_max .eq. 0.0d0) return
 
! Format Statements
! ===========================================================================
90      format (a25)
91      format (2x,a25)
100     format (4d18.10)
200     format (' Rc read in =', f8.4, ' abohr = ', f8.4, ' Angstroms')
300     format (' NORM (shell = ', i1, ') = ', f16.12)
 
        return
        end
