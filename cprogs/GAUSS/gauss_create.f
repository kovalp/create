c ==================== gauss_create.f ================================
c gauss_create.f:          Version of March 26, 2002.
c This version was ported over to the official Fireballs-2002
c program package on May 20, 2002 by Otto F. Sankey and
c James P. Lewis.
c We changed the program so that it computes both NA and XC. It does not
c matter if you want just one or the other. You always get both.
c ====================================================================
c Written by: 
c John Tomfohr concerning neutral atom potential.
c Hao Wang concerning electron density fitting
c Mail comments and questions to john.tomfohr@asu.edu hao.wang@asu.edu
c =====================================================================
c This program fits wavefunctions, NA potentials, and electron 
c densities to gaussians.
c It does this for the 3-center interactions. This removes the need
c for long tables of integrals in the standard method
c
c It will use the smallest number of gaussians--less than or equal 
c to "maxals"--that will give an error below the tolerance levels 
c maxals = the maximum number of gaussians.

c psi_tolerance (error allowed in wavefunction). 
c rna_tolerance (error allowed in neutral atom potential).
c den_tolerance (error allowed in electron density).

c Increasing maxals and decreasing the tolerances will give you 
c better fits BUT (!!!) large maxals and low tolerances will lead to
c a lot of closely spaced gaussians cancelling each other like this:
c    1000000000.*exp(-1.01*r**2) - 999999999.*exp(-1.02*r**2)
c This can lead to VERY bad roundoff errors in fireball.
c
c The best way to avoid this is to have maxals.le.6
c 
c Good values for tolerances are: 
c psi_tolerance = 1.E-3.
c rna_tolerance = 1.E-3. 
c den_tolerance = 1.E-3.
c
c Note: you may be able to get away with more gaussians/better fits
c but you better know what you're doing...

c these things are now determined by choosing the accuracy level...
c =================================================================
        subroutine gausscreate(naccuracy3Cgauss)
c =================================================================
        implicit double precision (a-h,o-z)
        integer naccuracy3Cgauss
        character*60 char1,char2,charsum1,charsum2,comparefile

        parameter (nshells_max = 8)
        parameter (nspec_max = 8)

        character*60 comparefilename(nspec_max, -1:4*nshells_max)
        integer iatomic(nspec_max)
        logical skip

        parameter (maxmaxals=10)

        parameter (rmax_coef=2000.)
        parameter (amount_beyond_rc_max=0.005)
        parameter (ncpnts=100)
        parameter (ndpnts=10)
c These parameters basically determine how refined the search for alphas is.
c Higher = finer search grids basically.
c ncpnts=100 and ndpnts=10 is good.

        parameter (maxpnts=10000)

        character*10 specfile
        character*25 nafile
        character*25 wavefile

        logical garbage
        logical warnings

        character ccontinue

        dimension numshl(12)

        dimension psi(maxpnts),sumfit(100)
        dimension f(maxmaxals),coef(maxmaxals)
        dimension G(maxmaxals,maxmaxals)
c        dimension Gtemp(maxmaxals,maxmaxals)
        dimension alpha(maxmaxals)
c workspace matrices for inversion subroutine
        dimension M(maxmaxals,maxmaxals),C(maxmaxals,maxmaxals)

        real*8 xnoccin

! the lagerest r in NA potential for given shell
        real R_na
c Proceedure
c =================================================================

c some constants:
        bohr=0.529177
        pi=4.d0*atan(1.d0)

c initialize warnings
        warnings=.false.

c WANG gauss, April 17, 2002.
c In case of density fit, we may get unphysical fitting results of  
c electron density with negative values. 
c To avoid this happen and improve the fitting, we introduce the a
c parameter called weight factor to carry out weighting fit.
        weight = 0.d0

        write(*,*)' ================================================== '
        write(*,*)'              WELCOME TO  GAUSSIAN CREATOR          '
        write(*,*)'                 (May 8, 2005 Edition )             '
        write(*,*)
        write(*,*)'          What level of accuracy would you like?    '
        write(*,*)' (Note: fireball run time increases dramatically as '
        write(*,*)' the accuracy is increased. So try to get away with '
        write(*,*)' the lower accuracy if you can.)        '
        write(*,*)
        write(*,*)'       1 = Low.'
        write(*,*)'       2 = Medium. (Default)'
        write(*,*)'       3 = High. (To be safe!)'
        write(*,*)'       4 = Very High (overkill).'
        write(*,*)' ================================================== '
        write(*,*)

        naccuracy=naccuracy3Cgauss
!        naccuracy = 4 
        if (naccuracy.eq.1) then
          psi_tolerance=0.001
          rna_tolerance=0.001
          den_tolerance=0.001
          maxals=3
        else if (naccuracy.eq.2) then
          psi_tolerance=0.0001
          rna_tolerance=0.001
          den_tolerance=0.001
          maxals=4
        else if (naccuracy.eq.3) then
          psi_tolerance=0.0001
          rna_tolerance=0.001
          den_tolerance=0.001
          maxals=5
        else if (naccuracy.eq.4) then
          psi_tolerance=0.00001
          rna_tolerance=0.0001
          den_tolerance=0.0001
          maxals=8
        end if

        if(naccuracy.gt.4.or.naccuracy.le.0)then
          stop 'naccuracy too small or too large. Fix gauss-create.'
        end if
c
        write(*,*)
        write(*,*)' Accuracy level = ',naccuracy
        write(*,*)
        write(*,*)' Under this accuracy level'
        write(*,*)' Maximum number of gaussians allowed is: ',maxals
        write(*,*)
        write(*,71)psi_tolerance
        write(*,72)rna_tolerance
        write(*,73)den_tolerance
c
71      format('  Wave-function          tolerance:',f12.6)
72      format('  Neutral atom potential tolerance:',f12.6)
73      format('  Density                tolerance:',f12.6)

        if (maxals.gt.maxmaxals) then
        write(*,*)' maxals.gt.maxmaxals'
        write(*,*)' Increase the parameter maxmaxals and recompile.'
        stop
        end if
c ========================================================================
c Open input file create.input.
        open(unit=1,file='create.input',status='unknown')
c
c -------------------------------------------------------------------
c     ish            Job
c -------------------------------------------------------------------
c     -1             Ntot                     total electron density
c      0             NA                       neutral atom potential
c      1             Psi(1)     and    n(1)   for   S_state        
c      2             Psi(2)     and    n(2)   for   P_state
c      3             Psi(3)     and    n(3)   for   D_state
c     ...            etc.
c -------------------------------------------------------------------
c output file = gauss.dat
        open(unit=9,file='./coutput/gauss.dat',status='unknown')

        read(1,*)nspecs
        write(*,*)' '
        write(*,*)' Number of species is:',nspecs
        write(9,*)' ===================================='
        write(9,*)' Gaussian fits: Accuracy level=',naccuracy
        write(9,*)
        write(9,*)' Number of species'
        write(9,*)nspecs
c	
c ************* ish0 is ALWAYS -1 *************
                    ish0 = -1
c ish = -1 means that you choose to fit everything, i.e., 
c electron density, NA potentials, and wavefunctions to gaussians.
c
c If you choose ish0 = 0, only wavefunction and NA will be fitted.
c
c =========================================================================
c Write ish0 to gauss.dat
        write(9,*)' Fitting selection: ish0 = 0 OR -1 '
        write(9,*) ish0
c =========================================================================
c Loop over all species.
        do 1000 ispec=1,nspecs
        read(1,526)specfile
526     format(a10)
        write(*,*)'  '
        write(*,*)' =================================================='

        if(ish0.eq.-1)then
c WANG density fit
c Electron density fitting
c This call calculate electron density (output called psi) using
c the wavefunctions of atom.
c
c We put this part outside of loop ish since in gauss_density() 
c we need to open and read data of wavefunction files of atom.
c All input files used in this call are closed, so we can use them 
c without problem in other parts of program.
c
        ish=ish0
        call gauss_density(ish, specfile, rc, npnts, rnorm, psi)
c
c Note even though the line above syas "psi" it is really total electron density. 
c Be careful about this. Stuff like this also happens lateer on. 
c This "psi" array holds different things.
        endif

        write(*,*)' Opening specfile:',specfile
        open(unit=2,file=specfile,status='unknown')
c
        read(2,*)
        read(2,*)n_atomic_number
        write(9,*)' ===================================='
        write(9,*)' Species ',ispec,' atomic number'
        write(9,*) n_atomic_number
        iatomic(ispec)=n_atomic_number

        write(*,*)' Begin Species:',ispec,' Atomic Number:',
     1              n_atomic_number

        read(2,*)
        read(2,*)
        read(2,527)wavefile
        read(2,*)nshells
        write(9,*)' Number of shells'
        write(9,*) nshells
c
        numshl(ispec)=nshells
c
c Old statement
c        do ish=ish0,nshells
c
c MHL to include spherical approximation (psi -> sqrt(psi**2))
        do 1500 ish = ish0, nshells*2
c
        do ial = 1, maxals
           alpha (ial) = 0.d0
            coef (ial) = 0.d0
        end do

c default lsh value...
        lsh=0

c Set up compare file name.
        char1='./coutput/compare_'
        open(unit=99,file='junk',status='unknown')
        write(99,*)n_atomic_number
        rewind(99)
        read(99,*)char2
        close(unit=99)
        call concat(char1,char2,charsum1)

c MHL added
        if (ish.le.nshells) then
c WANG
        if(ish.eq.-1)char2='_edens'
        if(ish.eq.0) char2='_napot'
        if(ish.eq.1) char2='_wave1'   
        if(ish.eq.2) char2='_wave2'   
        if(ish.eq.3) char2='_wave3'
        if(ish.eq.4) char2='_wave4'   
        if(ish.eq.5) char2='_wave5'
        if(ish.eq.6) char2='_wave6'
        if(ish.eq.7) char2='_wave7'
        if(ish.eq.8) char2='_wave8'
        if(ish.gt.8) then
        write(*,*)' '
        write(*,*)' The shell numbers exceed default setting'
        write(*,*)' Need to define more compare filename' 
        write(*,*)' ' 
        stop 
        endif

c MHL added
        else 
        if(ish.eq.nshells+1) char2='_wave1S'
        if(ish.eq.nshells+2) char2='_wave2S'
        if(ish.eq.nshells+3) char2='_wave3S'
        if(ish.eq.nshells+4) char2='_wave4S'
        if(ish.eq.nshells+5) char2='_wave5S'
        if(ish.eq.nshells+6) char2='_wave6S'
        if(ish.eq.nshells+7) char2='_wave7S'
        if(ish.eq.nshells+8) char2='_wave8S'
        end if
        if (ish.eq.nshells+1) then
        rewind(2)
        do i=1,6
        read(2,*)
        end do
        end if

        call concat(charsum1,char2,charsum2)
        char2='.dat'
        call concat(charsum2,char2,comparefile)
c
        comparefilename(ispec,ish)=comparefile
        open(unit=55,file=comparefile,status='unknown')

c WANG, density fit
        write(*,*)' '
        if (ish.eq.-1) 
     1  write(*,*)' Fit the electron density for the given species '
        if (ish.eq.-1) icounter=0
        if (ish.eq.-1) go to 8888

        if (ish.gt.0) then 
          read(2,*)lsh
          read(2,*)occupation
          read(2,*)rc
          read(2,527)wavefile
527       format(a25)
          read(2,*)
        end if
        if(ish.eq.0) write(*,*)' Fit neutral atom potential: ',wavefile
        if(ish.gt.0) write(*,*)' Fit wavefunction: ',wavefile

c now open wavefile and fit to gaussians...
        open(unit=3,file=wavefile,status='unknown')
c if it's a wavefunction...
        if (ish.gt.0) then
          read(3,*)
          read(3,*) natomicnumber
c atomic number check
          if (natomicnumber .ne. n_atomic_number) then
           write(*,*)' ERROR!!!!!!!!'
           write(*,*)' Atomic number in ',specfile
           write(*,*)' and ',wavefile
           write(*,*)' do not agree...'
           stop
          end if
          read(3,*)npnts
          read(3,*)rc_in_wavefile, rc_max_wf, xoccup_wf
c rc check 
         if (abs(rc_in_wavefile-rc) .gt. 1.E-2) then
           write(*,*)' ERROR!!!!!!!!'
           write(*,*)' rc in ',specfile
           write(*,*)' and ',wavefile
           write(*,*)' for shell number ',ish
           write(*,*)' do not agree!!!'
           write(*,*)' rc : ',rc_in_wavefile,rc
           stop
         end if
c occupation check
         if (xoccup_wf .ne. occupation ) then
           write(*,*)' occupation do not agree!!!'
           write(*,*)' occupation = ',xoccup_wf, occupation
           stop
         end if
c convert rc to angstroms
         rc=rc*bohr
         read(3,*) lsh_in_wavefile
c l number check
         if (lsh_in_wavefile .ne. lsh) then
           write(*,*)' l number in ',specfile
           write(*,*)' and ',wavefile
           write(*,*)' do not agree!!!'
           write(*,*)' l = ',lsh,lsh_in_wavefile
           stop
         end if
         read(3,*)(psi(ix),ix=1,npnts)

c MHL   
         if (ish.gt.nshells) then
           do ix=1,npnts
             psi(ix) = sqrt(psi(ix)**2.0d0)
           end do
         end if

c else it's the neutral atom potential, ...
         else

         read(3,*)
         read(3,*)natomicnumber
         if (natomicnumber.ne.n_atomic_number) then
           write(*,*)' ERROR!!!!!!!!'
           write(*,*)' Atomic number in ',specfile
           write(*,*)' and ',wavefile
           write(*,*)' do not agree...'
           stop
         end if

         read(3,*)rc
         rc=rc*bohr
         read(3,*)npnts
         read(3,*)
         do i=1,npnts
         read(3,*)rjunk,psi(i)
         end do

         end if !  end if for ish.gt.nshells

         close(unit=3)
c check normalization: 
c ----------------BEGIN INTEGRAL--------------------	
c integrate from A to B  (B > A).
        A=0.
        B=rc

c MHL TEST for checking <psi|n|psi>
c        if (ish.ge.1) then
c       dx=(B-A)/(npnts-1.)
c	sum=((psi(1)**4)*A**2+(psi(npnts)**4)*B**2)/2.
c        do i=1,npnts-2
c        x=A+i*dx
c        sum=sum+(psi(i+1))**4*(x**2)
c        end do
c        rnorm=sum*dx
c        write(*,*)'MHL TEST INTEGRAL(psi1**4*r**2)dr ',rnorm
c        end if
c MHL TEST end

c adjust the number of points to suit your problem.
        dx=(B-A)/(npnts-1.)
        sum=((psi(1)*A)**2+(psi(npnts)*B)**2)/2.
        do i=1,npnts-2
        x=A+i*dx
        sum=sum+(psi(i+1)*x)**2
        end do
        rnorm=sum*dx
c -----------------END INTEGRAL---------------------

c if it's a wavefunction it should be normalized...
        if (ish.ne.0) then
        write(*,*)' Input norm = ',rnorm
         if (abs(rnorm-1.).gt.1.E-2) then
           write(*,*)' ERROR!!!!!!!!!'
           write(*,*)' wavefunction not normalized.'
           write(*,*)' wavefunction = ',wavefile
           write(*,*)' norm = ',rnorm
           stop
         end if
        end if
c
8888    continue
        if(ish.eq.-1)write(*,*)' weight = ',weight
c
c now ready to fit to gaussians
c
        ial=0

c smallest error = big number
        error_smallest=1000000.

c ttest put in param.inc file...
c want potential to be very well fit...
        if(ish.eq.-1)tolerance=den_tolerance
        if (ish.gt.0)tolerance=psi_tolerance
        if (ish.eq.0)tolerance=rna_tolerance

c ttest new new enw enw en wne wn
        c_min=1./rc
        c_max=20./rc

c	ncpnts=100

c JKT ttest	 was 1000. trying 100
C 100 is probably safer...
C	d_min=(c_max-c_min)/1000.
        d_min=(c_max-c_min)/50.
        d_max=(c_max-c_min)/10.

c	ndpnts=10

c MAIN LOOP.
        do 3000 while (error_smallest.gt.tolerance .and. ial.lt.maxals)

c MHL for spherical approximation.
c In spherical approximation we use psi=sqrt(|R|**2)*Y_00 
c "gelementsGS" does the main calculation <psi|n_i|psi>, and it is
c designed to compute x*Gaussian not r*Gaussian.
c It's easier to get a new fit in "gauss_create" rather than changing
c "gelementsGS", so we set lsh=0 for wavefunciton fit for spherical
c approximation.

        if (ish.gt.nshells) then
           lsh=0
        end if
c MHL end

c number of gaussians to fit to.
        ial=ial+1

        error_smallest=1000000.

c	write(*,*)' Trying ',ial,' gaussians...'

        do ic=0,ncpnts
        do id=ndpnts,0,-1

        cnum=c_min+(c_max-c_min)*ic/ncpnts
        dnum=d_min+(d_max-d_min)*id/ndpnts

        do i=1,ial
           alpha(i)=(cnum+i*dnum)**2
        end do

c get f(i)= integral of exp(-al*r**2)*psi(r)*r**(lsh+2) dr
        do j=1,ial
c ---------------BEGIN INTEGRAL--------------------
c integrate from A to B  (B > A).
        A=0.
        B=rc
c adjust the number of points to suit your problem.
        dx=(B-A)/(npnts-1.)
c        sum=(f(A)+f(B))/2.

c WANG  Introduce weight factor for density fit to avoid
c       unphysical fitting result, e.g., n_fit < 0
c       In case of density fit, ish = -1 !

        if(ish.eq.-1)then
        sum=exp(-alpha(j)*B**2)*psi(npnts)*
     1       (weight+B**2)*B**lsh/2.
        else
        sum=exp(-alpha(j)*B**2)*psi(npnts)*B**(lsh+2)/2.
        endif

        do i=1,npnts-2
        x=A+i*dx
c        sum=sum+f(x)

c WANG  Introduce weight factor for density fit to avoid
c       unphysical fitting result, e.g., n_fit < 0
c       In case of density fit, ish = -1 !

        if(ish.eq.-1)then
        sum=sum+exp(-alpha(j)*x**2)*psi(i+1)*
     1       (weight+x**2)*x**lsh
        else
        sum=sum+exp(-alpha(j)*x**2)*psi(i+1)*x**(lsh+2)
        endif

        end do

        sum=sum*dx
        f(j)=sum
c ----------------END INTEGRAL---------------------
        end do

c now get G(i,j) = integral of g(-(ali+alj)*r**2)*r**2 dr
c ttest must set this up so that it does it for any lsh!!!!!!!
c right now it's only lsh=0 !!!!!!!!
        if (lsh.eq.0) then
        do i=1,ial
        do j=1,ial

c WANG  Introduce weight factor for density fit to avoid
c       unphysical fitting result, e.g., n_fit < 0
c       In case of density fit, ish = -1 !

        if(ish.eq.-1)then
        G(i,j)=sqrt(pi/(alpha(i)+alpha(j)))/
     1         (alpha(i)+alpha(j))/4.+
     2         0.5*weight*sqrt(pi/(alpha(i)+alpha(j)))
 
        else
        G(i,j)=sqrt(pi/(alpha(i)+alpha(j)))/
     1         (alpha(i)+alpha(j))/4.
        endif 
        end do
        end do
        else if (lsh.eq.1) then
        do i=1,ial
        do j=1,ial
        G(i,j)=3.*sqrt(pi/(alpha(i)+alpha(j)))/
     1         (alpha(i)+alpha(j))**2/8.
        end do
        end do
        else if (lsh.eq.2) then
        do i=1,ial
        do j=1,ial
        G(i,j)=15.*sqrt(pi/(alpha(i)+alpha(j)))/
     1         (alpha(i)+alpha(j))**3/16.
        end do
        end do
        end if

C	do i=1,ial
C	do j=1,ial
C	Gtemp(i,j)=G(i,j)
C	end do
C	end do

c now invert G
        call INVERT(G,ial,maxmaxals,M,C)

c get coefs...
        garbage=.false.
        do i=1,ial
        sum=0.
        do j=1,ial
        sum=sum+G(i,j)*f(j)
        end do
        coef(i)=sum
        if (coef(i).gt.rmax_coef*rnorm) then
           garbage=.true.
        end if
        end do

c now compute error

c -------------------BEGIN INTEGRAL--------------------
c integrate from A to B  (B > A).
        A=0.
        B=rc
c adjust the number of points to suit your problem.
        dx=(B-A)/(npnts-1.)
c        sum=(f(A)+f(B))/2.
        fB=0.
        do i=1,ial
        fB=fB+coef(i)*exp(-alpha(i)*B**2)*B**lsh
        end do
        sum=((fB-psi(npnts))*B)**2/2.

        gA=0.
        do i=1,ial
          gA=gA+coef(i)*exp(-alpha(i)*A**2)*A**lsh
        end do
        sum2=(gA+fB)/2.

        do i=1,npnts-2
        x=A+i*dx
c        sum=sum+f(x)
        fx=0.
        do j=1,ial
        fx=fx+coef(j)*exp(-alpha(j)*x**2)*x**lsh
        end do

c WANG  Introduce weight factor for density fit to avoid
c       unphysical fitting result, e.g., n_fit < 0

        if(ish.eq.-1)then
        sum=sum+(fx-psi(i+1))**2*(weight+x**2)
        else
        sum=sum+((fx-psi(i+1))*x)**2
        endif

        sum2=sum2+(fx*x)**2
        end do
        sum=sum*dx
        sum2=sum2*dx

C sum2 is the amount of gauss wavefunction inside rc.
C sum3 is the amount of gauss wavefunction outside of rc.
        sum3=0.
        do i=npnts-2,3*npnts
        x=A+i*dx
        fx=0.
        do j=1,ial
        fx=fx+coef(j)*exp(-alpha(j)*x**2)*x**lsh
        end do
        sum3=sum3+(fx*x)**2
        end do
        sum3=sum3*dx

C	sumtest=0.
C	do i=1,ial
C	do j=1,ial
C	sumtest=sumtest+coef(i)*Gtemp(i,j)*coef(j)
C	end do
C	end do
c weight it...
        sum=(sum+sum3)/rnorm
c --------------------END INTEGRAL---------------------

        if (sum.lt.error_smallest .and. .not.garbage) then
           error_smallest=sum
           c_best=cnum
           d_best=dnum
           amount_beyond_rc=sum3/sum2
C ttest
C	write(*,*)' sum2,sum3,sum3/sum2 = ',sum2,sum3,sum3/sum2
        end if

c id 
        end do
c ic
        end do
c Electron Density
        if(ish.eq.-1)then
c Use the RMS value past Rc.
        amount_beyond_rc=sqrt(amount_beyond_rc)
        write(*,910)ial,error_smallest,amount_beyond_rc
910     format('  Ngauss:',i3,' Error=',f12.7,' RMS  past Rc=',f12.7)
        end if
c NA potential
        if(ish.eq.0)then
c Use the RMS value past Rc.
        amount_beyond_rc=sqrt(amount_beyond_rc)
        write(*,911)ial,error_smallest,amount_beyond_rc
911     format('  Ngauss:',i3,' Error=',f12.7,' RMS  past Rc=',f12.7)
        end if
        if(ish.gt.0)then
        write(*,912)ial,error_smallest,amount_beyond_rc
912     format('  Ngauss:',i3,' Error=',f12.7,' Norm past Rc=',f12.7)
        end if

c       write(*,*)' error =',error_smallest
c	write(*,*)' fraction beyond rc = ',amount_beyond_rc
c	write(*,*)
c ial
3000    continue

c END OF MAIN LOOP

        if (amount_beyond_rc.gt.amount_beyond_rc_max) then
        write(*,*)
        write(*,*)' WARNING WARNING WARNING!!!'
        write(*,*)
        warnings=.true.
        write(*,*)' Fraction of gauss fit beyond rc may be large.'
        write(*,*)' The fracion is ',amount_beyond_rc
        write(*,*)' If this is significantly bigger than'
        write(*,*)amount_beyond_rc_max
        write(*,*)' you should beware.'
        write(*,*)
        end if

c best alphas...
        do i=1,ial
          alpha(i)=(c_best+i*d_best)**2
        end do

c ====================GET COEFS==========================
c got the best alpha's. Now get the coefs and we're done.
c get f(i)= integral of exp(-al*r**2)*psi(r)*r**(lsh+2) dr
        do j=1,ial
c ---------------BEGIN INTEGRAL--------------------
c integrate from A to B  (B > A).
        A=0.
        B=rc
c adjust the number of points to suit your problem.
        dx=(B-A)/(npnts-1.)
c        sum=(f(A)+f(B))/2.

c WANG  Introduce weight factor for density fit to avoid
c       unphysical fitting result, e.g., n_fit < 0
c       In case of density fit, lsh=0 !

        if(ish.eq.-1)then
        sum=exp(-alpha(j)*B**2)*psi(npnts)*
     1      (weight+B**2)*B**lsh/2.
        else
        sum=exp(-alpha(j)*B**2)*psi(npnts)*B**(lsh+2)/2.
        endif
        do i=1,npnts-2
        x=A+i*dx
c        sum=sum+f(x)

c WANG  Introduce weight factor for density fit to avoid
c       unphysical fitting result, e.g., n_fit < 0
c       In case of density fit, lsh=0 !

        if(ish.eq.-1)then
        sum=sum+exp(-alpha(j)*x**2)*psi(i+1)*
     1      (weight+x**2)*x**lsh
        else
        sum=sum+exp(-alpha(j)*x**2)*psi(i+1)*x**(lsh+2)
        endif
        end do
        sum=sum*dx
        f(j)=sum

c ----------------END INTEGRAL---------------------
        end do

c now get G(i,j) = integral of g(-(ali+alj)*r**2)*r**2 dr
c ttest must set this up so that it does it for any lsh!!!!!!!
c right now it's only lsh=0 !!!!!!!!
        if (lsh.eq.0) then
        do i=1,ial
        do j=1,ial

c WANG  Introduce weight factor for density fit to avoid 
c       unphysical fitting result, e.g., n_fit < 0
c       In case of density fit, lsh=0 !

        if(ish.eq.-1)then
        G(i,j)=sqrt(pi/(alpha(i)+alpha(j)))/
     1         (alpha(i)+alpha(j))/4.+
     2         0.5*weight*sqrt(pi/(alpha(i)+alpha(j)))

        else
        G(i,j)=sqrt(pi/(alpha(i)+alpha(j)))/
     1         (alpha(i)+alpha(j))/4.
        endif
        end do
        end do
        else if (lsh.eq.1) then
        do i=1,ial
        do j=1,ial
          G(i,j)=3.*sqrt(pi/(alpha(i)+alpha(j)))/
     1         (alpha(i)+alpha(j))**2/8.
        end do
        end do
        else if (lsh.eq.2) then
        do i=1,ial
        do j=1,ial
          G(i,j)=15.*sqrt(pi/(alpha(i)+alpha(j)))/
     1                       (alpha(i)+alpha(j))**3/16.
        end do
        end do
        end if

c now invert G
        call INVERT(G,ial,maxmaxals,M,C)

c get coefs...
        do i=1,ial
        sum=0.
        do j=1,ial
        sum=sum+G(i,j)*f(j)
        end do
        coef(i)=sum
        end do
c ================END GET COEFS==========================

c8828    continue
c WANG  write out original and fitted function for comparison

c	iunphys=0
        do i=0,1.1*npnts
        x=i*dx
        fx=0.
        do j=1,ial
        fx=fx+coef(j)*exp(-alpha(j)*x**2)*x**lsh
        end do
c WANG
        if(ish.eq.-1.and.fx.lt.0.)then
        write(*,*)' '
        write(*,*)' Unphysical density fitting, n_fit < 0, detected !!!'
c        write(*,*)' Adjust weight factor and rerun program'
c        weight=0.25d0+weight
c	icounter=icounter+1
c 10 is the maximum number of iterations.
c        if(icounter.gt.10)then
c        write(*,*)'10 is the maximum number of iterations. We stop here'
c        write(*,*)' Please check your density fit (compare) to see '
c        write(*,*)' how negative the density is. Be careful.'
c        write(*,*)' This is the TOTAL density'
c	go to 8889
c	end if
c	go to 8888
c HAO add at April 26, 2005
        write(*,*)' set negative density to the original value '
        fx = psi (i)
        endif
c8889	continue
        weight = 0.
c WANG
        if (i+1.le.npnts) then
        write(55,*)x,psi(i+1),fx
        else
        write(55,*)x,0.,fx
        end if
        end do
        write(55,*)

c write to gauss.dat file...
c	write(*,*)
c	write(*,*)' writing to gauss.dat file.'
        if (ish.eq.-1) then
          write(9,*)' Number of gaussians for electron density'
        endif
        if(ish.eq.0) then
          write(9,*)' Number of gaussians for NA potential'
        endif
c Old statement	if(ish.gt.0) then
c MHL
        if (ish.gt.0 .and. ish.le.nshells) then
        write(9,*)' Number of gaussians for wavefunction of shell ',ish
        endif
        if (ish.gt.nshells) then
        write(9,*)' Number of gaussians for spher.approx. shell',       
     1            ish-nshells
        end if
        write(9,*) ial
        write(9,*)' alphas'
        write(9,*) (alpha(ialpha),ialpha=1,ial)
        write(9,*)' coefs'
        write(9,*) (coef(icoef),icoef=1,ial)
c	write(*,*)
        close(unit=55)
999     continue
c ish
1500    continue
        close(unit=2)
c ==============================================================
c MHL added for getting electron density for each shell.
        do 2500 ish=0,nshells-1
        icounter=0

c HAO add at May 8, 2005, initialize coefficents and psi
        do ial = 1, maxals
           alpha (ial) = 0.d0
            coef (ial) = 0.d0
        end do

        call gauss_density(ish, specfile, rc, npnts, rnorm, psi)
        lsh=0
c Set up compare file name.
        char1='./coutput/compare_'
        open(unit=99,file='junk',status='unknown')
        write(99,*)n_atomic_number
        rewind(99)
        read(99,*)char2
        close(unit=99)
        call concat(char1,char2,charsum1)
c WANG
        if(ish.eq.0) char2='_edens1'
        if(ish.eq.1) char2='_edens2'
        if(ish.eq.2) char2='_edens3'
        if(ish.eq.3) char2='_edens4'
        if(ish.eq.4) char2='_edens5'
        if(ish.eq.5) char2='_edens6'
        if(ish.eq.6) char2='_edens7'
        if(ish.eq.7) char2='_edens8'
        if(ish.gt.7) then
        write(*,*)' '
        write(*,*)' The shell numbers exceed default setting'
        write(*,*)' Need to define more compare filename'
        write(*,*)' '
        stop
        endif

        call concat(charsum1,char2,charsum2)
        char2='.dat'
        call concat(charsum2,char2,comparefile)
c
        comparefilename(ispec,ish+nshells*2+1)=comparefile
        open(unit=55,file=comparefile,status='unknown')

        write(*,*)' '
        write(*,*)' Fit Electron Density for shell: ',ish+1
        write(*,*)
        write(*,*)' Input rnorm = ',rnorm
5555    continue
c        write(*,*)' weight = ',weight
c
c now ready to fit to gaussians

        ial=0

c smallest error = big number
        error_smallest=1000000.

c ttest put in param.inc file...
c want potential to be very well fit...
        tolerance=den_tolerance

c ttest new new enw enw en wne wn
        c_min=1./rc
        c_max=20./rc

c       ncpnts=100

c JKT ttest      was 1000. trying 100
C 100 is probably safer...
C       d_min=(c_max-c_min)/1000.
        d_min=(c_max-c_min)/50.
        d_max=(c_max-c_min)/10.

c       ndpnts=10

c MAIN LOOP.
        do 1250 while (error_smallest.gt.tolerance .and. ial.lt.maxals)

c number of gaussians to fit to.
        ial=ial+1

        error_smallest=1000000.

c       write(*,*)' Trying ',ial,' gaussians...'

        do ic=0,ncpnts
        do id=ndpnts,0,-1

        cnum=c_min+(c_max-c_min)*ic/ncpnts
        dnum=d_min+(d_max-d_min)*id/ndpnts

        do i=1,ial
        alpha(i)=(cnum+i*dnum)**2
        end do

c get f(i)= integral of exp(-al*r**2)*psi(r)*r**(lsh+2) dr
        do j=1,ial
c ---------------BEGIN INTEGRAL--------------------
c integrate from A to B  (B > A).
        A=0.
        B=rc
c adjust the number of points to suit your problem.
        dx=(B-A)/(npnts-1.)
c        sum=(f(A)+f(B))/2.

c WANG  Introduce weight factor for density fit to avoid
c       unphysical fitting result, e.g., n_fit < 0
c       In case of density fit, lsh=0 !

        sum=exp(-alpha(j)*B**2)*psi(npnts)*
     1       (weight+B**2)*B**lsh/2.
        do i=1,npnts-2
        x=A+i*dx
c        sum=sum+f(x)

c WANG  Introduce weight factor for density fit to avoid
c       unphysical fitting result, e.g., n_fit < 0
c       In case of density fit, lsh=0 !

        sum=sum+exp(-alpha(j)*x**2)*psi(i+1)*
     1       (weight+x**2)*x**lsh
        end do

        sum=sum*dx
        f(j)=sum
c ----------------END INTEGRAL---------------------
        end do

c now get G(i,j) = integral of g(-(ali+alj)*r**2)*r**2 dr
c ttest must set this up so that it does it for any lsh!!!!!!!
c right now it's only lsh=0 !!!!!!!!
        do i=1,ial
        do j=1,ial

c WANG  Introduce weight factor for density fit to avoid
c       unphysical fitting result, e.g., n_fit < 0
c       In case of density fit, lsh=0 !

        G(i,j)=sqrt(pi/(alpha(i)+alpha(j)))/
     1         (alpha(i)+alpha(j))/4.+
     2         0.5*weight*sqrt(pi/(alpha(i)+alpha(j)))

        end do
        end do

C       do i=1,ial
C       do j=1,ial
C       Gtemp(i,j)=G(i,j)
C       end do
C       end do

c now invert G
        call INVERT(G,ial,maxmaxals,M,C)

c get coefs...
        garbage=.false.
        do i=1,ial
        sum=0.
        do j=1,ial
        sum=sum+G(i,j)*f(j)
        end do
        coef(i)=sum
        if (coef(i).gt.rmax_coef*rnorm) then
          garbage=.true.
        end if
        end do

c now compute error

c -------------------BEGIN INTEGRAL--------------------
c integrate from A to B  (B > A).
        A=0.
        B=rc
c adjust the number of points to suit your problem.
        dx=(B-A)/(npnts-1.)
c        sum=(f(A)+f(B))/2.
        fB=0.
        do i=1,ial
        fB=fB+coef(i)*exp(-alpha(i)*B**2)*B**lsh
        end do
        sum=((fB-psi(npnts))*B)**2/2.

        gA=0.
        do i=1,ial
        gA=gA+coef(i)*exp(-alpha(i)*A**2)*A**lsh
        end do
        sum2=(gA+fB)/2.

        do i=1,npnts-2
        x=A+i*dx
c        sum=sum+f(x)
        fx=0.
        do j=1,ial
        fx=fx+coef(j)*exp(-alpha(j)*x**2)*x**lsh
        end do

c WANG  Introduce weight factor for density fit to avoid
c       unphysical fitting result, e.g., n_fit < 0

        sum=sum+(fx-psi(i+1))**2*(weight+x**2)
        sum2=sum2+(fx*x)**2
        end do
        sum=sum*dx
        sum2=sum2*dx

C sum2 is the amount of gauss wavefunction inside rc.
C sum3 is the amount of gauss wavefunction outside of rc.
        sum3=0.
        do i=npnts-2,3*npnts
        x=A+i*dx
        fx=0.
        do j=1,ial
        fx=fx+coef(j)*exp(-alpha(j)*x**2)*x**lsh
        end do
        sum3=sum3+(fx*x)**2
        end do
        sum3=sum3*dx

C       sumtest=0.
C       do i=1,ial
C       do j=1,ial
C       sumtest=sumtest+coef(i)*Gtemp(i,j)*coef(j)
C       end do
C       end do

c weight it...
        sum=(sum+sum3)/rnorm
c --------------------END INTEGRAL---------------------

        if (sum.lt.error_smallest .and. .not.garbage) then
           error_smallest=sum
           c_best=cnum
           d_best=dnum
           amount_beyond_rc=sum3/sum2
C ttest
C       write(*,*)' sum2,sum3,sum3/sum2 = ',sum2,sum3,sum3/sum2
        end if

c id
        end do
c ic
        end do

c Use the RMS value past Rc.
        amount_beyond_rc=sqrt(amount_beyond_rc)
        write(*,910)ial,error_smallest,amount_beyond_rc
c 910     format('  Ngauss:',i3,' Error=',f12.7,' RMS  past Rc=',f12.7)

c        write(*,*)' error =',error_smallest
c       write(*,*)' fraction beyond rc = ',amount_beyond_rc
c       write(*,*)
c ial
1250    continue

c END OF MAIN LOOP

        if (amount_beyond_rc.gt.amount_beyond_rc_max) then
        write(*,*)
        write(*,*)' WARNING WARNING WARNING!!!'
        write(*,*)
        warnings=.true.
        write(*,*)' Fraction of gauss fit beyond rc may be large.'
        write(*,*)' The fracion is ',amount_beyond_rc
        write(*,*)' If this is significantly bigger than'
        write(*,*)amount_beyond_rc_max
        write(*,*)' you should beware.'
        write(*,*)
        end if

c best alphas...
        do i=1,ial
        alpha(i)=(c_best+i*d_best)**2
        end do

c ====================GET COEFS==========================
c got the best alpha's. Now get the coefs and we're done.
c get f(i)= integral of exp(-al*r**2)*psi(r)*r**(lsh+2) dr
        do j=1,ial
c ---------------BEGIN INTEGRAL--------------------
c integrate from A to B  (B > A).
        A=0.
        B=rc
c adjust the number of points to suit your problem.
        dx=(B-A)/(npnts-1.)
c        sum=(f(A)+f(B))/2.

c WANG  Introduce weight factor for density fit to avoid
c       unphysical fitting result, e.g., n_fit < 0
c       In case of density fit, lsh=0 !

        sum=exp(-alpha(j)*B**2)*psi(npnts)*
     1      (weight+B**2)*B**lsh/2.
        do i=1,npnts-2
        x=A+i*dx
c        sum=sum+f(x)

c WANG  Introduce weight factor for density fit to avoid
c       unphysical fitting result, e.g., n_fit < 0
c       In case of density fit, lsh=0 !

        sum=sum+exp(-alpha(j)*x**2)*psi(i+1)*
     1      (weight+x**2)*x**lsh
        end do
        sum=sum*dx
        f(j)=sum

c ----------------END INTEGRAL---------------------
       end do

c        now get G(i,j) = integral of g(-(ali+alj)*r**2)*r**2 dr
c ttest must set this up so that it does it for any lsh!!!!!!!
c right now it's only lsh=0 !!!!!!!!
        do i=1,ial
        do j=1,ial

c WANG  Introduce weight factor for density fit to avoid
c       unphysical fitting result, e.g., n_fit < 0
c       In case of density fit, lsh=0 !

        G(i,j)=sqrt(pi/(alpha(i)+alpha(j)))/
     1         (alpha(i)+alpha(j))/4.+
     2         0.5*weight*sqrt(pi/(alpha(i)+alpha(j)))

        end do
        end do

c now invert G
        call INVERT(G,ial,maxmaxals,M,C)

c get coefs...
        do i=1,ial
        sum=0.
        do j=1,ial
        sum=sum+G(i,j)*f(j)
        end do
        coef(i)=sum
        end do
c ================END GET COEFS==========================
c WANG  write out original and fitted function for comparison
c	iunphys=0
        do i=0,1.1*npnts
        x=i*dx
        fx=0.
        do j=1,ial
        fx=fx+coef(j)*exp(-alpha(j)*x**2)*x**lsh
        end do
c WANG
        if(fx.lt.0.)then
        write(*,*)' '
        write(*,*)' Unphysical density fitting, n_fit < 0, detected !!!'
c	write(*,*)' Adjust weight factor and rerun program'
c        weight=0.25d0+weight
c	icounter=icounter+1
c 10 is the maximum number of iterations.
c        if(icounter.gt.10)then
c        write(*,*)'10 is the maximum number of iterations. We stop here'
c        write(*,*)' Please check your density fit (compare) to see '
c        write(*,*)' how negative the density is. Be careful.'
c        write(*,*)' The shell is',ish
c	goto 5559
c	end if
c        goto 5555
c HAO add at April 26, 2005
        write(*,*)' set negative density to zero '
        fx = psi(i)
        endif
c5559	continue
c ==============================================================
        weight = 0.
c WANG
        if (i+1.le.npnts) then
        write(55,*)x,psi(i+1),fx
        else
        write(55,*)x,0.,fx
        end if
        end do
        write(55,*)
c MHL TEST
        write(*,*)'ish, weight ', ish,weight

c write to gauss.dat file...
c       write(*,*)
c       write(*,*)' writing to gauss.dat file.'
        write(9,*)' Number of gaussians for edens of shell ',ish+1
        write(9,*) ial
        write(9,*)' alphas'
        write(9,*) (alpha(ialpha),ialpha=1,ial)
        write(9,*)' coefs'
        write(9,*) (coef(icoef),icoef=1,ial)
c       write(*,*)

        close(unit=55)
c 999     continue

c ish
2500    continue
c MHL end added
c =========================================================================
c HAO add for fitting neutral atom potential by shell
c
        open(unit=42,file=specfile,status='unknown')
c
        read(42,*)
        read(42,*) n_atomic_number
        read(42,*)
        read(42,*)
        read(42,*)
        read(42,*) nshells_na
c l number check
         if (nshells_na .ne. nshells) then
         write (*,*) ' shell number in is not consistent with what we '
         write (*,*) ' have before : ', nshells, nshells_na
         stop
         end if

        do 4500 ish = 0, nshells_na-1
        icounter=0

        read(42,*)
        read(42,*)
        read(42,*)
        read(42,*)
        read(42,527)wavefile

c now open wavefile and fit to gaussians...
        open(unit = 43,file = wavefile, status = 'unknown')

         read(43,*)
         read(43,*) natomicnumber
         if (natomicnumber .ne. n_atomic_number) then
           write(*,*)' ERROR!!!!!!!!'
           write(*,*)' Atomic number in ',specfile
           write(*,*)' and ',wavefile
           write(*,*)' do not agree...'
           stop
         end if

         read(43,*) rc
         rc = rc*bohr
         read(43,*) npnts
         do i = 1, npnts
         read(43,*) rjunk, psi(i)
         end do

         R_na = rjunk
         rc = R_na

! the lagerest r in NA potential for given shell

         do i = 1, npnts
         psi(i) = psi (i) - 1./R_na
         end do

         close(unit=43)

c HAO add at May 8, 2005, initialize coefficents and psi
        do ial = 1, maxals
           alpha (ial) = 0.d0
            coef (ial) = 0.d0
        end do
c
        lsh=0
c Set up compare file name.
        char1='./coutput/compare_'
        open(unit=99,file='junk',status='unknown')
        write(99,*)n_atomic_number
        rewind(99)
        read(99,*)char2
        close(unit=99)
        call concat(char1,char2,charsum1)
c WANG
        if(ish.eq.0) char2='_napot1'
        if(ish.eq.1) char2='_napot2'
        if(ish.eq.2) char2='_napot3'
        if(ish.eq.3) char2='_napot4'
        if(ish.eq.4) char2='_napot5'
        if(ish.eq.5) char2='_napot6'
        if(ish.eq.6) char2='_napot7'
        if(ish.eq.7) char2='_napot8'
        if(ish.gt.7) then
        write(*,*)' '
        write(*,*)' The shell numbers exceed default setting'
        write(*,*)' Need to define more compare filename'
        write(*,*)' '
        stop
        endif

        call concat(charsum1,char2,charsum2)
        char2='.dat'
        call concat(charsum2,char2,comparefile)
c
        comparefilename(ispec,ish+nshells*3+1)=comparefile
        open(unit=55,file=comparefile,status='unknown')

        write(*,*) ' '
        write(*,*) ' Fit neutral atom potential for shell: ', ish+1
        write(*,*)
c
c now ready to fit to gaussians
        ial=0

c smallest error = big number
        error_smallest=1000000.

c ttest put in param.inc file...
c want potential to be very well fit...
        tolerance = rna_tolerance

c ttest new new 
        c_min=1./rc
        c_max=20./rc

c       ncpnts=100

c JKT ttest      was 1000. trying 100
C 100 is probably safer...
C       d_min=(c_max-c_min)/1000.
        d_min=(c_max-c_min)/50.
        d_max=(c_max-c_min)/10.

c       ndpnts=10

c MAIN LOOP.
        do 1450 while (error_smallest.gt.tolerance .and. ial.lt.maxals)

c number of gaussians to fit to.
        ial=ial+1

        error_smallest=1000000.

c       write(*,*)' Trying ',ial,' gaussians...'

        do ic=0,ncpnts
        do id=ndpnts,0,-1

        cnum=c_min+(c_max-c_min)*ic/ncpnts
        dnum=d_min+(d_max-d_min)*id/ndpnts

        do i=1,ial
        alpha(i)=(cnum+i*dnum)**2
        end do

c get f(i)= integral of exp(-al*r**2)*psi(r)*r**(lsh+2) dr
        do j=1,ial
c ---------------BEGIN INTEGRAL--------------------
c integrate from A to B  (B > A).
        A=0.
        B=rc
c adjust the number of points to suit your problem.
        dx=(B-A)/(npnts-1.)
c        sum=(f(A)+f(B))/2.

c WANG  Introduce weight factor for density fit to avoid
c       unphysical fitting result, e.g., n_fit < 0
c       In case of density fit, lsh=0 !

        sum=exp(-alpha(j)*B**2)*psi(npnts)*B**(2+lsh)/2.
        do i=1,npnts-2
        x=A+i*dx
c        sum=sum+f(x)

c WANG  Introduce weight factor for density fit to avoid
c       unphysical fitting result, e.g., n_fit < 0
c       In case of density fit, lsh=0 !

        sum=sum+exp(-alpha(j)*x**2)*psi(i+1)*x**(2+lsh)
        end do

        sum=sum*dx
        f(j)=sum
c ----------------END INTEGRAL---------------------
        end do

c now get G(i,j) = integral of g(-(ali+alj)*r**2)*r**2 dr
c ttest must set this up so that it does it for any lsh!!!!!!!
c right now it's only lsh=0 !!!!!!!!
        do i=1,ial
        do j=1,ial

c WANG  Introduce weight factor for density fit to avoid
c       unphysical fitting result, e.g., n_fit < 0
c       In case of density fit, lsh=0 !

        G(i,j)=sqrt(pi/(alpha(i)+alpha(j)))/
     1         (alpha(i)+alpha(j))/4.

        end do
        end do

C       do i=1,ial
C       do j=1,ial
C       Gtemp(i,j)=G(i,j)
C       end do
C       end do

c now invert G
        call INVERT(G,ial,maxmaxals,M,C)

c get coefs...
        garbage=.false.
        do i=1,ial
        sum=0.
        do j=1,ial
        sum=sum+G(i,j)*f(j)
        end do
        coef(i)=sum
        if (coef(i).gt.rmax_coef*rnorm) then
          garbage=.true.
        end if
        end do

c now compute error

c -------------------BEGIN INTEGRAL--------------------
c integrate from A to B  (B > A).
        A=0.
        B=rc
c adjust the number of points to suit your problem.
        dx=(B-A)/(npnts-1.)
c        sum=(f(A)+f(B))/2.
        fB=0.
        do i=1,ial
        fB=fB+coef(i)*exp(-alpha(i)*B**2)*B**lsh
        end do
        sum=((fB-psi(npnts))*B)**2/2.

        gA=0.
        do i=1,ial
        gA=gA+coef(i)*exp(-alpha(i)*A**2)*A**lsh
        end do
        sum2=(gA+fB)/2.

        do i=1,npnts-2
        x=A+i*dx
c        sum=sum+f(x)
        fx=0.
        do j=1,ial
        fx=fx+coef(j)*exp(-alpha(j)*x**2)*x**lsh
        end do

c WANG  Introduce weight factor for density fit to avoid
c       unphysical fitting result, e.g., n_fit < 0

        sum=sum+(fx-psi(i+1))**2*x**2
        sum2=sum2+(fx*x)**2
        end do
        sum=sum*dx
        sum2=sum2*dx

C sum2 is the amount of gauss wavefunction inside rc.
C sum3 is the amount of gauss wavefunction outside of rc.
        sum3=0.
        do i=npnts-2,3*npnts
        x=A+i*dx
        fx=0.
        do j=1,ial
        fx=fx+coef(j)*exp(-alpha(j)*x**2)*x**lsh
        end do
        sum3=sum3+(fx*x)**2
        end do
        sum3=sum3*dx

c weight it...
        sum=(sum+sum3)/rnorm
c --------------------END INTEGRAL---------------------

        if (sum.lt.error_smallest .and. .not.garbage) then
           error_smallest=sum
           c_best=cnum
           d_best=dnum
           amount_beyond_rc=sum3/sum2
C ttest
C       write(*,*)' sum2,sum3,sum3/sum2 = ',sum2,sum3,sum3/sum2
        end if

c id
        end do
c ic
        end do

c Use the RMS value past Rc.
        amount_beyond_rc=sqrt(amount_beyond_rc)
        write(*,910)ial,error_smallest,amount_beyond_rc
c 910     format('  Ngauss:',i3,' Error=',f12.7,' RMS  past Rc=',f12.7)

c        write(*,*)' error =',error_smallest
c       write(*,*)' fraction beyond rc = ',amount_beyond_rc
c       write(*,*)
c ial
1450    continue

c END OF MAIN LOOP

        if (amount_beyond_rc.gt.amount_beyond_rc_max) then
        write(*,*)
        write(*,*)' WARNING WARNING WARNING!!!'
        write(*,*)
        warnings=.true.
        write(*,*)' Fraction of gauss fit beyond rc may be large.'
        write(*,*)' The fracion is ',amount_beyond_rc
        write(*,*)' If this is significantly bigger than'
        write(*,*)amount_beyond_rc_max
        write(*,*)' you should beware.'
        write(*,*)
        end if

c best alphas...
        do i=1,ial
        alpha(i)=(c_best+i*d_best)**2
        end do

c ====================GET COEFS==========================
c got the best alpha's. Now get the coefs and we're done.
c get f(i)= integral of exp(-al*r**2)*psi(r)*r**(lsh+2) dr
        do j=1,ial
c ---------------BEGIN INTEGRAL--------------------
c integrate from A to B  (B > A).
        A=0.
        B=rc
c adjust the number of points to suit your problem.
        dx=(B-A)/(npnts-1.)
c        sum=(f(A)+f(B))/2.

c WANG  Introduce weight factor for density fit to avoid
c       unphysical fitting result, e.g., n_fit < 0
c       In case of density fit, lsh=0 !

        sum=exp(-alpha(j)*B**2)*psi(npnts)*
     1      (B**2)*B**lsh/2.
        do i=1,npnts-2
        x=A+i*dx
c        sum=sum+f(x)

c WANG  Introduce weight factor for density fit to avoid
c       unphysical fitting result, e.g., n_fit < 0
c       In case of density fit, lsh=0 !

        sum=sum+exp(-alpha(j)*x**2)*psi(i+1)*x**(2+lsh)
        end do
        sum=sum*dx
        f(j)=sum

c ----------------END INTEGRAL---------------------
       end do

c        now get G(i,j) = integral of g(-(ali+alj)*r**2)*r**2 dr
c ttest must set this up so that it does it for any lsh!!!!!!!
c right now it's only lsh=0 !!!!!!!!
        do i=1,ial
        do j=1,ial

c WANG  Introduce weight factor for density fit to avoid
c       unphysical fitting result, e.g., n_fit < 0
c       In case of density fit, lsh=0 !

        G(i,j)=sqrt(pi/(alpha(i)+alpha(j)))/
     1         (alpha(i)+alpha(j))/4.

        end do
        end do

c now invert G
        call INVERT(G,ial,maxmaxals,M,C)

c get coefs...
        do i=1,ial
        sum=0.
        do j=1,ial
        sum=sum+G(i,j)*f(j)
        end do
        coef(i)=sum
        end do
c ================END GET COEFS==========================
c WANG  write out original and fitted function for comparison
c	iunphys=0
        do i=0,1.1*npnts
        x=i*dx
        fx=0.
        do j=1,ial
        fx=fx+coef(j)*exp(-alpha(j)*x**2)*x**lsh
        end do
c ==============================================================
        if (i+1.le.npnts) then
        write(55,*)x,psi(i+1),fx
        else
        write(55,*)x,0.,fx
        end if
        end do
        write(55,*)

c write to gauss.dat file...
        write(9,*)' Number of gaussians for napot of shell ',ish + 1
        write(9,*) ial
        write(9,*)' alphas '
        write(9,*) (alpha(ialpha),ialpha=1,ial)
        write(9,*)' coefs '
        write(9,*) (coef(icoef),icoef=1,ial)
        write(9,*) ' R_na '
        write(9,*) R_na
! the lagerest r in NA potential for given shell
c       write(*,*)

        close(unit=55)

c ish
4500    continue

        close (42)
c HAO end added for neutral atom potential
c ==============================================================
c
1000    continue  ! end loop of ispec

        write(9,*)' ===================================='

c OFS ofs added
        close(unit=1)
        close(unit=9)
c OFS end added

        write(*,*)'  '
        write(*,*)' =====================================',
     1            '===================='
        write(*,*)' '
        write(*,*)' Done complete.'
        write(*,*)' '
        write(*,*)' Summary:'
        write(*,*)' Accuracy level:', naccuracy
        write(*,*)' nspec=',nspecs
        write(*,610)(iatomic(i),i=1,nspecs)
610     format('  Atomic numbers are:',I4)
        write(*,*)' '
        write(*,*)' You can check the fits by plotting the data'
        write(*,*)' in compare file: ./coutput/compare_X.dat.'
        write(*,*)'    1st column = r'
        write(*,*)'    2nd column = exact or original'
        write(*,*)'    3rd column = fit with guassians'
        write(*,*)' The compare file names are:'
        do i=1,nspecs
c
        do issh=ish0,numshl(i)*4

        write(*,59)comparefilename(i,issh)
59      format('         ',a60)

        end do
        write(*,*)
        end do

        write(*,*)' compare file name+S is for spherical approximation'
        write(*,*)' Only gauss.dat file is needed by fireball '
        write(*,*)' which is written to ./coutput/gauss.dat'

        if (warnings) then
        write(*,*)
        write(*,*)' There were one or more WARNINGS triggered'
        write(*,*)' during your run. Please beware.'
        end if
        write(*,*)
        write(*,*)' Bye from gauss_create!'
        write(*,*)' '
        write(*,*)' '
c
        return
c        stop
        end
c ====================================================================
