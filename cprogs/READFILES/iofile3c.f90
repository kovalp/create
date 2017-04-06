! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! University of Utah - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! University of Regensburg - Juergen Fritsch
! Universidad de Madrid - Jose Ortega

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio State University - Dave Drabold

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! iofile3c.f90
! Program Description
! ===========================================================================
!       This subroutine takes a given root and suffix with three index 
! numbers.  This information is then combined to give an input or output 
! filename.  This file is then opened according to the input unit device 
! given.
!
! ===========================================================================
! Code rewritten by:
! James P. Lewis
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! (801) 585-1078 (office)      email: lewis@hec.utah.edu
! (801) 581-4353 (fax)         Web Site: http://
 
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine iofile3c (root, suffix, ith, isorp, index1, index2,       &
     &                       index3, iunit, filename, skip)
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: index1
        integer, intent (in) :: index2
        integer, intent (in) :: index3
        integer, intent (in) :: isorp
        integer, intent (in) :: ith
        integer, intent (in) :: iunit
 
        character (len=12), intent (in) :: root
        character (len=3), intent (in) :: suffix
 
! Output
        character (len=40), intent (out) :: filename

        logical, intent (out) :: skip
 
! Local Parameters and Data Declaration
! ===========================================================================
        integer, parameter :: nplace = 2
 
! Local Variable Declaration and Description
! ===========================================================================
        integer i, j
        integer lr
        integer ls
 
        character (len=80) form
 
! Procedure
! ===========================================================================
! Find true length of root and suffix
        lr = len(root)
        do i = lr, 1, -1
         if (root(i:i) .ne. ' ') exit 
        end do
        lr = i
        ls = len(suffix)
        do i = ls, 1, -1
         if (suffix(i:i) .ne. ' ') exit 
        end do
        ls = i
        if (lr + ls + 5*(nplace + 1) .gt. len(filename)) then
         write (*,'('' filename too short for input in iofile'')')
         write (*,*) root, suffix, index1, index2, index3, iunit
         stop ' error in iofile.f90 '
        end if
 
! Form file name
! ***************************************************************************
! Fill in root here
        do i = 1, lr
         filename(i:i) = root(i:i)
        end do
 
        if  (index1 .lt. 0 .or. index2 .lt. 0 .or. index3 .lt. 0 .or.        &
     &       index1 .gt. 10**nplace .or. index1 .gt. 10**nplace .or.         &
     &       index1 .gt. 10**nplace) then
         write (*,'('' index out of range in iofile'')')
         write (*,*) root, suffix, index1, index2, index3, iunit
         stop ' error in iofile.f90 '
        end if
 
! Write the format to form this will be (i3.3) for nplace = 3
        write (form,'(''(i'',i1,''.'',i1,'')'')') nplace, nplace
 
! Write out the index numbers with a '.' between each number
        filename(lr+1:lr+1) = '_'
        write (filename(lr+2:lr+nplace+1),form) ith
        filename(lr+nplace+2:lr+nplace+2) = '_'
        write (filename(lr+nplace+3:lr+2*nplace+2),form) isorp
        filename(lr+2*nplace+3:lr+2*nplace+3) = '.'
        write (filename(lr+2*nplace+4:lr+3*nplace+3),form) index1
        filename(lr+3*nplace+4:lr+3*nplace+4) = '.'
        write (filename(lr+3*nplace+5:lr+4*nplace+4),form) index2
        filename(lr+4*nplace+5:lr+4*nplace+5) = '.'
        write (filename(lr+4*nplace+6:lr+5*nplace+5),form) index3
        filename(lr+5*nplace+6:lr+5*nplace+6) = '.'
 
! Fill in the suffix here
        do i = 1, ls
         j = lr+5*nplace + 6 + i
         filename(j:j) = suffix(i:i)
        end do
 
! Fill rest with spaces in case there is garbage there now
        do i = lr+5*nplace+7 + ls, len(filename)
         filename(i:i) = ' '
        end do
 
! Filename is formed, open the file
        skip = .false.
        inquire (file = filename, exist = skip)
        if (.not. skip) open (unit = iunit, file = filename, status = 'unknown')
 
! Format Statements
! ===========================================================================
 
        return
        end