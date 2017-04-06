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
 
!
! RESTRICTED RIGHTS LEGEND
!
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.
 
! Program Description
! ===========================================================================
! Initializes the MPI space
! ===========================================================================
 
! Code written by:
! Kurt R. Glaesemann
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E. RM Dock
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine init_MPI(iammaster, iammpi, my_proc, nprocs)
        implicit none

! Argument Declaration and Description
! ======================================================================
! Output
        integer my_proc 
        integer nprocs

        logical iammaster
        logical iammpi

! Local Parameters and Data Declaration
! ======================================================================

! Local Variable Declaration and Description
! ======================================================================
! We are NOT doing MPI
        iammaster = .true.
        iammpi = .false.
        my_proc = 0
        nprocs = 1

! Format Statements
! ======================================================================

        return
        end
