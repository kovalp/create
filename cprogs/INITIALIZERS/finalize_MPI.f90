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

! finalize_MPI.f
! Program Description
! ===========================================================================
!       This file finalizes all the MPI processors. 
!
! ===========================================================================
! Code rewritten by:
! James P. Lewis
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
      	subroutine finalize_MPI 
      	implicit none

        include 'mpif.h'

! Argument Declaration and Description
! ===========================================================================

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer ierr_mpi

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
!  MPI Broadcast signature
        call MPI_Finalize (ierr_mpi)

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================

	return
        end
