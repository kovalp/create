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
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine init_MPI (iammaster, iammpi, my_proc, nproc)
        implicit none
 
        include 'mpif.h'
 
! Argument Declaration and Description
! ===========================================================================
! Output
        integer, intent(out) :: my_proc
        integer, intent(out) :: nproc
 
        logical, intent(out) :: iammaster
        logical, intent(out) :: iammpi
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ierr_mpi
 
! Procedure
! ===========================================================================
        call MPI_Init (ierr_mpi)
        if (ierr_mpi .ne. 0) then
          write(*,*) 'mpi error MPI_Init'
          stop
        end if
        call MPI_Comm_rank (MPI_COMM_WORLD, my_proc, ierr_mpi)
        if (ierr_mpi .ne. 0) then
          write(*,*) 'mpi error MPI_Comm_rank'
          stop
        end if
        call MPI_Comm_size (MPI_COMM_WORLD, nproc, ierr_mpi)
        if (ierr_mpi .ne. 0) then
          write(*,*) 'mpi error MPI_Comm_size'
          stop
        end if
        iammaster = .false.
        if (my_proc .eq. 0) iammaster = .true.
        iammpi = .true.
 
! Format Statements
! ===========================================================================
 
        return
        end

! Program Declaration
! ===========================================================================
        subroutine signature_MPI (signature)
        implicit none

        include 'mpif.h'

! Argument Declaration and Description
! ===========================================================================
        character*30 signature

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer ierr_mpi

! Procedure
! ===========================================================================
        call MPI_BCAST (signature, 30, MPI_CHARACTER, 0, MPI_COMM_WORLD,     &
     &                  ierr_mpi)

! Format Statements
! ===========================================================================

        return
        end

! Program Declaration
! ===========================================================================
        subroutine Finalize_MPI()
        implicit none
        include 'mpif.h'

! Argument Declaration and Description
! ===========================================================================

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer ierr_mpi

! Procedure
! ===========================================================================
        call MPI_barrier(MPI_COMM_WORLD, ierr_mpi)
        call MPI_Finalize(ierr_mpi)
! Format Statements
! ===========================================================================

        return
        end

