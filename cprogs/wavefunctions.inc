! These common blocks store information about the wavefunctions
! You should not need to change any of this. But here it all is anyway.

        integer npoints (nsh_max, nspec_max)
        integer npoints_rho (nspec_max)
 
        integer nnrho
        integer nnz
 
        real*8 drho
        real*8 dz
        real*8 rhomax
        real*8 rhomin
        real*8 zmax
        real*8 zmin
 
        real*8 drr (nsh_max, nspec_max)               ! spacing between points
        real*8 drr_rho (nspec_max)                    ! spacing between points
        real*8 rr (wfmax_points, nsh_max, nspec_max)
        real*8 rr_rho (wfmax_points, nspec_max)
        real*8 rrc (nsh_max, nspec_max)               ! cutoff points
        real*8 rrc_rho (nspec_max)                    ! maximum cutoff points
        real*8 psi (wfmax_points, nsh_max, nspec_max)
        real*8 psi_spline (4,wfmax_points, nsh_max, nspec_max)
        real*8 rho2c (nrho_points, nz_points)
        real*8 rhop2c (nrho_points, nz_points)
        real*8 rhopp2c (nrho_points, nz_points)
        real*8 rhopz2c (nrho_points, nz_points)
        real*8 rhoz2c (nrho_points, nz_points)
        real*8 rhozz2c (nrho_points, nz_points)
 
        common /density_limits_rho/ drho, rhomin, rhomax, nnrho
        common /density_limits_z/ dz, zmin, zmax, nnz
        common /densities/ rho2c, rhop2c, rhopp2c, rhopz2c,
     1                            rhoz2c, rhozz2c
        common /densparam/ drr_rho, rr_rho, rrc_rho, npoints_rho
        common /wavefunctions/ psi, psi_spline
        common /waveparam/ drr, rr, rrc, npoints
 
!$OMP   THREADPRIVATE(/density_limits_rho/,/density_limits_z/,
!$OMP1                /densities/)
 
