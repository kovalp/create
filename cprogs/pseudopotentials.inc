! These common blocks store information about the (non)-neutral potential
! You should not need to change any of this. But here it all is anyway.

        integer npoints_pp (0:nsh_max, nspec_max)
 
        real*8 cl_pp (nsh_max, nspec_max)
        real*8 drr_pp (nsh_max, nspec_max)           ! spacing between points
        real*8 rr_pp (max_points_pp, nsh_max, nspec_max)
        real*8 rrc_pp (nsh_max, nspec_max)
        real*8 vpp (max_points_pp, nsh_max, nspec_max)

        real*8 vpp_spline (4,max_points_pp, nsh_max, nspec_max)
        real*8 rcPP (nspec_max)                                                                                             
        common /pp/ rr_pp, vpp, vpp_spline
        common /pp_param/ cl_pp, drr_pp, rrc_pp, npoints_pp, rcPP

