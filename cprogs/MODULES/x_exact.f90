        module x_exact
         use precision

         integer nnrhop

         real*8 drhop

         real*8, dimension (:), allocatable :: rpoint
         real*8, dimension (:,:,:), allocatable :: rprime
         real*8, dimension (:,:,:,:), allocatable :: rprime_spline
        end module
