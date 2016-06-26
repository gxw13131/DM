!=====================================================================================      
      subroutine solver
!     *****************
!
      use main
      real*8 :: time1,time2
!*****************************************************************************
      CALL CPU_TIME ( time1 )
!     compute the area-vectors and interpolation coordinates  
        
!*****************************************************************************
!     set ghost-cell at boundary for various boundary conditions
      call bc_ghostcell_value

!     interpolate the primitive variables to cell interfaces
      call interpolation_to_interface

!     MUSCL interpolation at cell interface except boundaries
      call muscl_interpolation 
!     set muscl_interpolated value at the cell interfaces on the boundary
      call Reconstruct
!!!!  TURBULENCE interpolate the interface variables
      call SST_Reconstruct_MUSCL    !T
!      fix BC flux to guarantee BC
      call BC_FIX
!     evaluate the inviscid fluxes 
      call RiemannSolver
!!!! TURBULENCE convective item
      call SST_Flux_Lax      !T
      
!       compute the viscous flux
      call viscous_Flux
!!!! TURBULENCE compute SS,Vorticity,dKdO 
      call Derivative2rd    !T
      ! Diffusion+Production+Destruction item
      call SST_RHS          !T

      call pSolver
!!!! TURBULENCE update k and Omega
 !     call SST_LUSGS        !T
     
!     update all variables
      call update_variables
!     ===================== 
      
      CALL CPU_TIME ( time2 )
      
      !write(*,*) 'time spent'
      !write(*,*) time2-time1
    !  pause
      return
      end
      