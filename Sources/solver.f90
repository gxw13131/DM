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

!     evaluate the inviscid fluxes 
      call RiemannSolver
      
!       compute the viscous flux
      call viscous_Flux

!*****************************************************************************
!     Runge-Kutta time stepping
   !   call Runge_Kutta
  !    call LUSGS
      call pSolver

!    Turbulence model SST      
      call SST_Reconstruct_MUSCL   
      call SST_Flux_Lax
      call SST_LUSGS
!*****************************************************************************
!     update variables
      call update_variables
!     ===================== 
      
      CALL CPU_TIME ( time2 )
      
      !write(*,*) 'time spent'
      !write(*,*) time2-time1
    !  pause
      return
      end
      