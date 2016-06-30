!==================================================================!
!==================================================================!
!
!
      program Euler
     !======================
     

     use main

     integer :: IS_END
     real*8 :: TIME1,TIME2
     integer :: N_loop,II
!
      open (2,file='Euler.his',status='unknown')
!     bl.his:file that store the computation history
      
!
!     initialization
!     --------------
      write(*,*) '1'
!
      call startup
      ! compute the wall distance
      call DistanceWall
      write(*,*) '2'
!
!     start of main iteration loop
!----------------------------
      !iTimeMarch=0 ! LUSGS is the default time marching method
      select case (iTimeMarch)
      case (1)
       pSolver=>RK2
       N_loop=2
       NameTime='RK2'
      case (2)
       pSolver=>RK4
       N_loop=4
       NameTime='RK4'
      case default
       pSolver=>LUSGS
       N_loop=StepSub
       NameTime='LUSGS'
      end select
!--------------------------------------------
!       select the Reconstruction method
!       1: MUSCL; 2: WENO3; 3: WENO5
!       default is MUSCL
      !iReconstruct=0
      select case (iReconstruct)
      case (1)
        Reconstruct => WENO3_interpolation
        NameRec='WENO3'
      case (2)
        Reconstruct => WENO5_interpolation
        NameRec='WENO5'
      case default
        Reconstruct => muscl_interpolation
        NameRec='MUSCL'
      end select
!---------------------------------------------
!       select the Riemann solver
!       1: Lax ; 2: AUSMPW+ ; others: Roe
!       default is Roe 
      
      select case (irsolver)
      case (1)
        RiemannSolver => inviscid_fluxes_lax
        NameRiemann='Lax'
      case (2)
        RiemannSolver => inviscid_fluxes_AUSMPWplus
        NameRiemann='AUSMPWplus'
      case default
        RiemannSolver => inviscid_fluxes_roe
        NameRiemann='Roe'
      end select
!================================================================================== 
!==================================================================================
      is_End=0
!     if is_End==1, the computation stops
      
      istart=0
!     istart:store the advancing time steps for every new start up
!
!     start a new time step
!     ---------------------
!=====================================================================
! begin the main loop

      Nprint=0
      do while(is_End.ne.1)
100   is_print=0
!     if is_print==1 output the results

      
!     n: accumulated time step of the computation
      istart=istart+1
      
      if(mod(n,10).eq.0) then
          write(*,*) 'now, computing the ',n,'th step'
          write(2,*) 'now, computing the ',n,'th step'
      end if
      if((n.eq.nmax).or.(Nprint.eq.10)) is_End=1
      
      is_print=0
      !if(mod(n,iprint).eq.0) is_print=1
      if(ttime.gt.dt_print*Nprint) then
      is_print=1
      !read(*,*)
      Nprint=Nprint+1 ! Nprint is initialized to 0, initial flowfield will be printed
      end if
!     compute the minimum time step to guarantee CFL condition
      call radius
      call time_step

!     solve  equations
!     ------------------
      do ii=1,N_loop !RK-2 needs 2 loops, RK-4 needs 4 loops, LUSGS needs StepSub loops
       call solver
       n=n+1
      end do
      if(.not.steady1)then
          if(iTimeMarch.eq.0) then
          ttime=ttime+dt_physics
          else
           ttime=ttime+step(ib,jb)
           end if
      end if
!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! save solution
      !     store variables of previous time steps
      !     ---------------------------------------
      ! print the rms of equations
      call rms
      
      call saveStep
      
      call output
      if((.not.steady1).and.(is_print.eq.1)) call outputuns
            
      end do !
      Nprint=100
      call outputuns
!=========================================================
!     end of main loop
!========================================================
      close (2)
      write(*,*) 'Normal Termination' 
      write(*,*) time2-time1
      
      
      end program
!
!=================================================
      subroutine saveStep()
      use main
      implicit none
        do i=ib,im
        do j=jb,jm 
!       save last step solution
            Rho_m2(i,j)=Rho_m1(i,j)
            Rho_vx_m2(i,j)=rho_vx_m1(i,j)
            Rho_vy_m2(i,j)=rho_vy_m1(i,j)
            Rho_Et_m2(i,j)=rho_Et_m1(i,j)
!       save current step solution      
            Rho_m1(i,j)=rho(i,j)
            Rho_vx_m1(i,j)=rho_vx(i,j)
            Rho_vy_m1(i,j)=rho_vy(i,j)
            Rho_Et_m1(i,j)=rho_Et(i,j)
!
      end do
      end do
            
      end subroutine

!   compute the RMS and print to screen and file      
      subroutine rms
      use main
      implicit none
      REAL*8 ::  rms1,rms2,rms3,rms4
      REAL*8 ::  r1,r2,r3,r4
     integer :: imax1,jmax1, &
     &            imax2,jmax2, &
     &            imax3,jmax3, &
     &            imax4,jmax4
      real*8 :: Tinv
       rms1=0.
       rms2=0.
       rms3=0.
       rms4=0.
      
      
      do i=ib,im-1
      do j=jb,jm-1
      
      tinv=1./step(i,j)
      
      r1=abs(rho(i,j)-rho_m1(i,j))*tinv
      r2=abs(rho_vx(i,j)-Rho_vx_m1(i,j))*tinv
      r3=abs(rho_vy(i,j)-Rho_vy_m1(i,j))*tinv
      r4=abs(rho_Et(i,j)-Rho_Et_m1(i,j))*tinv
      
      if(r1>rms1) then
        rms1=r1
        imax1=i
        jmax1=j
      end if 
      
      if(r2>rms2) then
        rms2=r2
        imax2=i
        jmax2=j
      end if 
      
      if(r3>rms3) then
        rms3=r3
        imax3=i
        jmax3=j
      end if
      
      if(r4>rms4) then
        rms4=r4
        imax4=i
        jmax4=j
      end if  
      
      end do
      end do
      
      if(mod(n,10).eq.0) then
          write(*,*) 'time', ttime*L_ref/V_ref
          write(*,233) rms1,imax1,jmax1
          write(*,233) rms2,imax2,jmax2
          write(*,233) rms3,imax3,jmax3
          write(*,233) rms4,imax4,jmax4
        
          write(2,*) rms1,imax1,jmax1
          write(2,*) rms2,imax2,jmax2
          write(2,*) rms3,imax3,jmax3
          write(2,*) rms4,imax4,jmax4
          endif
233     FORMAT(5X,"RMS=",1X,E10.3,1X,"at(",I3,",",I3,")")
        end subroutine
          