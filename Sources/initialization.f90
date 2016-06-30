!==============================================================================================      
      subroutine startup
      !=================
!     
      use main
!
      integer :: imc,jmc
      real*8 :: Xcar,Ycar,flag,xx,yy
      real*8 :: e,TT,KE
      real*8 :: temp
      integer :: IO_STATUS
      
      namelist /step_/ nmax,StepSub,dt_physics,dt_print,StepSub,is_restart,iprint,steady1
      nmax=10000     !         nmax:  maximum advancing steps;
      StepSub=20       ! LUSGS setup
      dt_physics =0.001
      dt_print =0.02
      is_restart=0  !         is_restart: start from very beginning(0)                      or start from a previous computation(1);
      iprint=1000    !         iprint: output every iprint steps
      steady1=.false.       !         steady1: logical variable, 
                            !.false. the flow is unsteady
      icgl=0        !         icgl=1: gloal time step; =0 local time step
      cfl=0.1   ! CFL number
      namelist /solver_/ irsolver,iTimeMarch,iReconstruct
      irsolver=0 ! Roe
      iTimeMarch=0 ! LUSGS
      iReconstruct=0 ! MUSCL
      
       
      
      open (4,file='D:\\CFD_DM\\DM\\Input\\euler.in',mode='READ',status='unknown')
!
!     read in main integer/logical control variables
!     ----------------------------------------------
!
      read(4,NML=step_,IOSTAT=IO_STATUS) 
      read(4,NML=solver_,IOSTAT=IO_STATUS)  
       close(4)
      
      !read(4,*) cp,Gamma,Pr_L
!         cp: specfic heat at constant pressure
!         Gamma: ratio of specific heat
!         Pr_L:Prandt.al number
      R_air=287.0
      Gamma=1.4
      Pr_L=0.72
      Pr_T=0.9
      


     

!     read inlet/outlet parameters: uniform inlet, need to be modified.
!     -----------------------------------------------------------------
!     read(4,*) Vx_inlet,Vy_inlet,p_inlet,T_inlet,p_out
!============================
      ! boundary conditions:
      
      !Vx_inlet=1000.0
      !Vy_inlet=0.0
      !p_inlet=101325.0
      !T_inlet=300.0
      !KT_inlet=9.0e-9 !
      !OmegaT_inlet=1.0e-6 !
      !rho_inlet=p_inlet/T_inlet/R_air
      !Mu_inlet=1.458d-6*abs(T_inlet)**1.5/(T_inlet+110.4)
      !P_out=0.0
      !===============================

!===========end of the input file bl.in=========================
!
!     read grids 
!     ----------
!
        open (5,file='D:\\CFD_DM\\DM\\Input\\grid.in',status='unknown')
        read(5,*) imc,jmc
        ib=3
        jb=3
        im=imc+ib
        jm=jmc+jb
        do i=ib,im
        do j=jb,jm
        
            read(5,*) x(i,j),y(i,j)
        end do
        end do
        
      close(5)
      
      
      
      i=ib
      do j=jb,jm
      x(i-1,j)=2.0*x(i,j)-x(i+1,j)
      y(i-1,j)=y(i,j)
      end do
      
      i=im
      do j=jb,jm
      x(i+1,j)=2.0*x(i,j)-x(i-1,j)
      y(i+1,j)=y(i,j)
      end do
      
      j=jb
      do i=ib-1,im+1
      x(i,j-1)=x(i,j)
      y(i,j-1)=2.0*y(i,j)-y(i,j+1)
      end do
      
      j=jm
      do i=ib-1,im+1
      x(i,j+1)=x(i,j)
      y(i,j+1)=2.0*y(i,j)-y(i,j-1)
      end do
!
!     end of data input section
!     -------------------------
!

      open (30,file='grid.dat')
      write(30,*) 'title="contour"'
      write(30,*) 'variables="x","y"'
     
          write(30,*) 'zone i=',im-ib+1,' j=',jm-jb+1&
                    &,' f=point'
       
          do j=jb,jm
          do i=ib,im
             xcar=x(i,j)
             ycar=y(i,j) 
             write(30,*) xcar,ycar
          end do
          end do
        
       close(30)
            
!============================================================================      
!
!      non-dimensionalization
!   now all reference variables are set to 1, NO non-dimensionalization
      
       L_ref=1. !x(im,jb)-x(ib,jb)
       Rho_ref=1. !p_Inlet/(T_Inlet*(ga-1.)*cp)
       V_ref=1. !sqrt(p_Inlet/Rho_ref)
       T_ref=1. !T_Inlet 
       !basec=1. !sqrt((ga-1.)*cp*T_ref)
       Ma_ref=1. !sqrt(Vx_Inlet*Vx_Inlet+Vy_Inlet*Vy_Inlet)/basec
       
       
       write(*,*) 'baseparameters', L_ref,Rho_ref,V_ref,T_ref,Ma_ref
       
!      non-dim parameters

       cp=cp*T_ref/(V_ref*V_ref)
       
!      non-dim inlet/outlet parameters 
       
       Vx_Inlet=Vx_Inlet/V_ref
       Vy_Inlet=Vy_Inlet/V_ref
       
       p_Inlet=p_Inlet/(Rho_ref*V_ref*V_ref)
       T_Inlet=T_Inlet/T_ref
       p_Out=p_Out/(Rho_ref*V_ref*V_ref)
       
!      non-dim coordinates       
       
        do i=ib,im
        do j=jb,jm
           x(i,j)=x(i,j)/L_ref
           y(i,j)=y(i,j)/L_ref
        end do
        end do
!
!
!     set various constants
!     ---------------------
!
      Gamma1=Gamma-1
      Gamma2=(Gamma-1.0)/Gamma
      Gamma3=1.0/Gamma2
      cp=R_air*Gamma3
      

!
!     compute aera vectors
!     --------------------------------------
!

      call area
      call volume_static    
       
!     ==========
!
!
!
!     initialization
!     --------------
!
      if(is_restart.eq.0) n=0

!================================================
!     initial condition
!     -----------------
!
!!
!      do i=1,iq
!      do j=1,jq
!
!             ! xx=0.25*(x(i,j)+x(i,j+1)+x(i+1,j)+x(i+1,j+1))
!              !yy=0.25*(y(i,j)+y(i,j+1)+y(i+1,j)+y(i+1,j+1))
!              ! initialize the flow field with 
!                   rho(i,j)=rho_inlet
!                   vx(i,j)=Vx_inlet
!                   vy(i,j)=Vy_inlet
!                    p(i,j)=p_inlet
!                    
!                    T(i,j)=T_inlet
!                ! initialize turbulence variables
!                    KT(i,j)=KT_inlet*Rho_inlet
!                    OmegaT(i,j)=OmegaT_inlet*Rho_inlet
!                    Mu_L(i,j)=Mu_inlet
!                    Mu_T(i,j)=rho_inlet*KT_inlet/OmegaT_inlet
!                    Mu_E(i,j)=Mu_T(i,j)+Mu_L(i,j)
!      end do
!      end do
      do i=1,iq-1
      do j=1,jq-1
              xx=0.25*(x(i,j)+x(i,j+1)+x(i+1,j)+x(i+1,j+1))
              yy=0.25*(y(i,j)+y(i,j+1)+y(i+1,j)+y(i+1,j+1))
              flag=1.73205*(xx-0.1666667)-yy
              if(flag.gt.0) then
                   rho(i,j)=1.4
                   vx(i,j)=0.
                   vy(i,j)=0.
                    p(i,j)=1.0
              else
                    rho(i,j)=8.0
                    vx(i,j)=7.1447
                    vy(i,j)=-4.125
                     p(i,j)=116.5
              endif
      end do
      end do
      
      !initialize the last and 2nd last step solutions:
      do  i=1,iq
      do  j=1,jq
      
      e=.5*(vx(i,j)*vx(i,j)+vy(i,j)*vy(i,j))
      rho_Et(i,j)=p(i,j)/(Gamma-1)+(e)*rho(i,j)
      Ht(i,j)=Gamma*(rho_Et(i,j)/rho(i,j)-e)+e
      rho_vx(i,j)=rho(i,j)*vx(i,j)
      rho_vy(i,j)=rho(i,j)*vy(i,j)
      T(i,j)=p(i,j)/(rho(i,j)*R_air)
      Mu_L(i,j)=1.458d-6*abs(T(i,j))**1.5/(T(i,j)+110.4)
      Mu_E(i,j)=Mu_L(i,j)
      !!!!!! 
      KT(i,j)=1.0
      OmegaT(i,j)=1.0
      
      Rho_m1(i,j)=rho(i,j)
      Rho_Et_m1 (i,j)=rho_Et(i,j)
      Rho_vx_m1(i,j)=rho_vx(i,j)
      Rho_vy_m1(i,j)=rho_vy(i,j)
      !KT_m1(i,j)=KT(i,j)
      !OmegaT_m1(i,j)=OmegaT(i,j)
      
      Rho_m2(i,j)=rho(i,j)
      Rho_Et_m2(i,j)=rho_Et(i,j)
      Rho_vx_m2(i,j)=rho_vx(i,j)
      Rho_vy_m2(i,j)=rho_vy(i,j)
      !KT_m2(i,j)=KT(i,j)
      !OmegaT_m2(i,j)=OmegaT(i,j)
      end do
      end do
           
     
!
!     restart from backup intermediate results
!     --------------------------
!
      if(is_restart.eq.0) goto 10000
!
!      read:
!      velocities, static pressure, density,viscosity
!      ----------------------------------------------
!
      write(*,*) 'read  backuped intermediate results'
      open(5,file='bl.rd', form='unformatted')  !
      read(5) n,ttime
      
      do  i=ib,im
      do  j=jb,jm
             read(5) vx(i,j),vy(i,j)&
     &               ,p(i,j),rho(i,j)&
     &               ,Mu_L(i,j),Mu_T(i,j)&
     &               ,KT(i,j),OmegaT(i,j)&
     &               ,Rho_m1(i,j),Rho_Et_m1(i,j)&
     &               ,Rho_vx_m1(i,j),Rho_vy_m1(i,j)
     
      end do
      end do
!
      close(5)
      
      do i=ib,im-1
      do j=jb,jm-1
      TT=T_Inlet*T_ref
      Mu_L(i,j)=1.458*abs(TT)**1.5/(TT+110.4)*1.0d-6/(Rho_ref*V_ref*L_ref)
      Mu_E(i,j)=Mu_L(i,j)
      end do
      end do
      
      write(*,*) 'read  backuped intermediate results finish'
!
!      other variables
!      ---------------
!
      do  i=ib,im-1
      do  j=jb,jm-1
      
      KE=.5*(vx(i,j)*vx(i,j)+vy(i,j)*vy(i,j))
      rho_Et(i,j)=p(i,j)/Gamma1+KE*rho(i,j)
      Ht(i,j)=Gamma*(rho_Et(i,j)/rho(i,j)-KE)+KE
      rho_vx(i,j)=rho(i,j)*vx(i,j)
      rho_vy(i,j)=rho(i,j)*vy(i,j)
      T(i,j)=p(i,j)/(rho(i,j)*R_air)
      
      end do
      end do

10000 continue
!
      return
      end
