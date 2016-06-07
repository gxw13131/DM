!==============================================================================================      
      subroutine startup
      !=================
!     
      use main
!
      integer :: imc,jmc
      real*8 :: Xcar,Ycar
      real*8 :: e,TT,KE
      real*8 :: temp
      open (4,file='D:\\CFD_DM\\DM\\Input\\euler.in',status='unknown')
!
!     read in main integer/logical control variables
!     ----------------------------------------------
!
      read(4,*) nmax,is_restart,iprint 
!         nmax:  maximum advancing steps;
!         is_restart: start from very beginning(0) or start from a previous computation(1);
!         iprint: output every iprint steps
!

      read(4,*) steady1
!         steady1: logical variable, ==.true. compute steady flow; ==.false. the flow is unsteady

      read(4,*) icgl, irsolver  
!         icgl=1: gloal time step; =0 local time step
!         irsolver=0 Roe; =1 Lax -riemann solver
!     read constants and CFL number 
!     ------------------------------
!
      !read(4,*) cp,Gamma,Pr_L
!         cp: specfic heat at constant pressure
!         Gamma: ratio of specific heat
!         Pr_L:Prandt.al number
      !R_air=287.0
      !Gamma=1.4
      !Pr_L=0.72
      !Pr_E=Pr_L
      !Pr_T=Pr_L
      
      read(4,*) cfl
!         cfl: cfl number

      pi=3.1415926535897932384626

!     read inlet/outlet parameters: uniform inlet, need to be modified.
!     -----------------------------------------------------------------
!
!      read(4,*) Vx_Inlet,Vy_Inlet,p_Inlet,T_Inlet
!      read(4,*) p_Out
!
      close(4)
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
      
      !============================
      ! boundary conditions:
      R_air=287.0
      Gamma=1.4
      Pr_L=0.72
      Pr_E=Pr_L
      Pr_T=Pr_L
      
      Vx_inlet=1000.0
      Vy_inlet=0.0
      p_inlet=101325.0
      T_inlet=288.0
      rho_inlet=p_inlet/T_inlet/R_air
      Mu_inlet=1.458d-6*abs(T_inlet)**1.5/(T_inlet+110.4)
      P_out=0.0
      !===============================
      
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
!
      
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

!
!     initial condition
!     -----------------
!
!
      do i=1,iq
      do j=1,jq

             ! xx=0.25*(x(i,j)+x(i,j+1)+x(i+1,j)+x(i+1,j+1))
              !yy=0.25*(y(i,j)+y(i,j+1)+y(i+1,j)+y(i+1,j+1))
!==============================================================
              ! initilize the flow field with 
                   rho(i,j)=rho_inlet
                   vx(i,j)=Vx_inlet
                   vy(i,j)=Vy_inlet
                    p(i,j)=p_inlet
                    Mu_L(i,j)=Mu_inlet
                    Mu_E(i,j)=Mu_inlet
                    T(i,j)=T_inlet

      end do
      end do
      
      
      do  i=1,iq
      do  j=1,jq
      
      e=.5*(vx(i,j)*vx(i,j)+vy(i,j)*vy(i,j))
      rho_Et(i,j)=p(i,j)/(Gamma-1)+(e)*rho(i,j)
      Ht(i,j)=Gamma*(rho_Et(i,j)/rho(i,j)-e)+e
      rho_vx(i,j)=rho(i,j)*vx(i,j)
      rho_vy(i,j)=rho(i,j)*vy(i,j)
      T(i,j)=p(i,j)/(rho(i,j)*R_air)

      Rho_m1(i,j)=rho(i,j)
      Rho_Et_m1 (i,j)=rho_Et(i,j)
      Rho_vx_m1(i,j)=rho_vx(i,j)
      Rho_vy_m1(i,j)=rho_vy(i,j)
      
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
