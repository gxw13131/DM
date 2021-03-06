!==============================================================================================      
      subroutine startup
      !=================
!     
      use main
!
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
      read(4,*) cp,Gamma,Pr_L
!         cp: specfic heat at constant pressure
!         Gamma: ratio of specific heat
!         Pr_L:Prandt.al number
      Pr_E=Pr_L
      Pr_T=Pr_L
      
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
       basec=1. !sqrt((ga-1.)*cp*T_ref)
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
      cv=cp/Gamma
      R_air=cp-cv

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
      do i=ib,im-1
      do j=jb,jm-1
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
      
      
      do  i=ib,im-1
      do  j=jb,jm-1
      
      e=.5*(vx(i,j)*vx(i,j)+vy(i,j)*vy(i,j))
      rho_Et(i,j)=p(i,j)/(Gamma-1)+(e)*rho(i,j)
      Ht(i,j)=Gamma*(rho_Et(i,j)/rho(i,j)-e)+e
      ros=rho(i,j)
      rho_vx(i,j)=ros*vx(i,j)
      rho_vy(i,j)=ros*vy(i,j)
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
      temini=T_Inlet*T_ref
      Mu_L(i,j)=1.458*abs(temini)**1.5/(temini+110.4)*1.0d-6/&
         (Rho_ref*V_ref*L_ref)
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
      
      vsqh=.5*(vx(i,j)*vx(i,j)+vy(i,j)*vy(i,j))
      rho_Et(i,j)=p(i,j)/Gamma1+vsqh*rho(i,j)
      Ht(i,j)=Gamma*(rho_Et(i,j)/rho(i,j)-vsqh)+vsqh
      ros=rho(i,j)
      rho_vx(i,j)=ros*vx(i,j)
      rho_vy(i,j)=ros*vy(i,j)
      T(i,j)=p(i,j)/(rho(i,j)*R_air)
      
      end do
      end do

10000 continue
!
      return
      end
