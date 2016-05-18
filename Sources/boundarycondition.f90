!=========================================================================================
     
      subroutine bc_ghostcell_value
     !===================================
      use main
      real*8 :: xhalf
      ! set ghost-cell at boundary for various boundary conditions
      
      
      !I-boundary
      
      !inlet
      
      ! 方法：在i=ib-1 这个虚拟网格上给定进口条件 
      
      !暂时用均匀进口条件
      
      i=ib-1
     
      do j=jb,jm-1
  
      rho(i,j)=rho_inlet
      vx(i,j)=Vx_inlet
      vy(i,j)=Vy_inlet
      p(i,j)=p_inlet
      Mu_L(i,j)=Mu_inlet
      Mu_E(i,j)=Mu_L(i,j)
       !    rho(i,j)=max(2.*rho(i+1,j)-rho(i+2,j),1.0e-10)
       !    vx(i,j)=2.*vx(i+1,j)-vx(i+2,j)
       !    vy(i,j)=2.*vy(i+1,j)-vy(i+2,j)
       !     p(i,j)=max(2.*p(i+1,j)-p(i+2,j),1.0e-10)
       !  Mu_E(i,j)=2.*Mu_E(i+1,j)-Mu_E(i+2,j)
       !  Mu_L(i,j)=2.*Mu_L(i+1,j)-Mu_L(i+2,j)
       L_Cell_x(i,j)=L_Cell_x(i+1,j)
     
      end do
     
  
     !outlet，外推
      i=im
     
      do j=jb,jm-1
     
      p (i,j)=max(2.*p(i-1,j)-p(i-2,j),1.0e-10) 
      vx(i,j)=2.*vx(i-1,j)-vx(i-2,j)
      vy(i,j)=2.*vy(i-1,j)-vy(i-2,j)
      rho(i,j)=max(2.*rho(i-1,j)-rho(i-2,j),1.0e-10)
    Mu_E(i,j)=2.*Mu_E(i-1,j)-Mu_E(i-2,j)
    Mu_L(i,j)=2.*Mu_L(i-1,j)-Mu_L(i-2,j)
  L_Cell_x(i,j)=L_Cell_x(i-1,j)
  
      end do
      
     
      
      
      !j-boundary
      
      ! solid wall boundary
      
      j=jb
      
      do i=ib,im-1
      
      xhalf=0.5*(x(i,j)+x(i+1,j))
      if(xhalf>=1.) then  !x<1,均匀来流条件，否则固壁条件
      vx(i,j-1)=vx(i,j)
      vy(i,j-1)=-vy(i,j)
      p(i,j-1)=p(i,j)
      rho(i,j-1)=rho(i,j)
      Mu_E(i,j-1)=Mu_E(i,j)
      Mu_L(i,j-1)=Mu_L(i,j)
      L_Cell_y(i,j-1)=L_Cell_y(i,j)
      else
      vx(i,j-1)=2*vx(i,j)-vx(i,j+1)
      vy(i,j-1)=2*vy(i,j)-vy(i,j+1)
      p (i,j-1)=2*p (i,j)-p (i,j+1)
      rho(i,j-1)=2*rho(i,j)-rho(i,j+1)
      Mu_E(i,j-1)=2*Mu_E(i,j)-Mu_E(i,j+1)
      Mu_L(i,j-1)=2*Mu_L(i,j)-Mu_L(i,j+1)
      L_Cell_y(i,j-1)=L_Cell_y(i,j)
      end if
      end do
      
      j=jm
      
      do i=ib,im-1
      
      vx(i,j+1)=vx(i,j)
      vy(i,j+1)=-vy(i,j)
      p(i,j+1)=p(i,j)
      rho(i,j+1)=rho(i,j)
      Mu_E(i,j+1)=Mu_E(i,j)
      Mu_L(i,j+1)=Mu_L(i,j)
      L_Cell_y(i,j+1)=L_Cell_y(i,j)
      
      end do
       
      return
      end

!================================================================================== 

      subroutine bc_muscl_interpolation
      !================================
      use main
      
      ! set muscl_interpolated value at the cell interfaces on the boundary
      
      
      !I-boundary
      
      !inlet
      
      ! 
      i=ib
      
      do j=jb,jm-1
     
          
          uil(i,j)=Vx_inlet
          vil(i,j)=Vy_inlet
          pil(i,j)=p_inlet
          ril(i,j)=rho_inlet
      
      end do
      
      !outlet 
      i=im
      
    !  do j=jb,jm-1
       
    !  end do
     
      !j-boundary
      
      ! solid wall boundary
      
      j=jb
      
      do i=ib,im-1
      	
    xhalf=0.5*(x(i,j)+x(i+1,j))
      if(xhalf>=1.) then  !x<1,均匀来流条件，否则固壁条件  
      ujl(i,j)=0.0
      vjl(i,j)=0.0
      pjl(i,j)=0.5*(p(i,j)+p(i,j-1))
      rjl(i,j)=0.5*(rho(i,j)+rho(i,j-1))
      else
      ujl(i,j)=0.5*(vx(i,j)+vx(i,j-1))
      vjl(i,j)=0.5*(vy(i,j)+vy(i,j))
      pjl(i,j)=0.5*(p(i,j)+p(i,j-1))
      rjl(i,j)=0.5*(rho(i,j)+rho(i,j-1))
      end if
      
      end do
      
      j=jm
      
! upper boundary 时间相关的边界条件（先判断激波位置，再决定用激波前或者激波后的参数）
      
      do i=ib,im-1
      
         
            
      
      ujr(i,j)=0.5*(vx(i,j)+vx(i,j+1))
      vjr(i,j)=0.5*(vy(i,j)+vy(i,j+1))
      pjr(i,j)=0.5*(p(i,j)+p(i,j+1))
      rjr(i,j)=0.5*(rho(i,j)+rho(i,j+1))
 
      
      end do

       
      return
      end
     