!=========================================================================================
     
      subroutine bc_ghostcell_value
     !===================================
      use main
      
      ! set ghost-cell at boundary for various boundary conditions
      
      
      !I-boundary
      
      !inlet
      
      ! ��������i=ib-1 ������������ϸ����������� 
      
      !��ʱ�þ��Ƚ�������
      
      i=ib-1
     
      do j=jb,jm-1
  
           rho(i,j)=max(2.*rho(i+1,j)-rho(i+2,j),1.0e-10)
           vx(i,j)=2.*vx(i+1,j)-vx(i+2,j)
           vy(i,j)=2.*vy(i+1,j)-vy(i+2,j)
            p(i,j)=max(2.*p(i+1,j)-p(i+2,j),1.0e-10)
         Mu_E(i,j)=2.*Mu_E(i+1,j)-Mu_E(i+2,j)
         Mu_L(i,j)=2.*Mu_L(i+1,j)-Mu_L(i+2,j)
       L_Cell_x(i,j)=L_Cell_x(i+1,j)
     
      end do
     
  
     !outlet������
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
      if(xhalf>=1./6.) then  !x<1/6,������������������̱�����
      vx(i,j-1)=vx(i,j)
      vy(i,j-1)=-vy(i,j)
      p (i,j-1)=p (i,j)
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
      
      vx(i,j)=2*vx(i,j-1)-vx(i,j-2)
      vy(i,j)=2*vy(i,j-1)-vy(i,j-2)
      p (i,j)=2*p (i,j-1)-p (i,j-2)
      rho(i,j)=2*rho(i,j-1)-rho(i,j-2)
     Mu_E(i,j)=2*Mu_E(i,j-1)-Mu_E(i,j-2)
      Mu_L(i,j)=2*Mu_L(i,j-1)-Mu_L(i,j-2)
      L_Cell_y(i,j)=L_Cell_y(i,j-1)
      
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
     
          a=8.0
          b=7.1447
          c=-4.125
          d=116.5
          uil(i,j)=b
          vil(i,j)=c
          pil(i,j)=d
          ril(i,j)=a
      
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
      if(xhalf>=1./6.) then
     
      ujr(i,j)=ujl(i,j)
      vjr(i,j)=vjl(i,j)
      pjr(i,j)=pjl(i,j)
      rjr(i,j)=rjl(i,j)
      else
          a=8.0
          b=7.1447
          c=-4.125
          d=116.5
      ujl(i,j)=b
      vjl(i,j)=c
      pjl(i,j)=d
      rjl(i,j)=a
      end if
 
      end do
      
      j=jm
      
! upper boundary ʱ����صı߽����������жϼ���λ�ã��پ����ü���ǰ���߼�����Ĳ�����
      
      do i=ib,im-1
      
         xx=0.5*(x(i,j)+x(i+1,j))
         yy=0.5*(y(i,j)+y(i+1,j))
         x0=1/6.+1.73205/3.+ttime*20*1.73205/3.
         flag=1.73205*(xx-x0)-yy+1.
         
         if(flag.gt.0) then

             a=1.4
             b=0.
             c=0.
             d=1.0

           else

             a=8.0
             b=7.1447
             c=-4.125
             d=116.5

           endif
            
      
      ujr(i,j)=b
      vjr(i,j)=c
      pjr(i,j)=d
      rjr(i,j)=a
 
      
      end do

       
      return
      end
     