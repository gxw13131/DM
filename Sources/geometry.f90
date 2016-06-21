
!==============================================================================
      subroutine area
!     ***************
!     area vectors
!     ***************
      
!
      use main
!
!
!     this subrouitne computes cell face areas and interpolation coefficients.
!
!
!
!=====compute the coordinates of cell centers      
      
      real*8 :: Xw,Yw,Xe,Ye,Xn,Yn,Xs,Ys
 
      do i=ib-1,im
      do j=jb-1,jm
      
          xc(i,j)=0.25*(x(i,j)+x(i+1,j)&
     &                      +x(i,j+1)+x(i+1,j+1))
          yc(i,j)=0.25*(y(i,j)+y(i+1,j)&
     &                      +y(i,j+1)+y(i+1,j+1))       
      end do
      end do
     
!
!     i face: surface vector components, grid velocity, shifted volume
!
      do i=ib-1,im+1
      do j=jb-1,jm
      	
          Sn_X_i(i,j)=(y(i,j+1)-y(i,j))  !Ãæ»ýÊ¸Á¿
          Sn_Y_i(i,j)=-(x(i,j+1)-x(i,j))
          vib(i,j)=0. !grid velocity, no use for now
          
          S_i(i,j)=sqrt(Sn_X_i(i,j)**2+Sn_Y_i(i,j)**2)
     
      end do
      end do
!
!     j face: surface vector components, grid velocity, shifted volume
!
      do i=ib-1,im
      do j=jb-1,jm+1
          
          Sn_X_j(i,j)=-(y(i+1,j)-y(i,j))
          Sn_Y_j(i,j)=(x(i+1,j)-x(i,j))
          vjb(i,j)=0.0
          
          S_j(i,j)=sqrt(Sn_X_j(i,j)**2+Sn_Y_j(i,j)**2)
      end do
      end do
      
      
!     interpolation coefficient

      do i=ib-1,im
      do j=jb-1,jm
 !      interface pisition     
      xw=0.5*(x(i,j)+x(i,j+1))
      yw=0.5*(y(i,j)+y(i,j+1))
     
      xe=0.5*(x(i+1,j)+x(i+1,j+1))
      ye=0.5*(y(i+1,j)+y(i+1,j+1))
     
      xs=0.5*(x(i,j)+x(i+1,j))
      ys=0.5*(y(i,j)+y(i+1,j))
     
      xn=0.5*(x(i,j+1)+x(i+1,j+1)) 
      yn=0.5*(y(i,j+1)+y(i+1,j+1)) 
     
 !      distance between two parallel interface    
      L_Cell_x(i,j)=sqrt((xw-xe)**2+(yw-ye)**2)
       
      L_Cell_y(i,j)=sqrt((xs-xn)**2+(ys-yn)**2)
       
      
      end do
      end do
     
      return
      end
     
!==============================================================================
      subroutine volume_static
!     **********************************
!     compute volume through coordinates
!     **********************************      
!

      use main

      do i=ib-1,im
      do j=jb-1,jm
      
          Vcell(i,j)=0.5*((x(i+1,j+1)-x(i,j))*&
     &                    (y(i,j+1)-y(i+1,j))-&
     &                    (y(i+1,j+1)-y(i,j))*&
     &                    (x(i,j+1)-x(i+1,j)))
          
          if(Vcell(i,j).lt.0.) then
                 write(*,*) 'volume.lt.0',i,j
                 stop
          end if
          
      end do
      end do
!
      return
      end
      
