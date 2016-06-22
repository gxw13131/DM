      
      subroutine outputuns
!     ********************
      use main
      character*7 nam
      integer :: k
      
      real*8 :: Xcar, Ycar,Ucar,Vcar,Pcar,Rcar,entr,Amach,Vmur,TT
        if(mod(n,iprint)==0) then
        k=(n/iprint)
        if(k<1000) then
         
            if(k<10) then
               write(nam,'(a,i1,a)') '00',k,'.dat'
             else if(k<100) then
               write(nam,'(a,i2,a)') '0',k,'.dat'
             else
               write(nam,'(i3,a)') k,'.dat'
            end if
            
            open(30,file=nam)
            write(30,*) 'title="contour"'
            write(30,*) 'variables="x","y","u","v",'
            write(30,*)'"p","den","T","mach","entropy","mu","k","Omega"'
     
            write(30,*) 'zone i=',im-ib,' j=',jm-jb &
                   &,' f=point'
       
          do j=jb,jm-1
          do i=ib,im-1
             xcar=xc(i,j)*L_ref
             ycar=yc(i,j)*L_ref
             
             ucar=vx(i,j)*V_ref
             vcar=vy(i,j)*V_ref
             
             pcar=p(i,j)*Rho_ref*V_ref*V_ref
             rcar=rho(i,j)*Rho_ref
             TT=T(i,j)*T_ref
             KT_out=KT(i,j)
             OmegaT_out=OmegaT(i,j)
             entr=pcar/rcar**Gamma
             amach=sqrt(ucar*ucar+vcar*vcar)&
     &             /sqrt(Gamma*pcar/rcar)
             vmur=Mu_E(i,j)*Rho_ref*V_ref*L_ref
             
             write(30,"(E10.3,1X,  E10.3,1X,   E10.3,1X,  E10.3,1X,E10.3,1X,&
     &             E10.3,1X,  E10.3,1X,  E10.3,1X,  E10.3, 1X, E10.3, 1X, E10.3,1X, E10.3, 1X, E10.3)") xcar,ycar&
     &                  ,ucar,vcar&
     &                  ,pcar,rcar,TT,amach,entr,vmur,&
     &                  KT_out,OmegaT_out
           end do
           end do
           close(30)
        
        
        end if
        end if
      
      return
      end
      
      
      subroutine output
!     *****************
      use main
!
!     save data for restart
!     ---------------------
!
      if(is_print.eq.1) then
      write(*,*) 'now,backup the intermediate results'
      open(1,file='bl.sav',form='unformatted')
      write(1) n,ttime
      
      do  i=ib,im
      do  j=jb,jm
          write(1)    vx(i,j),vy(i,j)&
     &               ,p(i,j),rho(i,j)&
     &               ,Mu_L(i,j),Mu_E(i,j)&
     &               ,KT(i,j),OmegaT(i,j)&
     &               ,Rho_m1(i,j),Rho_Et_m1(i,j)&
     &               ,Rho_vx_m1(i,j),Rho_vy_m1(i,j)
      end do
      end do
     
      close(1)
      write(*,*) 'backup finish'
      end if
!
      return
      end

      
