subroutine LUSGS()
!     *****************
!
      use main
      
      real*8 :: dt_physics=1E20  !physics time step
      real*8 :: Ng(iq,jq),dt(iq,jq)
      real*8 :: temp
      real*8 :: AdRho,AdRhoVx,AdRhoVy,AdRhoEt 
      real*8 :: dRhoI,dRhoVxI,dRhoVyI,dRhoEtI
      real*8 :: dRhoJ,dRhoVxJ,dRhoVyJ,dRhoEtJ
      real*8 :: dRhoM,dRhoVxM,dRhoVyM,dRhoEtM
      real*8 :: SDx,SDy
      real*8 :: Vx_g,Vy_g,p_g,rho_g !temp variables, updated each loop
      do j=jb,jm-1
      do i=ib,im-1
      
      dt(i,j)=1.0/(1.0/step(i,j)+1.5/dt_physics)
      Ng(i,j)=1.0+(radius_i(i,j)+radius_j(i,j))*dt(i,j)/Vcell(i,j)
      
      F_Rho(i,j)    =F_Rho(i,j)/Vcell(i,j)&
                       &-(1.5*Rho(i,j)      -2.0*Rho_m1(i,j)    +0.5*Rho_m2(i,j))/dt_physics
      F_Rho_vx(i,j) =F_Rho_vx(i,j)/Vcell(i,j)&
                       &-(1.5*Rho_vx(i,j)   -2.0*Rho_vx_m1(i,j) +0.5*Rho_vx_m2(i,j))/dt_physics
      F_Rho_vy(i,j) =F_Rho_vy(i,j)/Vcell(i,j)&
                       &-(1.5*Rho_vy(i,j)   -2.0*Rho_vy_m1(i,j) +0.5*Rho_vy_m2(i,j))/dt_physics
      F_Rho_Et(i,j) =F_Rho_Et(i,j)/Vcell(i,j)&
                       &-(1.5*Rho_Et(i,j)   -2.0*Rho_Et_m1(i,j) +0.5*Rho_Et_m2(i,j))/dt_physics
      end do
      end do
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   LOWER 
      
      do j=jb,jm-1
      do i=ib,im-1
      if(i==ib) then
      dRhoI=0.0
      dRhoVxI=0.0
      dRhoVyI=0.0
      dRhoEtI=0.0
      else
      SDx=0.5*(Sn_X_i(i,j)+Sn_X_i(i-1,j))
      SDy=0.5*(Sn_Y_i(i,j)+SN_Y_i(i-1,j))
      
      Vx_g  =Vx(i-1,j)
      Vy_g  =Vy(i-1,j)
      p_g   =p(i-1,j)
      rho_g =rho(i-1,j)
      
      dRhoM     =F_Rho(i-1,j)
      dRhoVxM   =F_Rho_Vx(i-1,j)
      dRhoVyM   =F_Rho_Vy(i-1,j)
      dRhoEtM   =F_Rho_Et(i-1,j)
      
      call AdU(Vx_g,Vy_g,SDx,SDy,p_g,rho_g,dRhoM,dRhoVxM,dRhoVyM,dRhoEtM)
      dRhoI     =0.5*(AdRho     +radius_i(i-1,j)*dRhoM)
      dRhoVxI   =0.5*(AdRhoVx   +radius_i(i-1,j)*dRhoVxM)
      dRhovyI   =0.5*(AdRhoVy   +radius_i(i-1,j)*dRhoVyM)
      dRhoEtI   =0.5*(AdRhoEt   +radius_i(i-1,j)*dRhoEtM)
      end if
      
      if(j==jb) then
      dRhoJ     =0.0
      dRhoVxJ   =0.0
      dRhoVyJ   =0.0
      dRhoEtJ   =0.0
      else 
      SDx=0.5*(Sn_X_j(i,j)+Sn_X_j(i,j-1))
      SDy=0.5*(Sn_Y_j(i,j)+Sn_Y_j(i,j-1))
      
      Vx_g  =Vx(i,j-1)
      Vy_g  =Vy(i,j-1)
      p_g   =p(i,j-1)
      rho_g =rho(i,j-1)
      
      dRhoM     =F_Rho(i,j-1)
      dRhoVxM   =F_Rho_Vx(i,j-1)
      dRhoVyM   =F_Rho_Vy(i,j-1)
      dRhoEtM   =F_Rho_Et(i,j-1)
      
      call AdU(Vx_g,Vy_g,SDx,SDy,p_g,rho_g,dRhoM,dRhoVxM,dRhoVyM,dRhoEtM)
      
      dRhoJ     =0.5*(AdRho     +radius_j(i,j-1)*dRhoM)
      dRhoVxJ   =0.5*(AdRhoVx   +radius_j(i,j-1)*dRhoVxM)
      dRhoVyJ   =0.5*(AdRhoVy   +radius_j(i,j-1)*dRhoVyM)
      dRhoEtJ   =0.5*(AdRhoEt   +radius_j(i,j-1)*dRhoEtM)
      end if
      temp=dt(i,j)/Ng(i,j)
      
      F_Rho(i,j)    =(F_Rho(i,j)    +(dRhoI+dRhoJ)      /Vcell(i,j))*temp
      F_Rho_vx(i,j) =(F_Rho_vx(i,j) +(dRhoVxI+dRhoVxJ)  /Vcell(i,j))*temp
      F_Rho_vy(i,j) =(F_Rho_vy(i,j) +(dRhoVyI+dRhoVyJ)  /Vcell(i,j))*temp
      F_Rho_Et(i,j) =(F_Rho_Et(i,j) +(dRhoEtI+dRhoEtJ)  /Vcell(i,j))*temp
      end do
      end do
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   UPPER 
     
      do j=jm-1,jb
      do i=im-1,ib
      
      
      if(i==im-1)then
      dRhoI=0.0
      dRhoVxI=0.0
      dRhoVyI=0.0
      dRhoEtI=0.0
      else
      SDx=0.5*(Sn_X_i(i+2,j)+Sn_X_i(i+1,j))
      SDy=0.5*(Sn_Y_i(i+2,j)+SN_Y_i(i+1,j))
      
      Vx_g  =Vx(i+1,j)
      Vy_g  =Vy(i+1,j)
      p_g   =p(i+1,j)
      rho_g =rho(i+1,j)
      
      dRhoM     =F_Rho(i+1,j)
      dRhoVxM   =F_Rho_Vx(i+1,j)
      dRhoVyM   =F_Rho_Vy(i+1,j)
      dRhoEtM   =F_Rho_Et(i+1,j)
      
      call AdU(Vx_g,Vy_g,SDx,SDy,p_g,rho_g,dRhoM,dRhoVxM,dRhoVyM,dRhoEtM)
      
      dRhoI     =0.5*(AdRho     -radius_i(i+1,j)*dRhoM)
      dRhoVxI   =0.5*(AdRhoVx   -radius_i(i+1,j)*dRhoVxM)
      dRhovyI   =0.5*(AdRhoVy   -radius_i(i+1,j)*dRhoVyM)
      dRhoEtI   =0.5*(AdRhoEt   -radius_i(i+1,j)*dRhoEtM)
      
      end if    !i direction
      
      if(j==jm-1)then 
      dRhoJ=0.0
      dRhoVxJ=0.0
      dRhoVyJ=0.0
      dRhoEtJ=0.0
      else
      SDx=0.5*(Sn_X_j(i,j+1)+Sn_X_j(i,j+2))
      SDy=0.5*(Sn_Y_j(i,j+1)+Sn_Y_j(i,j+2))
      
      Vx_g  =Vx(i,j+1)
      Vy_g  =Vy(i,j+1)
      p_g   =p(i,j+1)
      rho_g =rho(i,j+1)
      
      dRhoM     =F_Rho(i,j+1)
      dRhoVxM   =F_Rho_Vx(i,j+1)
      dRhoVyM   =F_Rho_Vy(i,j+1)
      dRhoEtM   =F_Rho_Et(i,j+1)
      call AdU(Vx_g,Vy_g,SDx,SDy,p_g,rho_g,dRhoM,dRhoVxM,dRhoVyM,dRhoEtM)
      
      dRhoJ     =0.5*(AdRho     -radius_j(i,j+1)*dRhoM)
      dRhoVxJ   =0.5*(AdRhoVx   -radius_j(i,j+1)*dRhoVxM)
      dRhoVyJ   =0.5*(AdRhoVy   -radius_j(i,j+1)*dRhoVyM)
      dRhoEtJ   =0.5*(AdRhoEt   -radius_j(i,j+1)*dRhoEtM)
      
      end if    !j direction
      
      temp=dt(i,j)/Vcell(i,j)/Ng(i,j)
      
      F_Rho(i,j)    =F_Rho(i,j)    -(dRhoI+dRhoJ)      *temp
      F_Rho_vx(i,j) =F_Rho_vx(i,j) -(dRhoVxI+dRhoVxJ)  *temp
      F_Rho_vy(i,j) =F_Rho_vy(i,j) -(dRhoVyI+dRhoVyJ)  *temp
      F_Rho_Et(i,j) =F_Rho_Et(i,j) -(dRhoEtI+dRhoEtJ)  *temp
      
      end do
      end do 
      
      do j=jb,jm
      do i=ib,im
      
      Rho(i,j)      =Rho(i,j)   +F_Rho(i,j)
      Rho_Vx(i,j)   =Rho_Vx(i,j)+F_Rho_Vx(i,j)
      Rho_Vy(i,j)   =Rho_Vy(i,j)+F_Rho_Vy(i,j)
      Rho_Et(i,j)   =Rho_Et(i,j)+F_Rho_Et(i,j)
      end do
      end do
      
      
      
      
 contains
 subroutine AdU(Vx_,Vy_,SDx_,SDy_,p_,rho_,dRhoM_,dRhoVxM_,dRhoVyM_,dRhoEtM_)
 implicit none
 real*8 :: Vx_,Vy_,SDx_,SDy_,p_,rho_
 real*8 :: dRhoM_,dRhoVxM_,dRhoVyM_,dRhoEtM_
 
 real*8 :: Vm2,Vn
 real*8 :: Gm,Gm1,Gm2,phi,e1,e2
 
 Gm=Gamma
 Gm1=Gm-1.0
 Gm2=2.0-Gm
 
 Vn=Vx_*SDx_+Vy_*SDy_ ! normal velocity
 Vm2=Vx_*Vx_+Vy_*Vy_ !  square of velocity
 phi=0.5*Gm1*Vm2
 e1=-0.5*Gm2*Vm2-Gm*p_/rho_/Gm1
 e2=0.5*Vm2+Gm*p_/rho_/Gm1
 
 AdRho=      dRhoM_     *(0.0)+&
            dRhoVxM_    *(SDx_)+ &
            dRhoVyM_    *(SDy_)+ &
            dRhoEtM_    *(0.0)
 AdRhoVx=    dRhoM_     *(SDx_*phi-Vx_*Vn)+&
            dRhoVxM_    *(Vn+SDx_*Gm2*Vx_)+ &
            dRhoVyM_    *(SDy_*Vx_-SDx_*Gm1*Vy_)+ &
            dRhoEtM_    *(SDx_*Gm1)
 AdRhoVy=    dRhoM_     *(SDy_*phi-Vy_*Vn)+&
            dRhoVxM_    *(SDx_*Vy_-SDy_*Gm1*Vx_)+ &
            dRhoVyM_    *(Vn+SDy_*Gm2*Vy_)+ &
            dRhoEtM_    *(SDy_*Gm1)
 AdRhoEt=    dRhoM_     *(e1*Vn)+&
            dRhoVxM_    *(SDx_*e2-Gm1*Vx_*Vn)+ &
            dRhoVyM_    *(SDy_*e2-Gm1*Vy_*Vn)+ &
            dRhoEtM_    *(Gm*Vn)
 
 end subroutine Adu
      
end subroutine