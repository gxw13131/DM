subroutine LUSGS()
!     *****************
!
      use main
      
      real*8 :: dt_lusgs=0.01  !vitual time step
      real*8 :: N(iq,jq),dt(iq,jq)
      real*8 :: temp
      real*8 :: AdRho,AdRhoVx,AdRhoVy,AdRhoEt 
      do j=jb,jm-1
      do i=ib,im-1
      
      dt(i,j)=1.0/(1.0/dt_lusgs+1.5/step(i,j))
      N(i,i)=1.0+(radius_i(i,j)+radius_j(i,j))/dt(i,j)/Vcell(i,j)
      F_Rho(i,j)=F_Rho(i,j)/Vcell(i,j)&
                       &-(1.5*Rho(i,j)-2.0*Rho_m1(i,j)+Rho_m2(i,j))/step(i,j)
      F_Rho_vx(i,j)=F_Rho_vx(i,j)/Vcell(i,j)&
                       &-(1.5*Rho_vx(i,j)-2.0*Rho_vx_m1(i,j)+Rho_vx_m2(i,j))/step(i,j)
      F_Rho_vy(i,j)=F_Rho_vy(i,j)/Vcell(i,j)&
                       &-(1.5*Rho_vy(i,j)-2.0*Rho_vy_m1(i,j)+Rho_vy_m2(i,j))/step(i,j)
      F_Rho_Et(i,j)=F_Rho_Et(i,j)/Vcell(i,j)&
                       &-(1.5*Rho_Et(i,j)-2.0*Rho_Et_m1(i,j)+Rho_Et_m2(i,j))/step(i,j)
      end do
      end do
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   LOWER 
      
      do j=jb,jm
      do i=ib,im
      if(i==ib) then
      dRhoI=0.0
      dRhoVxI=0.0
      dRhoVyI=0.0
      dRhoEtI=0.0
      else
      Vx=Vx(i-1,j)
      Vy=Vy(i-1,j)
      SDx=0.5*(Sn_X_i(i,j)+Sn_X_i(i-1,j))
      SDy=0.5*(Sn_Y_i(i,j)+SN_Y_i(i-1,j))
      
      dRhoM=
      
      call AdU(Vx,Vy,SDx,SDy,p,rho,dRhoM,dRhoVxM,dRhoVyM,dRhoEtM)
      dRhoI=0.5*(AdRho+radius_i(i-1,j)*F_rho(i-1,j))
      dRhoVxI=0.5*(AdRhoVx+radius_i(i-1,j)*F_rho_vx(i-1,j))
      dRhovyI=0.5*(AdRhoVy+radius_i(i-1,j)*F_rho_vy(i-1,j))
      dRhoEtI=0.5*(AdRhoEt+radius_i(i-1,j)*F_rho_Et(i-1,j))
      end if
      
      if(j==jb) then
      dRhoJ=0.0
      dRhoVxJ=0.0
      dRhoVyJ=0.0
      dRhoEtJ=0.0
      else 
      
      dRhoJ=0.5*(AdRho+radius_j(i-1,j)*F_rho(i,j-1))
      dRhoVxJ=0.5*(AdRhoVx+radius_j(i-1,j)*F_rho_vx(i,j-1))
      dRhoVyJ=0.5*(AdRhoVy+radius_j(i-1,j)*F_rho_vy(i,j-1))
      dRhoEtJ=0.5*(AdRhoEt+radius_j(i-1,j)*F_rho_Et(i,j-1))
      end if
      temp=dt(i,j)/N(i,j)
      F_Rho(i,j)=(F_Rho(i,j)+(dRhoI+dRhoJ)/Vcell(i,j))*temp
      F_Rho_vx(i,j)=(F_Rho_vx(i,j)+(dRhoVxI+dRhoVxJ)/Vcell(i,j))*temp
      F_Rho_vy(i,j)=(F_Rho_vy(i,j)+(dRhoVyI+dRhoVyJ)/Vcell(i,j))*temp
      F_Rho_Et(i,j)=(F_Rho_Et(i,j)+(dRhoEtI+dRhoEtJ)/Vcell(i,j))*temp
      end do
      end do
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   LOWER 
     
      do j=jm-1,jb
      do i=im-1,ib
      
      rho_m1(i,j)=0
      rho_vx_m1(i,j)=0
      rho_vy_m1(i,j)=0
      rho_Et_m1(i,j)=0
      
      
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
 
 dRho=      dRhoM_       *(0.0)+&
            dRhoVxM_     *(SDx_)+ &
            dRhoVyM_     *(SDy_)+ &
            dRhoEtM_     *(0.0)
 dRhoVx=    dRhoM_       *(SDx_*phi-Vx_*Vn)+&
            dRhoVxM_     *(Vn+SDx_*Gm2*Vx_)+ &
            dRhoVyM_     *(SDy_*Vx_-SDx_*Gm1*Vy_)+ &
            dRhoEtM_     *(SDx_*Gm1)
 dRhoVy=    dRhoM_       *(SDy_*phi-Vy_*Vn)+&
            dRhoVxM_     *(SDx_*Vy_-SDy_*Gm1*Vx_)+ &
            dRhoVyM_     *(Vn+SDy_*Gm2*Vy_)+ &
            dRhoEtM_     *(SDy_*Gm1)
 dRhoEt=    dRhoM_       *(e1*Vn)+&
            dRhoVxM_     *(SDx_*e2-Gm1*Vx_*Vn)+ &
            dRhoVyM_     *(SDy_*e2-Gm1*Vy_*Vn)+ &
            dRhoEtM_     *(Gm*Vn)
 
 end subroutine Adu
      
end subroutine