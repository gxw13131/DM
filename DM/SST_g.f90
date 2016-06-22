subroutine SST_Reconstruct_MUSCL
    use main
    
    REAL*8 :: dl,dr
    ! interpolate the turbulence variables (k/w) to cell interfaces
      
      ! iface
     
      do j=jb,jm-1
      do i=ib,im
      ! linear average to interpolate,
      ! distance weight is used to deal with non-homogeneous grid  
          alsinv=1./(L_Cell_x(i,j)+L_Cell_x(i-1,j))
          KTil (i,j)=(L_Cell_x(i,j)*KT(i-1,j)+L_Cell_x(i-1,j)*KT(i,j))&
     &                  *alsinv
          OmegaTil (i,j)=(L_Cell_x(i,j)*OmegaT(i-1,j)+L_Cell_x(i-1,j)*OmegaT(i,j))&
     &                  *alsinv
      
     KTir(i,j)=KTil(i,j)
     OmegaTir(i,j)=OmegaTil(i,j)
      end do
      end do
!
!     j face
!
     
      do j=jb,jm
      do i=ib,im-1
          alsinv=1./(L_Cell_y(i,j)+L_Cell_y(i,j-1))
          KTjl(i,j)=(L_Cell_y(i,j)*KT(i,j-1)+L_Cell_y(i,j-1)*KT(i,j))&
     &                  *alsinv
          OmegaTjl(i,j)=(L_Cell_y(i,j)*OmegaT(i,j-1)+L_Cell_y(i,j-1)*OmegaT(i,j))&
     &                  *alsinv
     
     KTjr(i,j)=KTjl(i,j)
     OmegaTjr(i,j)=OmegaTjl(i,j)
      end do
      end do
  
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      ! limiting the KT/OmegaT
      do j=jb,jm-1
       do i=ib,im-1
      
          !limiting KT
      
          dr=(KTil(i+1,j)-  KT(i,j))*2.
          dl=(KT(i,j)    -KTir(i,j))*2. !KTir(i,j)=KTil(i,j)
      
          KTir(i,j)  =KT(i,j)-0.5*alimiter(dl,dr)
          KTil(i+1,j)=KT(i,j)+0.5*alimiter(dr,dl)
      
          dr=(KTjl(i,j+1)-  KT(i,j))*2.
          dl=(KT(i,j)    -KTjr(i,j))*2.
      
          KTjr(i,j)  =KT(i,j)-0.5*alimiter(dl,dr)
          KTjl(i,j+1)=KT(i,j)+0.5*alimiter(dr,dl)
      
      
          !limiting OmegaT
      
          dr=(OmegaTil(i+1,j)-  OmegaT(i,j))*2.
          dl=(OmegaT(i,j)    -OmegaTir(i,j))*2.
      
          OmegaTir(i,j)  =OmegaT(i,j)-0.5*alimiter(dl,dr)
          OmegaTil(i+1,j)=OmegaT(i,j)+0.5*alimiter(dr,dl)
      
          dr=(OmegaTjl(i,j+1)-  OmegaT(i,j))*2.
          dl=(OmegaT(i,j)    -OmegaTjr(i,j))*2.
      
          OmegaTjr(i,j)  =OmegaT(i,j)-0.5*alimiter(dl,dr)
          OmegaTjl(i,j+1)=OmegaT(i,j)+0.5*alimiter(dr,dl)
       end do
      end do
         
 ! limiter, the same 3rd MUSCL     
contains
      function  alimiter(var1,var2)
      real*8 :: var1,var2
      
      real*8 :: alimiter
       alimiter=(var1*(var2*var2+2.*epsm*epsm)+&
     &                     var2*(2.*var1*var1+epsm*epsm))/&
     &            (2.*var1*var1-var1*var2+2.*var2*var2+3.*epsm*epsm) &
     &            *0.5*abs(sign(1.0,var1)+sign(1.0,var2)) 
     
      end function alimiter  
      
      
end subroutine  SST_Reconstruct_MUSCL     

subroutine SST_Flux_Lax
use main
REAL*8 :: ul,ur,vl,vr,KTl,KTr,OmegaTl,OmegaTr,rl,rr
REAL*8 :: sav1,sav2
! evaluate the  fluxes   
!
! -   i direction
      do j=jb,jm-1
       do i=ib,im    
   
!     left&right states      
      ul=uil(i,j)
      ur=uir(i,j) 
     
      vl=vil(i,j)
      vr=vir(i,j)
     
      KTl=max(KTil(i,j),EPSILON)
      KTr=max(KTir(i,j),EPSILON)
      
      OmegaTl=max(OmegaTil(i,j),EPSILON)
      OmegaTr=max(OmegaTir(i,j),EPSILON)
      
      rl=max(ril(i,j),EPSILON)
      rr=max(rir(i,j),EPSILON)
     
!     surface area vectors
!
      sav1=Sn_X_i(i,j)
      sav2=Sn_Y_i(i,j)
      

      uln=(sav1*ul+sav2*vl+vib(i,j))
      urn=(sav1*ur+sav2*vr+vib(i,j))
      
    
      eigenT=0.5*max(abs(uln),abs(urn)) !
    
!     flux terms
!
      dKT       =-0.5*(uln*KTl*rl       +urn*rr*KTr     -eigenT*(rr*KTr-rl*KTl))
      dOmegaT   =-0.5*(uln*OmegaTl*rl   +urn*rr*OmegaTr -eigenT*(rr*OmegaTr-rl*OmegaTl))
      
     
      F_KT(i,j)=F_KT(i,j)-dKT
      F_OmegaT(i,j)=F_OmegaT(i,j)-dOmegaT
     
      F_KT(i-1,j)=F_KT(i-1,j)+dKT
      F_OmegaT(i-1,j)=F_OmegaT(i-1,j)+dOmegaT
      
      end do
      end do
      
! -   j direction
     
!
      do j=jb,jm
       do i=ib,im-1
      
      
!     left&right states      
      ul=ujl(i,j)
      ur=ujr(i,j) 
     
      vl=vjl(i,j)
      vr=vjr(i,j)
     
      KTl=max(KTjl(i,j),EPSILON)
      KTr=max(KTjr(i,j),EPSILON)
      
      OmegaTl=max(OmegaTjl(i,j),EPSILON)
      OmegaTr=max(OmegaTjr(i,j),EPSILON)
      
      rl=max(rjl(i,j),EPSILON)
      rr=max(rjr(i,j),EPSILON)
     
!     surface area vectors
!
      sav1=Sn_X_j(i,j)
      sav2=Sn_Y_j(i,j)
      
      uln=(sav1*ul+sav2*vl+vjb(i,j))
      urn=(sav1*ur+sav2*vr+vjb(i,j))
          
      eigenT=0.5*max(abs(uln),abs(urn)) !
    
!     flux terms
!
      dKT       =-0.5*(uln*KTl*rl       +urn*rr*KTr     -eigenT*(rr*KTr-rl*KTl))
      dOmegaT   =-0.5*(uln*OmegaTl*rl   +urn*rr*OmegaTr -eigenT*(rr*OmegaTr-rl*OmegaTl))
      
     
      F_KT(i,j)=F_KT(i,j)-dKT
      F_OmegaT(i,j)=F_OmegaT(i,j)-dOmegaT
     
      F_KT(i,j-1)=F_KT(i,j-1)+dKT
      F_OmegaT(i,j-1)=F_OmegaT(i,j-1)+dOmegaT
      
      end do
      end do     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      return
end subroutine SST_Flux_Lax

subroutine SST_Diffusion
use main
real*8 :: KT_Df,OmegaT_Df
real*8 :: sigmaK=1/0.5
real*8 :: sigmaO=1/0.5
! -   i direction
      do j=jb,jm-1
      do i=ib,im
     
      sav1=Sn_X_i(i,j) !dy
      sav2=Sn_Y_i(i,j) !dx
      
      KT_Df=(Mu_L(i,j)+Mu_T(i,j)/sigmaK)*&
      &(sav1*dx_i(KT,i,j)+sav2*dy_i(KT,i,j))
      OmegaT_Df=(Mu_L(i,j)+Mu_T(i,j)/sigmaO)*&
      &(sav1*dx_i(KT,i,j)+sav2*dy_i(KT,i,j))
      
      F_KT(i,j)=F_KT(i,j)-KT_Df
      F_OmegaT(i,j)=F_OmegaT(i,j)-OmegaT_Df
      
      F_KT(i-1,j)=F_KT(i-1,j)+KT_Df
      F_OmegaT(i-1,j)=F_OmegaT(i-1,j)+OmegaT_Df
      end do
      end do

! -   j direction
      do j=jb,jm
      do i=ib,im-1 
      sav1=Sn_X_j(i,j)
      sav2=Sn_Y_j(i,j)  
      
      KT_Df=(Mu_L(i,j)+Mu_T(i,j)/sigmaK)*&
      &(sav1*dx_j(KT,i,j)+sav2*dy_j(KT,i,j))
      
      OmegaT_Df=(Mu_L(i,j)+Mu_T(i,j)/sigmaO)*&
      &(sav1*dx_j(OmegaT,i,j)+sav2*dy_j(OmegaT,i,j))
      
      F_KT(i,j)=F_KT(i,j)-KT_Df
      F_OmegaT(i,j)=F_OmegaT(i,j)-OmegaT_Df
      
      F_KT(i,j-1)=F_KT(i,j-1)+KT_Df
      F_OmegaT(i,j-1)=F_OmegaT(i,j-1)+OmegaT_Df
      end do
      end do


end subroutine SST_Diffusion

subroutine SST_LUSGS
use main
 use main
      
      real*8 :: dt_physics=1E30  !physics time step
      real*8 :: Ng(iq,jq),dt(iq,jq),beta(iq,jq)
      real*8 :: temp
      real*8 :: dKTI,dKTJ,dOmegaTI,dOmegaTJ
      real*8 :: dKTm,dOmegaTm
      real*8 :: SDx,SDy
      real*8 :: Vx_g,Vy_g,Vn !temp variables, updated each loop
      do j=jb,jm-1
      do i=ib,im-1
      
      dt(i,j)=1.0/(1.0/step(i,j)+1.5/dt_physics)
      beta(i,j)=dt(i,j)/Vcell(i,j)
      Ng(i,j)=1.0+(radius_i(i,j)+radius_j(i,j))*beta(i,j)
      
      !F_* means -Rt(U_(i,j))
      F_KT(i,j)    =   F_KT(i,j)&
                       &-(1.5*KT(i,j)      -2.0*KT_m1(i,j)    +0.5*KT_m2(i,j))&
                       &*Vcell(i,j)/dt_physics
      F_OmegaT(i,j) =   F_OmegaT(i,j)&
                       &-(1.5*OmegaT(i,j)   -2.0*OmegaT_m1(i,j) +0.5*OmegaT_m2(i,j))&
                       &*Vcell(i,j)/dt_physics
       !  Source Item
                       
        
        
       !  Dissipation Item 
                       
                       
      end do
      end do
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   LOWER 
      
      ! boundary DU=0
      j=jb
      do i=ib,im-1
      temp=beta(i,j)/Ng(i,j)
      
      F_KT(i,j)    =F_KT(i,j)    *temp
      F_OmegaT(i,j) =F_OmegaT(i,j) *temp
      
      end do
      
      i=ib
      do j=jb,jm-1
      temp=beta(i,j)/Ng(i,j)
      
      F_KT(i,j)    =F_KT(i,j)    *temp
      F_OmegaT(i,j) =F_OmegaT(i,j) *temp
      
      end do
      
      
      do j=jb+1,jm-1
      do i=ib+1,im-1
      ! I direction
      SDx=0.5*(Sn_X_i(i,j)+Sn_X_i(i-1,j))
      SDy=0.5*(Sn_Y_i(i,j)+SN_Y_i(i-1,j))
      
      Vx_g  =Vx(i-1,j)
      Vy_g  =Vy(i-1,j)
      
      Vn=SDx*Vx_g+SDy*Vy_g
      
      dKTm     =F_KT(i-1,j)
      dOmegaTm   =F_OmegaT(i-1,j)
      
      dKTI     =0.5*(Vn*dKTm     +radius_i(i-1,j)*dKTm)
      dOmegaTI   =0.5*(VN*dOmegaTm   +radius_i(i-1,j)*dOmegaTm)
      
     ! J direction
      SDx=0.5*(Sn_X_j(i,j)+Sn_X_j(i,j-1))
      SDy=0.5*(Sn_Y_j(i,j)+Sn_Y_j(i,j-1))
      
      Vx_g  =Vx(i,j-1)
      Vy_g  =Vy(i,j-1)
      
      Vn=SDx*Vx_g+SDy*Vy_g
      
      dKTm     =F_KT(i,j-1)
      dOmegaTm   =F_OmegaT(i,j-1)
      
      dKTJ     =0.5*(Vn*dKTm      +radius_j(i,j-1)*dKTm)
      dOmegaTJ   =0.5*(Vn*dOmegaTm   +radius_j(i,j-1)*dOmegaTm)
      
      temp=beta(i,j)/Ng(i,j)
      
      F_KT(i,j)    =(F_KT(i,j)    +(dKTI+dKTJ)      )*temp
      F_OmegaT(i,j) =(F_OmegaT(i,j) +(dOmegaTI+dOmegaTJ)  )*temp
      
      end do
      end do
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   UPPER 
     
      !jm=jm-1
      !do i=im-1,ib
      !
      !end do
      !
      !i=im-1
      !do j=jm-1,jb
      !
      !end do
      
      do j=jm-2,jb
      do i=im-2,ib
       ! I direction
      SDx=0.5*(Sn_X_i(i+2,j)+Sn_X_i(i+1,j))
      SDy=0.5*(Sn_Y_i(i+2,j)+SN_Y_i(i+1,j))
      
      Vx_g  =Vx(i+1,j)
      Vy_g  =Vy(i+1,j)
      
      Vn=SDx*Vx_g+SDy*Vy_g
      
      dKTm     =F_KT(i+1,j)
      dOmegaTm  =F_OmegaT(i+1,j)
      
      dKTI     =0.5*(Vn*dKTm     -radius_i(i+1,j)*dKTm)
      dOmegaTI   =0.5*(Vn*dOmegaTm   -radius_i(i+1,j)*dOmegaTm)
      
       ! J direction
      SDx=0.5*(Sn_X_j(i,j+1)+Sn_X_j(i,j+2))
      SDy=0.5*(Sn_Y_j(i,j+1)+Sn_Y_j(i,j+2))
      
      Vx_g  =Vx(i,j+1)
      Vy_g  =Vy(i,j+1)
      
      Vn=SDx*Vx_g+SDy*Vy_g
      
      dKTm     =F_KT(i,j+1)
      dOmegaTm  =F_OmegaT(i,j+1)
      
      
      dKTJ     =0.5*(Vn*dKTm     -radius_j(i,j+1)*dKTm)
      dOmegaTJ   =0.5*(Vn*dOmegaTm   -radius_j(i,j+1)*dOmegaTm)
     
    
      temp=beta(i,j)/Ng(i,j)
      
      F_KT(i,j)    =F_KT(i,j)    -(dKTI+dKTJ)      *temp
      F_OmegaT(i,j) =F_OmegaT(i,j) -(dOmegaTI+dOmegaTJ)  *temp
      
      end do
      end do 
      
      do j=jb,jm
      do i=ib,im
      
      KT(i,j)       =KT(i,j)   +F_KT(i,j)
      OmegaT(i,j)   =OmegaT(i,j)+F_OmegaT(i,j)
      
      end do
      end do
      


end subroutine SST_LUSGS
