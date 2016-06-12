subroutine inviscid_fluxes_AUSMPWplus
!     =========================
      
      use main
 
      real*8 :: uL,uR,vL,vR,pL,pR,rhoL,rhoR
      real*8 :: SNx,SNy
      real*8 :: dRho,dRhoVx,dRhoVy,dRhoEt
 
      
 !=======================================
      ! reset the flux to zero     
      do i=ib-1,im+1
      do j=jb-1,jm+1
      
      F_rho(i,j)=0.
      F_rho_Et(i,j)=0.
      F_rho_vx(i,j)=0.
      F_rho_vy(i,j)=0.     
      end do
      end do
     
! evaluate the inviscid fluxes   
!
! -   i direction
!
!
      do j=jb,jm-1
      do i=ib,im
        
!     left&right states
     
      uL=uil(i,j)
      uR=uir(i,j) 
     
      vL=vil(i,j)
      vR=vir(i,j)
     
      pL=max(pil(i,j),EPSILON)
      pR=max(pir(i,j),EPSILON)
     
      rhoL=max(ril(i,j),EPSILON)
      rhoR=max(rir(i,j),EPSILON)
      
!     surface area vectors
!
      SNx=Sn_X_i(i,j)
      SNy=Sn_Y_i(i,j)
             
          
!     flux terms
    
      call AUSMPWplus(uL,vL,pL,rhoL,uR,vR,pR,rhoR,SNx,SNy,&
                        &dRho,dRhoVx,dRhoVy,dRhoEt)
     
      F_rho(i,j)=F_rho(i,j)-dRho
      F_rho_Et(i,j)=F_rho_Et(i,j)-dRhoEt
      F_rho_vx(i,j)=F_rho_vx(i,j)-dRhoVx
      F_rho_vy(i,j)=F_rho_vy(i,j)-dRhoVy
     
      F_rho(i-1,j)=F_rho(i-1,j)+dRho
      F_rho_Et(i-1,j)=F_rho_Et(i-1,j)+dRhoEt
      F_rho_vx(i-1,j)=F_rho_vx(i-1,j)+dRhoVx
      F_rho_vy(i-1,j)=F_rho_vy(i-1,j)+dRhoVy
      
      
      
      end do
      end do
      
!
! -   j direction
!
      
!
      do j=jb,jm
      do i=ib,im-1
      
      
!     left&right states
! 
      uL=ujl(i,j)
      uR=ujr(i,j) 
     
      vL=vjl(i,j)
      vR=vjr(i,j)
     
      pL=max(pjl(i,j),EPSILON)
      pR=max(pjr(i,j),EPSILON)
     
      rhoL=max(rjl(i,j),EPSILON)
      rhoR=max(rjr(i,j),EPSILON)
          

      !     surface area vectors
      SNx=Sn_X_j(i,j)
      SNy=Sn_Y_j(i,j)
      
      call AUSMPWplus(uL,vL,pL,rhoL,uR,vR,pR,rhoR,SNx,SNy,&
                        &dRho,dRhoVx,dRhoVy,dRhoEt)
      
      F_rho(i,j)=F_rho(i,j)-dRho
      F_rho_Et(i,j)=F_rho_Et(i,j)-dRhoEt
      F_rho_vx(i,j)=F_rho_vx(i,j)-dRhoVx
      F_rho_vy(i,j)=F_rho_vy(i,j)-dRhoVy
      
      F_rho(i,j-1)=F_rho(i,j-1)+dRho
      F_rho_Et(i,j-1)=F_rho_Et(i,j-1)+dRhoEt
      F_rho_vx(i,j-1)=F_rho_vx(i,j-1)+dRhoVx
      F_rho_vy(i,j-1)=F_rho_vy(i,j-1)+dRhoVy
      
      
      
      end do
      end do     
!  

contains 
    subroutine AUSMPWplus(uL,vL,pL,rhoL,uR,vR,pR,rhoR,SNx,SNy,&
    &dRho,dRhoVx,dRhoVy,dRhoEt)
        implicit none
        real*8 :: EPSILON=1E-10
        real*8 :: alpha=3.0/16.0
        real*8 :: uL,vL,pL,rhoL,uR,vR,pR,rhoR
        real*8 :: SNx,SNy,SNxN,SNyN

        real*8 :: Vm2L,Vm2R,HL,HR,VnL,VnR
        real*8 :: Gm,Gm1,Gm2
        
        !intermediate variables
        real*8 :: aStarL,aStarR,aTildeL,aTildeR,aS
        real*8 :: MpL,MnR,PpL,pnR,MpL_,MnR_
        real*8 :: fL,fR
        real*8 :: ps,mHalf,wp
        
        Gm=Gamma
        Gm1=Gm/(Gm-1.0)
        Gm2=2.0*(Gm-1.0)/(Gm+1.0)
        ! magnitude^2 of velocity
        Vm2L=uL**2+vL**2 
        Vm2R=uR**2+vR**2
        ! enthalpy
        HL=Gm1*pL/rhoL+0.5*Vm2L
        HR=Gm1*pR/rhoR+0.5*Vm2R
        !sound speed
        aL=sqrt(Gm*pL/rhoL)
        aR=sqrt(Gm*pR/rhoR)
        ! normalized SD
        SNxN=SNx/sqrt(SNx**2+SNy**2)
        SNyN=SNy/sqrt(SNx**2+SNy**2)
        ! normal speed on interface
        VnL=(uL*SNxN+vL*SNyN)
        VnR=(uR*SNxN+vR*SNyN)
        ! critical sound speed
        aStarL=sqrt(Gm2*HL)
        aStarR=sqrt(Gm2*HR)
        ! numerical sound speed
        !aTildeL=aStarL**2/max(sqrt(Vm2L),aStarL)
        !aTildeR=aStarR**2/max(sqrt(Vm2R),aStarR)
        
        aTildeL=aStarL**2/max(abs(VnL),aStarL)
        aTildeR=aStarR**2/max(abs(VnR),aStarR)
        
        ! sound speed on the interface 
        aS=min(aTildeL,aTildeR)
        
        ! Mach number
        MaL=VnL/aS
        maR=VnR/aS
        
        !split function
        if(abs(MaL)>1.0)then
        MpL=0.5*(MaL+abs(MaL))
        !PpL=0.5*(1.0+sign(M))
        PpL=0.5*(MaL+abs(MaL))/(MaL+EPSILON)
        
        else
        MpL=0.25*(MaL+1.0)**2
        PpL=0.25*(MaL+1.0)**2*(2.0-MaL)+alpha*(MaL**2-1.0)**2
        
        endif
        
        if(abs(MaR)>1.0)then
        MnR=0.5*(MaR-abs(MaR))
        !PpL=0.5*(1.0+sign(M))
        PnR=0.5*(MaR-abs(MaR))/(MaR+EPSILON)
        else
        MnR=-0.25*(MaR-1.0)**2
        PnR=0.25*(MaR-1.0)**2*(2.0+MaR)-alpha*(MaR**2-1.0)**2
        endif
        
        pS=PpL*pL+pnR*pR
        mHalf=MpL+MnR
        
        wp=1.0-min(pL/pR,pR/pL)**3
        
        if(abs(MaL)<1.0)then
        fL=pL/(pS+EPSILON)-1.0
        else 
        fL=0.0
        endif

        if(abs(MaR)<1.0)then
        fR=pR/(pS+EPSILON)-1.0
        else
        fR=0.0
        endif
        
        if(mHalf>=0)then
        MpL_=MpL+MnR*((1.0-wp)*(1.0+fR)-fL)
        MnR_=MnR*wp*(1.0+fR)
        else
        MpL_=MnR*wp*(1.0+fL)
        MnR_=MnR+MpL*((1.0-wp)*(1.0+fL)-fR)
        endif
        
        dRho=MpL_*aS*rhoL+MnR_*aS*rhoR
        dRhoVx=MpL_*aS*rhoL*uL+MnR_*aS*rhoR*uR+(PpL*pL+PnR*pR)*SNxN
        dRhoVy=MpL_*aS*rhoL*vL+MnR_*aS*rhoR*vR+(PpL*pL+PnR*pR)*SNyN
        dRhoEt=MpL_*aS*rhoL*HL+MnR_*aS*rhoR*HR

    end

      
      end 
 