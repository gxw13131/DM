subroutine inviscid_fluxes_AUSMPWplus
!     =========================
      
      use main
 
      real*8 :: uL,uR,vL,vR,pL,pR,rhoL,rhoR
      real*8 :: SNx,SNy,S
      real*8 :: dRho,dRhoVx,dRhoVy,dRhoEt
 
      
 !=======================================
      ! reset the flux to zero     
      do j=jb-1,jm+1
      do i=ib-1,im+1
      
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
      S=S_i(i,j)       
          
!     flux terms
    
      call AUSMPWplus(uL,vL,pL,rhoL,uR,vR,pR,rhoR,SNx,SNy)
      !if(i==100.AND.j==3)then
      !write(*,*) "pause"
      !endif
      F_rho(i,j)    =F_rho(i,j)     +dRho*S
      F_rho_Et(i,j) =F_rho_Et(i,j)  +dRhoEt*S
      F_rho_vx(i,j) =F_rho_vx(i,j)  +dRhoVx*S
      F_rho_vy(i,j) =F_rho_vy(i,j)  +dRhoVy*S
     
      F_rho(i-1,j)      =F_rho(i-1,j)   -dRho*S
      F_rho_Et(i-1,j)   =F_rho_Et(i-1,j)-dRhoEt*S
      F_rho_vx(i-1,j)   =F_rho_vx(i-1,j)-dRhoVx*S
      F_rho_vy(i-1,j)   =F_rho_vy(i-1,j)-dRhoVy*S
          
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
      S=S_j(i,j)
      
      call AUSMPWplus(uL,vL,pL,rhoL,uR,vR,pR,rhoR,SNx,SNy)
      !if(i==100.AND.j==4)then
      !write(*,*) "pause"
      !endif
      F_rho(i,j)    =F_rho(i,j)     +dRho*S
      F_rho_Et(i,j) =F_rho_Et(i,j)  +dRhoEt*S
      F_rho_vx(i,j) =F_rho_vx(i,j)  +dRhoVx*S
      F_rho_vy(i,j) =F_rho_vy(i,j)  +dRhoVy*S
      
      F_rho(i,j-1)      =F_rho(i,j-1)   -dRho*S
      F_rho_Et(i,j-1)   =F_rho_Et(i,j-1)-dRhoEt*S
      F_rho_vx(i,j-1)   =F_rho_vx(i,j-1)-dRhoVx*S
      F_rho_vy(i,j-1)   =F_rho_vy(i,j-1)-dRhoVy*S
      
      
      
      end do
      end do     
!  

contains 
    subroutine AUSMPWplus(uL_,vL_,pL_,rhoL_,uR_,vR_,pR_,rhoR_,SNx_,SNy_)
        implicit none
        real*8 :: EPSILON=1E-10
        real*8 :: alpha=3.0/16.0
        real*8 :: uL_,vL_,pL_,rhoL_,uR_,vR_,pR_,rhoR_
        real*8 :: SNx_,SNy_,SNxN,SNyN

        real*8 :: Vm2L,Vm2R,HL,HR,VnL,VnR
        real*8 :: Gm,Gm1,Gm2
        
        !intermediate variables
        real*8 :: aL,aR,aStarL,aStarR,aTildeL,aTildeR,aS
        real*8 :: MaL,MaR,MpL,MnR,PpL,pnR,MpL_,MnR_
        real*8 :: fL,fR
        real*8 :: ps,mHalf,wp
        
        Gm=Gamma
        Gm1=Gm/(Gm-1.0)
        Gm2=2.0*(Gm-1.0)/(Gm+1.0)
        ! magnitude^2 of velocity
        Vm2L=uL_**2+vL_**2 
        Vm2R=uR_**2+vR_**2
        ! enthalpy
        HL=Gm1*pL_/rhoL_+0.5*Vm2L
        HR=Gm1*pR_/rhoR_+0.5*Vm2R
        !sound speed
        aL=sqrt(Gm*pL_/rhoL_)
        aR=sqrt(Gm*pR_/rhoR_)
        ! normalized SD
        SNxN=SNx/sqrt(SNx_**2+SNy_**2)
        SNyN=SNy/sqrt(SNx_**2+SNy_**2)
        ! normal speed on interface
        VnL=(uL_*SNxN+vL_*SNyN)
        VnR=(uR_*SNxN+vR_*SNyN)
        ! critical sound speed
        aStarL=sqrt(Gm2*HL)
        aStarR=sqrt(Gm2*HR)
        ! numerical sound speed
        aTildeL=aStarL**2/max(sqrt(Vm2L),aStarL)
        aTildeR=aStarR**2/max(sqrt(Vm2R),aStarR)
        
        !aTildeL=aStarL**2/max(abs(VnL),aStarL)
        !aTildeR=aStarR**2/max(abs(VnR),aStarR)
        
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
        PpL=0.25*(MaL+1.0)**2*(2.0-MaL)+alpha*MaL*(MaL**2-1.0)**2
        
        endif
        
        if(abs(MaR)>1.0)then
        MnR=0.5*(MaR-abs(MaR))
        !PpL=0.5*(1.0+sign(M))
        PnR=0.5*(MaR-abs(MaR))/(MaR+EPSILON)
        else
        MnR=-0.25*(MaR-1.0)**2
        PnR=0.25*(MaR-1.0)**2*(2.0+MaR)-alpha*MaR*(MaR**2-1.0)**2
        endif
        
        pS=PpL*pL_+PnR*pR_
        mHalf=MpL+MnR
        
        wp=1.0-min(pL_/pR_,pR_/pL_)**3
        
        if(abs(MaL)<1.0)then
        fL=pL_/(pS+EPSILON)-1.0
        else 
        fL=0.0
        endif

        if(abs(MaR)<1.0)then
        fR=pR_/(pS+EPSILON)-1.0
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
        
        dRho  = MpL_*aS*rhoL_      +MnR_*aS*rhoR_
        dRhoVx= MpL_*aS*rhoL_*uL_  +MnR_*aS*rhoR_*uR_  +(PpL*pL+PnR*pR)*SNxN
        dRhoVy= MpL_*aS*rhoL_*vL_  +MnR_*aS*rhoR_*vR_  +(PpL*pL+PnR*pR)*SNyN
        dRhoEt= MpL_*aS*rhoL_*HL   +MnR_*aS*rhoR_*HR

    end

      
      end 
 