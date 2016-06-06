subroutine viscous_Flux
    use main 
    ! compute the viscous flux 
    ! add to SumFlux()
    real*8 :: sav1,sav2
    real*8 :: Tau_xx,Tau_xy,Tau_yy !stress
    real*8 :: LFlux_rho_vx,LFlux_rho_vy,LFlux_Et(iq,jq)
    real*8 :: k_Fourier=1.0 ! rate of heat conduction
    real*8 :: DTdx_I(iq,jq),DTdy_I(iq,jq),DTdx_J(iq,jq),DTdy_J(iq,jq)
      
    
    ! -   i direction
      do j=jb,jm-1
      do i=ib,im
      
      
      !Stokes hypothesis: lambda=2*mu/3
      !mu on the interface is given by linear average: 0.5*(Mu_L(i,j)+Mu_L(i+1,j))
      Tau_xx=(Mu_E(i,j)+Mu_E(i-1,j))/3.0*(2.0*dx_i(vx,i,j)-dy_i(vy,i,j))
      !write(*,*) dx_i(vx,i,j),dy_i(vx,i,j)
      Tau_xy=0.5*(Mu_E(i,j)+Mu_E(i-1,j))*(dx_i(vy,i,j)+dy_i(vx,i,j))
      Tau_yy=(Mu_E(i,j)+Mu_E(i-1,j))/3.0*(2.0*dy_i(vy,i,j)-dx_i(vx,i,j))
!     surface area vectors
!
      sav1=Sn_X_i(i,j) !dy
      sav2=Sn_Y_i(i,j) !dx
!
      k_Fourier=Mu_L(i,j)*cp/Pr_L+Mu_T(i,j)*cp/Pr_L
      
      LFlux_rho_vx=Tau_xx*sav1+Tau_xy*sav2
      LFlux_rho_vy=Tau_xy*sav1+Tau_yy*sav2
      !DTdx_I(i,j)=dx_i(T,i,j)
      !DTdy_I(i,j)=dy_i(T,i,j)
      LFlux_Et(i,j)=(uil(i,j)*Tau_xx+vil(i,j)*Tau_xy+k_Fourier*dx_i(T,i,j))*sav1 &
      & +(uil(i,j)*Tau_xy+vil(i,j)*Tau_yy+k_Fourier*dy_i(T,i,j))*sav2
      
      !SumFlux_rho(i,j)=SumFlux_rho(i,j)  ! no Flux to add
      SumFlux_rho_Et(i,j)=SumFlux_rho_Et(i,j)-LFlux_Et(i,j)
      SumFlux_rho_vx(i,j)=SumFlux_rho_vx(i,j)-LFlux_rho_vx
      SumFlux_rho_vy(i,j)=SumFlux_rho_vy(i,j)-LFlux_rho_vy
     
      !SumFlux_rho(i-1,j)=SumFlux_rho(i-1,j) ! no Flux to add
      !RFlux(i-1,j)=LFlux(i,j)
      SumFlux_rho_Et(i-1,j)=SumFlux_rho_Et(i-1,j)+LFlux_Et(i,j) 
      SumFlux_rho_vx(i-1,j)=SumFlux_rho_vx(i-1,j)+LFlux_rho_vx
      SumFlux_rho_vy(i-1,j)=SumFlux_rho_vy(i-1,j)+LFlux_rho_vy
      end do
      end do

!
! -   j direction
!
      
!
      do j=jb,jm
      do i=ib,im-1
      !Stokes hypothesis: lambda=2*mu/3
      !mu on the interface is given by linear average:
      !     0.5*(Mu_L(i,j)+Mu_L(i-1,j))
      Tau_xx=(Mu_E(i,j)+Mu_E(i,j-1))/3.0*(2.0*dx_j(vx,i,j)-dy_j(vy,i,j))
      Tau_xy=0.5*(Mu_E(i,j)+Mu_E(i,j-1))*(dx_j(vy,i,j)+dy_j(vx,i,j))
     
      Tau_yy=(Mu_E(i,j)+Mu_E(i,j-1))/3.0*(2.0*dy_j(vy,i,j)-dx_j(vx,i,j))
      !write(*,*) dy_j(vx,i,j),dx_j(vy,i,j),Tau_xy
!     surface area vectors
!
      sav1=Sn_X_j(i,j)
      sav2=Sn_Y_j(i,j)      
      
      
      k_Fourier=Mu_L(i,j)*cp/Pr_L
      
      LFlux_rho_vx=Tau_xx*sav1+Tau_xy*sav2
      LFlux_rho_vy=Tau_xy*sav1+Tau_yy*sav2
      !DTdx_J(i,j)=dx_J(T,i,j)
      !DTdy_J(i,j)=dy_J(T,i,j)
      LFlux_Et(i,j)=(ujl(i,j)*Tau_xx+vjl(i,j)*Tau_xy+k_Fourier*dx_j(T,i,j))*sav1 &
      & +(ujl(i,j)*Tau_xy+vjl(i,j)*Tau_yy+k_Fourier*dy_j(T,i,j))*sav2
      
      
      !SumFlux_rho(i,j)=SumFlux_rho(i,j) ! no Flux to add
      SumFlux_rho_Et(i,j)=SumFlux_rho_Et(i,j)-LFlux_Et(i,j)
      SumFlux_rho_vx(i,j)=SumFlux_rho_vx(i,j)-LFlux_rho_vx
      SumFlux_rho_vy(i,j)=SumFlux_rho_vy(i,j)-LFlux_rho_vy
      
      !SumFlux_rho(i,j-1)=SumFlux_rho(i,j-1)    ! no Flux to add
      SumFlux_rho_Et(i,j-1)=SumFlux_rho_Et(i,j-1)+LFlux_Et(i,j)
      SumFlux_rho_vx(i,j-1)=SumFlux_rho_vx(i,j-1)+LFlux_rho_vx
      SumFlux_rho_vy(i,j-1)=SumFlux_rho_vy(i,j-1)+LFlux_rho_vy
      end do
      end do
    
end subroutine

