      module main
      implicit none
      integer :: i,j
      integer :: iq,jq,ngmax
      real :: EPSILON,EPSILON1,INFINITY,pi
      parameter(iq=405,jq=105,ngmax=max(iq,jq)+5)
      parameter(EPSILON=1.0d-20,EPSILON1=1.0d-10,INFINITY=1.0d20)
      parameter(pi=3.1415926535897932384626)
     
      
      interface
        subroutine LUSGS
        end subroutine LUSGS
        
        subroutine RK2
        end subroutine RK2
        
        subroutine RK4
        end subroutine RK4
        
        subroutine inviscid_fluxes_lax
        end subroutine inviscid_fluxes_lax
        
        subroutine inviscid_fluxes_AUSMPWplus
        end subroutine inviscid_fluxes_AUSMPWplus
     
        subroutine inviscid_fluxes_roe
        end subroutine inviscid_fluxes_roe
        
        subroutine muscl_interpolation
        end subroutine muscl_interpolation
        
        subroutine WENO3_interpolation
        end subroutine WENO3_interpolation
        
        subroutine WENO5_interpolation
        end subroutine WENO5_interpolation
      end interface
      
! function pointer
     procedure(LUSGS),pointer:: pSolver=>NULL()
     procedure(inviscid_fluxes_roe),pointer:: RiemannSolver=>NULL()
     procedure(muscl_interpolation),pointer:: Reconstruct=>NULL()
      logical :: steady1
      
      integer :: &
     & ib,jb &
     &,im,jm,irsolver,iReconstruct
     
      REAL*8 :: &
     & L_ref,V_ref,Rho_ref,T_ref,Ma_ref
      !============================
      ! boundary conditions:
     REAL*8 :: &
     & Vx_inlet,Vy_inlet,p_inlet,T_inlet,Rho_inlet,Mu_inlet,P_out,&
     & KT_inlet,OmegaT_inlet
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !CURRENT STEP SOLUTION variables
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! conservative variables
      REAL*8 :: &
     & rho(iq,jq)&
     &,rho_Et(iq,jq)&
     &,rho_vx(iq,jq),rho_vy(iq,jq)
     
     !  turbulence variables
     REAL*8 :: KT(iq,jq),OmegaT(iq,jq) ! turbulent kinetic energy,Omega
     !  VORTICITY AND SHEAR STRAIN
     REAL*8 :: Vort2(iq,jq),SSM2(iq,jq),DIV(iq,jq),dKdO(iq,jq)
     REAL*8 :: DistW(iq,jq) !distance to the nearest wall
     ! primitive variables
     REAL*8 :: p(iq,jq),T(iq,jq),vx(iq,jq),vy(iq,jq),Ht(iq,jq)
     
     ! SUM Flux
     REAL*8 :: &
     &F_rho(iq,jq),F_rho_Et(iq,jq)&
     &,F_rho_vx(iq,jq),F_rho_vy(iq,jq)
     
     ! TURBULENCE Flux
     REAL*8 :: F_KT(iq,jq),F_OmegaT(iq,jq)
     ! TURBULENCE PRODUCTION AND DISSIPATION ITEM
     REAL*8 :: SPk(iq,jq),SDk(iq,jq)
     REAL*8 :: SPomega(iq,jq),SDomega(iq,jq)
     real*8 :: F2(iq,jq)
     
     ! LAST STEP SOLUTION
     REAL*8 :: &
     &Rho_m1(iq,jq),Rho_Et_m1(iq,jq)&
     &,Rho_vx_m1(iq,jq),Rho_vy_m1(iq,jq)&
     &,KT_m1(iq,jq),OmegaT_m1(iq,jq)
     
  
     ! LAST  2 STEP SOLUTION
     REAL*8 :: &
     &Rho_m2(iq,jq),Rho_Et_m2(iq,jq)&
     &,Rho_vx_m2(iq,jq),Rho_vy_m2(iq,jq)& 
     &,KT_m2(iq,jq),OmegaT_m2(iq,jq)
     
     REAL*8 :: radius_i(iq,jq),radius_j(iq,jq) !SPECTRAL RADIUS
     
     !===================================================
     !      geometry variables
      REAL*8 :: &
     & x(iq,jq),y(iq,jq)&
     &,xc(iq,jq),yc(iq,jq),step(iq,jq)&
     &,Vcell(iq,jq)

!     
      REAL*8 :: &
     &            Sn_X_i(iq,jq),Sn_Y_i(iq,jq)&
     &           ,Sn_X_j(iq,jq),Sn_Y_j(iq,jq)&
     &           ,S_i(iq,jq),S_j(iq,jq)&
     &           ,L_Cell_x(iq,jq),L_Cell_y(iq,jq)
!   
      REAL*8 :: &
!
     & vib(iq,jq),vjb(iq,jq)

     
     
     ! interface variables      
      REAL*8 :: &
!      
     &uil (iq,jq),uir(iq,jq),&
     &vil (iq,jq),vir(iq,jq),&
     &pil (iq,jq),pir(iq,jq),&
     &ril (iq,jq),rir(iq,jq),&

     &ujl (iq,jq),ujr(iq,jq),&
     &vjl (iq,jq),vjr(iq,jq),&
     &pjl (iq,jq),pjr(iq,jq),&
     &rjl (iq,jq),rjr(iq,jq)

     REAL*8 :: &
     & KTil(iq,jq),KTir(iq,jq),&
     & OmegaTil(iq,jq),OmegaTir(iq,jq),&
     
     & KTjl(iq,jq),KTjr(iq,jq),&
     & OmegaTjl(iq,jq),OmegaTjr(iq,jq)
     
    integer :: TimeMarch ! select time marching method 0:LUSGS; 1:RK2; 2:RK4
    integer :: &
     &nmax,n,is_restart,&
     &ioup,iinp,icgl,iimplicit,idgstart&
     &,nts1p,nunsloop,ninnerloop,is_print,iprint,istart,nb

      REAL*8 :: &
!
     & cfl,ttime,tstep,epsall
!     
      REAL*8 :: &
!
     & cp,Gamma,Gamma2,Gamma1,rcp,R_air,Gamma3,cv,Pr_L,Pr_T,Pr_E
!     
!
      REAL*8 :: &
!
     &Mu_L(iq,jq),Mu_T(iq,jq),Mu_E(iq,jq)


	  !REAL*8 :: &
   ! & Vx_Inlet,Vy_Inlet,T_Inlet,p_Inlet,p_Out
    
     

    !================ definations of some DIFF functions =============
    !       diff on the &&  LEFT/SOUTH  && interface
    contains

    function dx_i(Var,ii,jj)
    implicit none
    integer :: ii,jj
    real*8 :: dx_i,DI,DJ
    real*8 :: Var(iq,jq)
    
    DI=Var(ii,jj)-Var(ii-1,jj)
    DJ=0.25*(Var(ii-1,jj+1)+Var(ii,jj+1)-Var(ii-1,jj-1)-Var(ii,jj-1))
    
    dx_i=( DI*(y(ii,jj+1)-y(ii,jj)) &
         &-DJ*(yc(ii,jj)-yc(ii-1,jj))) &
         & / &
         &((xc(ii,jj)-xc(ii-1,jj))*(y(ii,jj+1)-y(ii,jj)) &
         &-(yc(ii,jj)-yc(ii-1,jj))*(x(ii,jj+1)-x(ii,jj)) )
       
    end function
    
    function dy_i(Var,ii,jj)
    implicit none
    integer :: ii,jj
    real*8 :: dy_i,DI,DJ
    real*8 :: Var(iq,jq)
    
    DI=Var(ii,jj)-Var(ii-1,jj)
    DJ=0.25*(Var(ii-1,jj+1)+Var(ii,jj+1)-Var(ii-1,jj-1)-Var(ii,jj-1))
    
    dy_i=( DJ*(xc(ii,jj+1)-xc(ii,jj)) &
         &-DI*(y(ii,jj+1)-y(ii,jj))) &
         & / &
         &((xc(ii,jj)-xc(ii-1,jj))*(y(ii,jj+1)-y(ii,jj)) &
         &-(yc(ii,jj)-yc(ii-1,jj))*(x(ii,jj+1)-x(ii,jj)) )
       
    end function
    
    function dx_j(Var,ii,jj)
    implicit none
    integer :: ii,jj
    real*8 :: dx_j,DI,DJ
    real*8 :: Var(iq,jq)
    
    DJ=Var(ii,jj)-Var(ii,jj-1)
    DI=0.25*(Var(ii+1,jj)+Var(ii+1,jj-1)-Var(ii-1,jj)-Var(ii-1,jj-1))
    
    dx_j=( DI*(yc(ii,jj)-yc(ii,jj-1)) &
         &-DJ*(y(ii+1,jj)-y(ii,jj))) &
         & / &
         &((x(ii+1,jj)-x(ii,jj))*(yc(ii,jj)-yc(ii,jj-1)) &
         &-(y(ii+1,jj)-y(ii,jj))*(xc(ii,jj)-xc(ii,jj-1)) )
       
    end function
    
    function dy_j(Var,ii,jj)
    implicit none
    integer :: ii,jj
    real*8 :: dy_j,DI,DJ
    real*8 :: Var(iq,jq)
    
    DJ=Var(ii,jj)-Var(ii,jj-1)
    DI=0.25*(Var(ii+1,jj)+Var(ii+1,jj-1)-Var(ii-1,jj)-Var(ii-1,jj-1))
    
    dy_j=( DJ*(x(ii+1,jj)-x(ii,jj)) &
         &-DI*(xc(ii,jj)-xc(ii,jj-1))) &
         & / &
         &((x(ii+1,jj)-x(ii,jj))*(yc(ii,jj)-yc(ii,jj-1)) &
         &-(y(ii+1,jj)-y(ii,jj))*(xc(ii,jj)-xc(ii,jj-1)) )
       
    end function
    ! get the directional derivative by Green method
    ! interface flux is to be used
    function dx(FW,FE,FN,FS,ii,jj)
    implicit none
    real*8 :: FW(iq,jq),FE(iq,jq),FN(iq,jq),FS(iq,jq)
    real*8 :: dx
    integer :: ii,jj
    dx=(FE(i+1,j)*SN_x_i(i+1,j)-FW(i,j)*SN_x_i(i,j)&
    &     +FN(i,j+1)*SN_x_j(i,j+1)-FS(i,j)*SN_x_j(i,j))&
            &/Vcell(i,j)
    end function
    
    function dy(FW,FE,FN,FS,ii,jj)
    implicit none
    real*8 :: FW(iq,jq),FE(iq,jq),FN(iq,jq),FS(iq,jq)
    real*8 :: dy
    integer :: ii,jj
    dy=(FE(i+1,j)*SN_Y_i(i+1,j)-FW(i,j)*SN_Y_i(i,j)&
    &     +FN(i,j+1)*SN_Y_j(i,j+1)-FS(i,j)*SN_Y_j(i,j))&
            &/Vcell(i,j)
    end function
    
      end module main
      
!
