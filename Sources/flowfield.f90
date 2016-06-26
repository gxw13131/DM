!====================================================================================
    subroutine DistanceWall
    ! compute the distance to the nearest wall
    ! for SST turbulence model
    use main
    implicit none
    real*8 :: xx,yy,temp
    integer :: ii
      do j=jb,jm-1
      do i=ib,im-1
       xx=x(i,j)
       yy=y(i,j)
       distW(i,j)=1e10
           do ii=ib,im-1
           ! j=jb is set to the wall by default
           temp=sqrt((xx-y(ii,jb))**2+(yy-y(ii,jb))**2)  
           if (temp.lt.distW(i,j)) then
            distW(i,j)=temp
           end if
          end do
      end do
      end do
    end subroutine
    
      subroutine radius
      !compute the spectrum radius
      !MIND THAT: !!!the result is coupled with J
      use main  
!     compute spectral radii
!     ----------------------
!
      real*8 :: Aco2,Aco,Arx,Ary,Vibc,Sii,Armx,Velpro,Acoa,Vjbc,Sjj,Army,Ste
      real*8 :: ttg
      do j=jb,jm-1
      do i=ib,im-1
          
           aco2=abs(Gamma*p(i,j)/rho(i,j))
           aco=sqrt(aco2) 
          
!       i direction
          arx=0.5*(Sn_X_i(i,j)+Sn_X_i(i+1,j))
          ary=0.5*(Sn_Y_i(i,j)+Sn_Y_i(i+1,j))
          vibc=0.5*(vib(i,j)+vib(i+1,j))
          sii=arx*arx+ary*ary
          armx=sqrt(sii)
          velpro=vx(i,j)*arx+vy(i,j)*ary+vibc
          acoa=armx*aco
          radius_i(i,j)=abs(velpro)+acoa&
          &+2.0*sii*max(1.3333,Gamma/Pr_L)*Mu_E(i,j)/rho(i,j)/Vcell(i,j)

!       j direction
          arx=0.5*(Sn_X_j(i,j)+Sn_X_j(i,j+1))
          ary=0.5*(Sn_Y_j(i,j)+Sn_Y_j(i,j+1))
          vjbc=0.5*(vjb(i,j)+vjb(i,j+1))
          sjj=arx*arx+ary*ary
          army=sqrt(sjj)
          velpro=vx(i,j)*arx+vy(i,j)*ary+vjbc
          acoa=army*aco
          radius_j(i,j)=abs(velpro)+acoa&
          &+2.0*sjj*max(1.3333,Gamma/Pr_L)*Mu_E(i,j)/rho(i,j)/Vcell(i,j)
!
      end do
      end do
      
      end subroutine     
      
      subroutine time_step
      !===================
      
      use main  
      
!
!     compute spectral radii
!     ----------------------
!
      real*8 :: ttg
      do j=jb,jm-1
      do i=ib,im-1
      
          step(i,j)=cfl*Vcell(i,j)/(radius_i(i,j)+radius_j(i,j))
!
      end do
      end do
!
      if(icgl.eq.1) then
!
      ttg=1000.
      
     
      do j=jb,jm-1
       do i=ib,im-1
      if(ttg.gt.step(i,j)) ttg=step(i,j)
!       ttg is set to the minimum of step(i,j)
      end do
      end do
!
      
      do j=jb,jm-1
      do i=ib,im-1
!       all step(i,j) is set to ttg(the minimum of step(i,j) 
!       to guarantee CFL condition is not violated
      step(i,j)=ttg
!
      end do
      end do
!
      end if
      
      return
      end
!====================================================================================
      
     subroutine interpolation_to_interface
     !===================================
      use main
      
      real*8 :: AlsInv
      
      
      
      ! interpolate the primitive variables to cell interfaces
      
      ! iface
     
      do j=jb,jm-1
      do i=ib,im
      ! linear average to interpolate,
      ! distance weight is used to deal with non-homogeneous grid  
          alsinv=1./(L_Cell_x(i,j)+L_Cell_x(i-1,j))
          uil (i,j)=(L_Cell_x(i,j)*vx (i-1,j)+L_Cell_x(i-1,j)*vx (i,j))&
     &                  *alsinv
          vil (i,j)=(L_Cell_x(i,j)*vy (i-1,j)+L_Cell_x(i-1,j)*vy (i,j))&
     &                  *alsinv
          
          pil (i,j)=(L_Cell_x(i,j)*p  (i-1,j)+L_Cell_x(i-1,j)*p  (i,j))&
     &                  *alsinv
          ril (i,j)=(L_Cell_x(i,j)*rho(i-1,j)+L_Cell_x(i-1,j)*rho(i,j))&
     &                  *alsinv
     
      end do
      end do
!
!     j face
!
     
      do j=jb,jm
      do i=ib,im-1
          alsinv=1./(L_Cell_y(i,j)+L_Cell_y(i,j-1))
          ujl (i,j)=(L_Cell_y(i,j)*vx (i,j-1)+L_Cell_y(i,j-1)*vx (i,j))&
     &                  *alsinv
          vjl (i,j)=(L_Cell_y(i,j)*vy (i,j-1)+L_Cell_y(i,j-1)*vy (i,j))&
     &                  *alsinv
          
          pjl (i,j)=(L_Cell_y(i,j)*p  (i,j-1)+L_Cell_y(i,j-1)*p  (i,j))&
     &                  *alsinv
          rjl (i,j)=(L_Cell_y(i,j)*rho(i,j-1)+L_Cell_y(i,j-1)*rho(i,j))&
     &                  *alsinv
                 
      end do
      end do
      
      
      return
      end
      
!=====================================================================================      
      subroutine muscl_interpolation 
     !===============================
      use main 
      
      real*8 :: Aco2,Aco,Arx,Ary,Vibc,Sii,Armx,Velpro,Acoa,Eox1,Vjbc,Sjj,Army,Eoy1,Ste
      real*8 :: ttg, Alsinv,epsm
      real*8 :: DR,DL,UL,UR,VL,VR,PL,PR,RL,RR,VamL,VamR,HL,HR,AcL,AcR
      real*8 :: RRoRL,RRoRLp1
      real*8 :: Rm,Um,Vm,Hm,Vm2,Am2,Am
      real*8 :: sav1,sav2,sav,sav1N,sav2N
      real*8 :: DUnormal,Drou,Dp,DRU,DRV,DRE
      real*8 :: sos,sos1,QNo,QN
      
      ! MUSCL interpolation at cell interface except boundaries
     
     
   !   alimiter(var1,var2)=2./3.*var1+1./3.*var2
      epsm=EPSILON
      
    ! store cell interface value to facilitate the limiting process   
      
      
      do j=jb,jm-1
      do i=ib,im
      
         uir(i,j)=uil(i,j)
         vir(i,j)=vil(i,j)
         pir(i,j)=pil(i,j)
         rir(i,j)=ril(i,j)
      
      end do
      end do
      
     
      do j=jb,jm
       do i=ib,im-1
      
         ujr(i,j)=ujl(i,j)
         vjr(i,j)=vjl(i,j)
         pjr(i,j)=pjl(i,j)
         rjr(i,j)=rjl(i,j)
         
      end do
      end do
      
      
!
!
!     compute limited cell-interface value
!
      
      do j=jb,jm-1
      do i=ib,im-1
      
      !limiting u
      
      dr=(uil(i+1,j)-  vx(i,j))*2.
      dl=(vx(i,j)    -uir(i,j))*2. !uir(i,j)=uil(i,j)
      
      uir(i,j)  =vx(i,j)-0.5*alimiter(dl,dr)
      uil(i+1,j)=vx(i,j)+0.5*alimiter(dr,dl)
      
      dr=(ujl(i,j+1)-  vx(i,j))*2.
      dl=(vx(i,j)    -ujr(i,j))*2.
      
      ujr(i,j)  =vx(i,j)-0.5*alimiter(dl,dr)
      ujl(i,j+1)=vx(i,j)+0.5*alimiter(dr,dl)
      
      
      !limiting v
      
      dr=(vil(i+1,j)-  vy(i,j))*2.
      dl=(vy(i,j)    -vir(i,j))*2.
      
      vir(i,j)  =vy(i,j)-0.5*alimiter(dl,dr)
      vil(i+1,j)=vy(i,j)+0.5*alimiter(dr,dl)
      
      dr=(vjl(i,j+1)-  vy(i,j))*2.
      dl=(vy(i,j)    -vjr(i,j))*2.
      
      vjr(i,j)  =vy(i,j)-0.5*alimiter(dl,dr)
      vjl(i,j+1)=vy(i,j)+0.5*alimiter(dr,dl)
      
      !limiting p
          
      dr=(pil(i+1,j)-  p(i,j))*2.
      dl=(p(i,j)    -pir(i,j))*2.
      
      pir(i,j)  =p(i,j)-0.5*alimiter(dl,dr)
      pil(i+1,j)=p(i,j)+0.5*alimiter(dr,dl)
      
      dr=(pjl(i,j+1)-  p(i,j))*2.
      dl=(p(i,j)    -pjr(i,j))*2.
      
      pjr(i,j)  =p(i,j)-0.5*alimiter(dl,dr)
      pjl(i,j+1)=p(i,j)+0.5*alimiter(dr,dl)
      
      !limiting rou
      
      dr=(ril(i+1,j)-  rho(i,j))*2.
      dl=(rho(i,j)    -rir(i,j))*2.
      
      rir(i,j)  =rho(i,j)-0.5*alimiter(dl,dr)
      ril(i+1,j)=rho(i,j)+0.5*alimiter(dr,dl)
      
      dr=(rjl(i,j+1)-  rho(i,j))*2.
      dl=(rho(i,j)    -rjr(i,j))*2.
      
      rjr(i,j)  =rho(i,j)-0.5*alimiter(dl,dr)
      rjl(i,j+1)=rho(i,j)+0.5*alimiter(dr,dl)
      
      end do
      end do
     
contains
      function  alimiter(var1,var2)
      real*8 :: var1,var2
      
      real*8 :: alimiter
       alimiter=(var1*(var2*var2+2.*epsm*epsm)+&
     &                     var2*(2.*var1*var1+epsm*epsm))/&
     &            (2.*var1*var1-var1*var2+2.*var2*var2+3.*epsm*epsm) &
     &            *0.5*abs(sign(1.0,var1)+sign(1.0,var2)) 
     
      end function alimiter
      
      end
      
      subroutine WENO5_interpolation
      use main
      implicit none
      do j=jb,jm-1
      do i=ib,im
         
         uil(i,j)=WENO5(vx(i-3,j),vx(i-2,j),vx(i-1,j),vx(i,j),vx(i+1,j))
         uir(i,j)=WENO5(vx(i+2,j),vx(i+1,j),vx(i,j),vx(i-1,j),vx(i-2,j))
         
         vil(i,j)=WENO5(vy(i-3,j),vy(i-2,j),vy(i-1,j),vy(i,j),vy(i+1,j))
         vir(i,j)=WENO5(vy(i+2,j),vy(i+1,j),vy(i,j),vy(i-1,j),vy(i-2,j))
         
         pil(i,j)=WENO5(p(i-3,j),p(i-2,j),p(i-1,j),p(i,j),p(i+1,j))
         pir(i,j)=WENO5(p(i+2,j),p(i+1,j),p(i,j),p(i-1,j),p(i-2,j))
         
         ril(i,j)=WENO5(Rho(i-3,j),Rho(i-2,j),Rho(i-1,j),Rho(i,j),Rho(i+1,j))
         rir(i,j)=WENO5(Rho(i+2,j),Rho(i+1,j),Rho(i,j),Rho(i-1,j),Rho(i-2,j))
      
      end do
      end do
      
     
      do j=jb,jm
       do i=ib,im-1
         ujl(i,j)=WENO5(vx(i,j-3),vx(i,j-2),vx(i,j-1),vx(i,j),vx(i,j+1))
         ujr(i,j)=WENO5(vx(i,j+2),vx(i,j+1),vx(i,j),vx(i,j-1),vx(i,j-2))
         
         vjl(i,j)=WENO5(vy(i,j-3),vy(i,j-2),vy(i,j-1),vy(i,j),vy(i,j+1))
         vjr(i,j)=WENO5(vy(i,j+2),vy(i,j+1),vy(i,j),vy(i,j-1),vy(i,j-2))
         
         pjl(i,j)=WENO5(p(i,j-3),p(i,j-2),p(i,j-1),p(i,j),p(i,j+1))
         pjr(i,j)=WENO5(p(i,j+2),p(i,j+1),p(i,j),p(i,j-1),p(i,j-2))
         
         rjl(i,j)=WENO5(Rho(i,j-3),Rho(i,j-2),Rho(i,j-1),Rho(i,j),Rho(i,j+1))
         rjr(i,j)=WENO5(Rho(i,j+2),Rho(i,j+1),Rho(i,j),Rho(i,j-1),Rho(i,j-2))
         
         
      end do
      end do
      
contains
    function WENO5(v_3,v_2,v_1,v0,v1)
        ! 5 order WENO 
        ! use point i-3,i-2,i-1,i,i+1 to reconstruct interface i-1/2
        implicit none
        real*8 :: v_3,v_2,v_1,v0,v1
        real*8 :: WENO5
	    real*8 :: aw0 = 0.1, aw1 = 0.6, aw2 = 0.3
	
	    real*8 :: fS0, fS1, fS2
	    real*8 :: w0, w1, w2, w0_, w1_, w2_
	    real*8 :: vL0, vL1, vL2
	    real*8 :: EPS = 1e-7
	    !smooth factor
        !1.083333=13/12
	    fS0 = 1.08333*(v_3 - 2.0*v_2 + v_1)**2  + 0.25*(1.0*v_3 - 4.0*v_2 + 3.0*v_1)**2
	    fS1 = 1.08333*(v_2 - 2.0*v_1 + v0)**2   + 0.25*(-1.0*v_2  + 1.0*v0)**2
	    fS2 = 1.08333*(v_1 - 2.0*v0 + v1)**2    + 0.25*(-3.0*v_1 + 4.0*v0 - 1.0*v1)**2

	    ! weight factor
	    w0_ = aw0 / (EPS + fS0)**2
	    w1_ = aw1 / (EPS + fS1)**2
	    w2_ = aw2 / (EPS + fS2)**2
	    ! normalize
	    w0 = w0_ / (w0_ + w1_ + w2_)
	    w1 = w1_ / (w0_ + w1_ + w2_)
	    w2 = w2_ / (w0_ + w1_ + w2_)
	    vL0 = (2.0*v_3	- 7.0*v_2	+ 11.0*v_1) / 6.0
	    vL1 = (-1.0*v_2 + 5.0*v_1	+ 2.0*v0) / 6.0
	    vL2 = (2.0*v_1	+ 5.0*v0	- 1.0*v1) / 6.0
        ! return
	    WENO5 = vL0*w0 + vL1*w1 + vL2*w2
    end function WENO5   
      end subroutine
!====================================================================================    
      
      subroutine WENO3_interpolation
      use main
      implicit none
      do j=jb,jm-1
      do i=ib,im
         
         uil(i,j)=WENO3(vx(i-2,j),vx(i-1,j),vx(i,j))
         uir(i,j)=WENO3(vx(i+1,j),vx(i,j),vx(i-1,j))
         
         vil(i,j)=WENO3(vy(i-2,j),vy(i-1,j),vy(i,j))
         vir(i,j)=WENO3(vy(i+1,j),vy(i,j),vy(i-1,j))
         
         pil(i,j)=WENO3(p(i-2,j),p(i-1,j),p(i,j))
         pir(i,j)=WENO3(p(i+1,j),p(i,j),p(i-1,j))
         
         ril(i,j)=WENO3(Rho(i-2,j),Rho(i-1,j),Rho(i,j))
         rir(i,j)=WENO3(Rho(i+1,j),Rho(i,j),Rho(i-1,j))
      
      end do
      end do
      
     
      do j=jb,jm
       do i=ib,im-1
         ujl(i,j)=WENO3(vx(i,j-2),vx(i,j-1),vx(i,j))
         ujr(i,j)=WENO3(vx(i,j+1),vx(i,j),vx(i,j-1))
         
         vjl(i,j)=WENO3(vy(i,j-2),vy(i,j-1),vy(i,j))
         vjr(i,j)=WENO3(vy(i,j+1),vy(i,j),vy(i,j-1))
         
         pjl(i,j)=WENO3(p(i,j-2),p(i,j-1),p(i,j))
         pjr(i,j)=WENO3(p(i,j+1),p(i,j),p(i,j-1))
         
         rjl(i,j)=WENO3(Rho(i,j-2),Rho(i,j-1),Rho(i,j))
         rjr(i,j)=WENO3(Rho(i,j+1),Rho(i,j),Rho(i,j-1))        
      end do
      end do
      
contains
    function WENO3(v_2, v_1, v0)
        ! 3 order WENO 
        ! use point i-2,i-1,i to reconstruct interface i-1/2
        implicit none
        real*8 :: v_2,v_1,v0
        real*8 :: WENO3
    
	    real*8 :: aw0 = 1.0/3.0, aw1 = 2.0/3.0
	    real*8 :: fS0, fS1
	    real*8 :: w0, w1, w0_, w1_
	    real*8 :: vL0, vL1
	    real*8 :: EPS = 1e-8
	    !smooth factor
	    fS0 = (v_1-v_2)**2
	    fS1 = (v_1-v0 )**2 
	
	    ! weight factor
	    w0_ = aw0 / (EPS + fS0)**2
	    w1_ = aw1 / (EPS + fS1)**2
	
	    ! normalize
	    w0 = w0_ / (w0_ + w1_ )
	    w1 = w1_ / (w0_ + w1_ )
	
	    vL0 = 1.5*v_1-0.5*v_2
	    vL1 = 0.5*v_1+0.5*v0
	    ! return 
	    WENO3 = vL0*w0 + vL1*w1
    end function WENO3
    
 end subroutine WENO3_interpolation
 
      subroutine inviscid_fluxes_roe
!     =========================
      
      use main
      
      real*8 :: Aco2,Aco,Arx,Ary,Vibc,Sii,Armx,Velpro,Acoa,Eox1,Vjbc,Sjj,Army,Eoy1,Ste
      real*8 :: ttg, Alsinv,epsm
      real*8 :: DR,DL,UL,UR,VL,VR,PL,PR,RL,RR,VamL,VamR,HL,HR,AcL,AcR
      real*8 :: RRoRL,RRoRLp1
      real*8 :: Rm,Um,Vm,Hm,Vm2,Am2,Am
      real*8 :: sav1,sav2,sav,sav1N,sav2N
      real*8 :: DUnormal,Drou,Dp,DRU,DRV,DRE
      real*8 :: sos,sosI,qn,qn0
       
      real*8 :: DU,DV
      real*8 :: aqn,amn,am0,aam0,omam0,sosIP,sosp,coef,Dpc1,Dpc2,Dqnc1,Dqnc2
      real*8 :: Adr,AdrU,AdrE,ULnormal,Urnormal,RULnormal,RURnormal
      real*8 :: Droi,Dxi,Drei,Droj,Dxj,Dyj,Drej
      real*8 :: EL,ER
      real*8 :: Adrv,Dyi
      
      epsm=0.1
      
      do i=ib-1,im+1
!
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
!
      do i=ib,im
      
      
      
      
!    
!     left&right states
!
      
      
      ul=uil(i,j)
      ur=uir(i,j) 
     
      vl=vil(i,j)
      vr=vir(i,j)
     
     
      pl=max(pil(i,j),EPSILON)
      pr=max(pir(i,j),EPSILON)
     
      rl=max(ril(i,j),EPSILON)
      rr=max(rir(i,j),EPSILON)
           
      vaml=ul*ul+vl*vl
      vamr=ur*ur+vr*vr
      hl=pl*Gamma/(rl*Gamma1)+0.5*vaml
      hr=pr*Gamma/(rr*Gamma1)+0.5*vamr 
      acl=sqrt(Gamma*pl/rl)
      acr=sqrt(Gamma*pr/rr)
!
!     Roe average
!
      
      rrorl=sqrt(rr/rl)
      rrorlp1=1.+rrorl
      
      rm=sqrt(rr*rl)
      um=(ul+ur*rrorl)/rrorlp1
      vm=(vl+vr*rrorl)/rrorlp1
     
      hm=(hl+hr*rrorl)/rrorlp1
      vm2=um*um+vm*vm
      am2=Gamma1*abs(hm-0.5*vm2)
      am=sqrt(am2)
      
                  
! 
!     
!     surface area vectors
!
      sav1=Sn_X_i(i,j)
      sav2=Sn_Y_i(i,j)
      
      sav=S_i(i,j)+EPSILON
      sav1n=sav1/sav
      sav2n=sav2/sav
      
!
!     engenvalues
!      
         
           du=(ur-ul)
           dv=(vr-vl)
           
           dunormal=(sav1n*du+sav2n*dv)
           drou=rr-rl
           dp=pr-pl
           dru=rm*du+um*drou
           drv=rm*dv+vm*drou
          
           dre=dp/Gamma1+0.5*vm2*drou+rm*um*du+rm*vm*dv
           
           sos=am
           sosi=1./sos
      
           qn0=sav1n*um+sav2n*vm
           qn=qn0+vib(i,j)/sav
           aqn=abs(qn)
           amn=qn*sosi
           am0=min(abs(amn),1.)*sign(1.,amn)
           if(abs(amn-1.)<epsm) am0=(-(1.-epsm)**2+2.*(1.+epsm)*amn-amn*amn)/(4.*epsm)
           if(abs(amn+1.)<epsm) am0=((1.-epsm)**2+2.*(1.+epsm)*amn+amn*amn)/(4.*epsm)
           aam0=abs(am0)
          ! if(aam0<epsm) aam0=(epsm*epsm+aam0*aam0)/(2.*epsm)
          ! if(aqn<epsm*sos) aqn=(epsm*epsm*sos*sos+aqn*aqn)/(2.*epsm*sos)
           omam0=1.-aam0
           
           sosip=sosi
           sosp=sos 
           
           
!
!     inviscid flux
!      
      coef=S_i(i,j)
      dpc1=sosip*omam0*dp
      dpc2=am0*dp
      dqnc1=rm*am0*dunormal
      dqnc2=rm*sosp*omam0*dunormal
      
      
      adr =coef*(              dpc1              +dqnc1+aqn*drou)
      adru=coef*(sav1n*dpc2+um*dpc1+sav1n*dqnc2+um*dqnc1+aqn*dru)
      adrv=coef*(sav2n*dpc2+vm*dpc1+sav2n*dqnc2+vm*dqnc1+aqn*drv)  
      adre=coef*(  qn0*dpc2+hm*dpc1+  qn0*dqnc2+hm*dqnc1+aqn*dre)  
       
      
      ulnormal=(sav1*ul+sav2*vl+vib(i,j))
      urnormal=(sav1*ur+sav2*vr+vib(i,j))
      rulnormal=rl*ulnormal
      rurnormal=rr*urnormal
      
!
!     flux terms
!
      droi=-0.5*(rulnormal+rurnormal-adr)
      dxi=-0.5*(rulnormal*ul+rurnormal*ur&
     &            +sav1*(pl+pr)-adru) 
      dyi=-0.5*(rulnormal*vl+rurnormal*vr&
     &            +sav2*(pl+pr)-adrv)  
      drei=-0.5*(rulnormal*hl+rurnormal*hr&
     &            -vib(i,j)*(pl+pr)-adre)    
    
      
     
      F_rho(i,j)=F_rho(i,j)-droi
      F_rho_Et(i,j)=F_rho_Et(i,j)-drei
      F_rho_vx(i,j)=F_rho_vx(i,j)-dxi
      F_rho_vy(i,j)=F_rho_vy(i,j)-dyi
     
      F_rho(i-1,j)=F_rho(i-1,j)+droi
      F_rho_Et(i-1,j)=F_rho_Et(i-1,j)+drei
      F_rho_vx(i-1,j)=F_rho_vx(i-1,j)+dxi
      F_rho_vy(i-1,j)=F_rho_vy(i-1,j)+dyi
      
      
      
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
      ul=ujl(i,j)
      ur=ujr(i,j) 
     
      vl=vjl(i,j)
      vr=vjr(i,j)
     
     
      pl=max(pjl(i,j),EPSILON)
      pr=max(pjr(i,j),EPSILON)
     
      rl=max(rjl(i,j),EPSILON)
      rr=max(rjr(i,j),EPSILON)
          
      vaml=ul*ul+vl*vl
      vamr=ur*ur+vr*vr
      hl=pl*Gamma/(rl*Gamma1)+0.5*vaml
      hr=pr*Gamma/(rr*Gamma1)+0.5*vamr 
      acl=sqrt(Gamma*pl/rl)
      acr=sqrt(Gamma*pr/rr)
!
!     Roe average
!
      
      rrorl=sqrt(rr/rl)
      rrorlp1=1.+rrorl
      
      rm=sqrt(rr*rl)
      um=(ul+ur*rrorl)/rrorlp1
      vm=(vl+vr*rrorl)/rrorlp1
      hm=(hl+hr*rrorl)/rrorlp1
      vm2=um*um+vm*vm
      am2=Gamma1*abs(hm-0.5*vm2)
      am=sqrt(am2)
      
                  
! 
!     
!     surface area vectors
!
      sav1=Sn_X_j(i,j)
      sav2=Sn_Y_j(i,j)
      sav=S_j(i,j)+EPSILON
      sav1n=sav1/sav
      sav2n=sav2/sav
          
!
!     engenvalues
!      
           du=(ur-ul)
           dv=(vr-vl)
           dunormal=(sav1n*du+sav2n*dv)
           drou=rr-rl
           dp=pr-pl
           dru=rm*du+um*drou
           drv=rm*dv+vm*drou
           dre=dp/Gamma1+0.5*vm2*drou+rm*um*du+rm*vm*dv 
           
           sos=am
           sosi=1./sos
       
           qn0=sav1n*um+sav2n*vm
           qn=qn0+vjb(i,j)/sav
           aqn=abs(qn)
           amn=qn*sosi
           am0=min(abs(amn),1.)*sign(1.,amn)
           if(abs(amn-1.)<epsm) am0=(-(1.-epsm)**2+2.*(1.+epsm)*amn-amn*amn)/(4.*epsm)
           if(abs(amn+1.)<epsm) am0=((1.-epsm)**2+2.*(1.+epsm)*amn+amn*amn)/(4.*epsm)
           aam0=abs(am0)
           if(aam0<epsm) aam0=(epsm*epsm+aam0*aam0)/(2.*epsm)
           if(aqn<epsm*sos) aqn=(epsm*epsm*sos*sos+aqn*aqn)/(2.*epsm*sos)
           omam0=1.-aam0
           
           sosip=sosi
           sosp=sos 
           
          
           
           
           
           
!
!     inviscid flux
!      
      coef=S_j(i,j)
      dpc1=sosip*omam0*dp
      dpc2=am0*dp
      dqnc1=rm*am0*dunormal
      dqnc2=rm*sosp*omam0*dunormal
      
      adr =coef*(              dpc1              +dqnc1+aqn*drou)
      adru=coef*(sav1n*dpc2+um*dpc1+sav1n*dqnc2+um*dqnc1+aqn*dru)
      adrv=coef*(sav2n*dpc2+vm*dpc1+sav2n*dqnc2+vm*dqnc1+aqn*drv)  
      adre=coef*(  qn0*dpc2+hm*dpc1+  qn0*dqnc2+hm*dqnc1+aqn*dre)    
      
      ulnormal=(sav1*ul+sav2*vl+vjb(i,j))
      urnormal=(sav1*ur+sav2*vr+vjb(i,j))
      rulnormal=rl*ulnormal
      rurnormal=rr*urnormal
      
!
!     flux terms
!
      droj=-0.5*(rulnormal+rurnormal-adr)
      dxj=-0.5*(rulnormal*ul+rurnormal*ur&
     &            +sav1*(pl+pr)-adru) 
      dyj=-0.5*(rulnormal*vl+rurnormal*vr&
     &            +sav2*(pl+pr)-adrv) 
      drej=-0.5*(rulnormal*hl+rurnormal*hr&
     &            -vjb(i,j)*(pl+pr)-adre)
           
      
      F_rho(i,j)=F_rho(i,j)-droj
      F_rho_Et(i,j)=F_rho_Et(i,j)-drej
      F_rho_vx(i,j)=F_rho_vx(i,j)-dxj
      F_rho_vy(i,j)=F_rho_vy(i,j)-dyj
      
      F_rho(i,j-1)=F_rho(i,j-1)+droj
      F_rho_Et(i,j-1)=F_rho_Et(i,j-1)+drej
      F_rho_vx(i,j-1)=F_rho_vx(i,j-1)+dxj
      F_rho_vy(i,j-1)=F_rho_vy(i,j-1)+dyj
      
      
      
      end do
      end do     
!  
      return
      end 
     

!=====================================================================================
       subroutine inviscid_fluxes_lax
!     =========================
      
      use main
      real*8 :: Aco2,Aco,Arx,Ary,Vibc,Sii,Armx,Velpro,Acoa,Eox1,Vjbc,Sjj,Army,Eoy1,Ste
      real*8 :: ttg, Alsinv,epsm
      real*8 :: DR,DL,UL,UR,VL,VR,PL,PR,RL,RR,VamL,VamR,HL,HR,AcL,AcR
      real*8 :: RRoRL,RRoRLp1
      real*8 :: Rm,Um,Vm,Hm,Vm2,Am2,Am
      real*8 :: sav1,sav2,sav,sav1N,sav2N
      real*8 :: DUnormal,Drou,Dp,DRU,DRV,DRE
      real*8 :: sos,sos1,QNo,QN
      
       
      real*8 :: DU,DV
      real*8 :: aqn,amn,am0,aam0,omam0,sosIP,sosp,coef,Dpc1,Dpc2,Dqnc1,Dqnc2
      real*8 :: Adr,AdrU,AdrE,ULnormal,Urnormal,RULnormal,RURnormal
      real*8 :: Droi,Dxi,Drei,Droj,Dxj,Dyj,Drej
      real*8 :: EL,ER
      real*8 :: Adrv,Dyi
      real*8 :: coeX,coeY
      
      epsm=0.15
      
     
!
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
!
      do i=ib,im    
!    
!     left&right states
!
      
      
      ul=uil(i,j)
      ur=uir(i,j) 
     
      vl=vil(i,j)
      vr=vir(i,j)
     
     
      pl=max(pil(i,j),EPSILON)
      pr=max(pir(i,j),EPSILON)
     
      rl=max(ril(i,j),EPSILON)
      rr=max(rir(i,j),EPSILON)
           
      vaml=ul*ul+vl*vl
      vamr=ur*ur+vr*vr
      hl=pl*Gamma/(rl*Gamma1)+0.5*vaml
      hr=pr*Gamma/(rr*Gamma1)+0.5*vamr 
      el=pl/(rl*Gamma1)+0.5*vaml
      er=pr/(rr*Gamma1)+0.5*vamr 
      acl=sqrt(Gamma*pl/rl)
      acr=sqrt(Gamma*pr/rr)

                  
! 
!     
!     surface area vectors
!
      sav1=Sn_X_i(i,j)
      sav2=Sn_Y_i(i,j)
      
      sav=S_i(i,j)+EPSILON
      sav1n=sav1/sav
      sav2n=sav2/sav
      

      ulnormal=(sav1*ul+sav2*vl+vib(i,j))
      urnormal=(sav1*ur+sav2*vr+vib(i,j))
      
      

	  coex=max(abs(ulnormal),abs(urnormal))+max(acl,acr)*sav

      rulnormal=rl*ulnormal
      rurnormal=rr*urnormal


      
!
!     flux terms
!
      droi=-0.5*(rulnormal+rurnormal-coex*(rr-rl))
       dxi=-0.5*(rulnormal*ul+rurnormal*ur&
     &            +sav1*(pl+pr)-coex*(rr*ur-rl*ul)) 
       dyi=-0.5*(rulnormal*vl+rurnormal*vr&
     &            +sav2*(pl+pr)-coex*(rr*vr-rl*vl))  
      drei=-0.5*(rulnormal*hl+rurnormal*hr&
     &            -vib(i,j)*(pl+pr)-coex*(rr*er-rl*el))    
    
      
     
      F_rho(i,j)=F_rho(i,j)-droi
      F_rho_Et(i,j)=F_rho_Et(i,j)-drei
      F_rho_vx(i,j)=F_rho_vx(i,j)-dxi
      F_rho_vy(i,j)=F_rho_vy(i,j)-dyi
     
      F_rho(i-1,j)=F_rho(i-1,j)+droi
      F_rho_Et(i-1,j)=F_rho_Et(i-1,j)+drei
      F_rho_vx(i-1,j)=F_rho_vx(i-1,j)+dxi
      F_rho_vy(i-1,j)=F_rho_vy(i-1,j)+dyi
      
      
      
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
      ul=ujl(i,j)
      ur=ujr(i,j) 
     
      vl=vjl(i,j)
      vr=vjr(i,j)
     
     
      pl=max(pjl(i,j),EPSILON)
      pr=max(pjr(i,j),EPSILON)
     
      rl=max(rjl(i,j),EPSILON)
      rr=max(rjr(i,j),EPSILON) 
          
      vaml=ul*ul+vl*vl
      vamr=ur*ur+vr*vr
      hl=pl*Gamma/(rl*Gamma1)+0.5*vaml
      hr=pr*Gamma/(rr*Gamma1)+0.5*vamr 
	    el=pl/(rl*Gamma1)+0.5*vaml
      er=pr/(rr*Gamma1)+0.5*vamr 
      acl=sqrt(Gamma*pl/rl)
      acr=sqrt(Gamma*pr/rr)

                  
! 
!     
!     surface area vectors
!
      sav1=Sn_X_j(i,j)
      sav2=Sn_Y_j(i,j)
      sav=S_j(i,j)+EPSILON
      sav1n=sav1/sav
      sav2n=sav2/sav
          

!
!     inviscid flux
!      
       
      
      ulnormal=(sav1*ul+sav2*vl+vjb(i,j))
      urnormal=(sav1*ur+sav2*vr+vjb(i,j))
	    coey=max(abs(ulnormal),abs(urnormal))+max(acl,acr)*sav
      rulnormal=rl*ulnormal
      rurnormal=rr*urnormal
      
      
!
!     flux terms
!
      droj=-0.5*(rulnormal+rurnormal-coey*(rr-rl))
       dxj=-0.5*(rulnormal*ul+rurnormal*ur&
     &            +sav1*(pl+pr)-coey*(rr*ur-rl*ul)) 
       dyj=-0.5*(rulnormal*vl+rurnormal*vr&
     &            +sav2*(pl+pr)-coey*(rr*vr-rl*vl))  
      drej=-0.5*(rulnormal*hl+rurnormal*hr&
     &            -vjb(i,j)*(pl+pr)-coey*(rr*er-rl*el)) 
           
      
      F_rho(i,j)=F_rho(i,j)-droj
      F_rho_Et(i,j)=F_rho_Et(i,j)-drej
      F_rho_vx(i,j)=F_rho_vx(i,j)-dxj
      F_rho_vy(i,j)=F_rho_vy(i,j)-dyj
      
      F_rho(i,j-1)=F_rho(i,j-1)+droj
      F_rho_Et(i,j-1)=F_rho_Et(i,j-1)+drej
      F_rho_vx(i,j-1)=F_rho_vx(i,j-1)+dxj
      F_rho_vy(i,j-1)=F_rho_vy(i,j-1)+dyj
      
      
      
      end do
      end do     
!  
      return
      end

!==============================================================================
      subroutine RK2
      !=========================
      
      ! compute the right hand side of the Navier-Stokes equations
      
      use main
      
      integer :: nrk
      save nrk
      data nrk /1/
      REAL*8:: rkpa(2)
      real*8 :: temp
      rkpa(1)=1.0
      rkpa(2)=0.5
!      
!     update variables
!
      
      
      do  j=jb,jm-1
      do  i=ib,im-1
      
      temp=step(i,j)/Vcell(i,j)
      rho(i,j)=(1.-rkpa(nrk))*Rho_m1(i,j)+rkpa(nrk)*rho(i,j)+&
     &              rkpa(nrk)*F_rho(i,j)*temp
     
      rho_Et(i,j)=(1.-rkpa(nrk))*Rho_Et_m1(i,j)+rkpa(nrk)*rho_Et(i,j)+&
     &              rkpa(nrk)*F_rho_Et(i,j)*temp
     
      rho_vx(i,j)=(1.-rkpa(nrk))*Rho_vx_m1(i,j)+rkpa(nrk)*rho_vx(i,j)+&
     &              rkpa(nrk)*F_rho_vx(i,j)*temp
     
       rho_vy(i,j)=(1.-rkpa(nrk))*Rho_vy_m1(i,j)+rkpa(nrk)*rho_vy(i,j)+&
     &              rkpa(nrk)*F_rho_vy(i,j)*temp
     
      end do
      end do
      nrk=2-mod(nrk+1,2)
      
      return
      end
      
      !==============================================================================
      subroutine RK4
      !=========================
      
      ! compute the right hand side of the Navier-Stokes equations
      
      use main
      
      integer :: nrk
      save nrk
      data nrk /1/
      REAL*8:: rkpa(4)
      real*8 :: temp
      rkpa(1)=1.0/4.0
      rkpa(2)=1.0/3.0
      rkpa(3)=0.5
      rkpa(4)=1.0
!      
!     update variables
!
      
      
      do  j=jb,jm-1
      do  i=ib,im-1
      
      temp=step(i,j)/Vcell(i,j)
      rho(i,j)      =rho_m1(i,j)       +rkpa(nrk)*F_rho(i,j)*temp
     
      rho_Et(i,j)   =rho_Et_m1(i,j)    +rkpa(nrk)*F_rho_Et(i,j)*temp
     
      rho_vx(i,j)   =rho_vx_m1(i,j)    +rkpa(nrk)*F_rho_vx(i,j)*temp
     
       rho_vy(i,j)  =rho_vy_m1(i,j)    +rkpa(nrk)*F_rho_vy(i,j)*temp
     
      end do
      end do
      nrk=mod(nrk,4)+1
      
      return
      end
      
!====================================================================================
      subroutine update_variables
!     ===========================

      use main                
      real*8 ::TT,Vsrh,RhoInv,VolInv

!
!
!     update variables
!
!
      
      do  j=jb-1,jm
      do  i=ib-1,im
      
!
      RhoInv=1./rho(i,j)
      
      vx(i,j)=rho_vx(i,j)*RhoInv
      vy(i,j)=rho_vy(i,j)*RhoInv
      
      
      vsrh=.5*(vx(i,j)*vx(i,j)+vy(i,j)*vy(i,j))
      p(i,j)=max(0.0,Gamma1*(rho_Et(i,j)-vsrh*rho(i,j)))
      Ht(i,j)=Gamma*(rho_Et(i,j)/rho(i,j)-vsrh)+vsrh
      T(i,j)=p(i,j)*RhoInv/R_air
!
      TT=T(i,j)*T_ref
      Mu_L(i,j)=1.458d-6*abs(TT)**1.5/(TT+110.4)/(Rho_ref*V_ref*L_ref)
      
      Mu_T(i,j)=min(rho(i,j)*KT(i,j)/OmegaT(i,j),5*Mu_L(i,j))
      Mu_E(i,j)=Mu_L(i,j)+Mu_T(i,j) !for turbulence!!
      end do
      end do
     
!
      return
      end
