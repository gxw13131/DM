subroutine inviscid_fluxes_AUSMPWplus
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
     