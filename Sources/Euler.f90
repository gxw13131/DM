!==================================================================!
!======= 2d Euler solver ==========!
!==================================================================!
!                  V 1.0 by REN Yuxin, Oct. 12, 2007
!==================================================================!
!
!
      program Euler
     !======================
     

     use main

!
      open (2,file='Euler.his',status='unknown')
!     bl.his:file that store the computation history
      
!
!     initialization
!     --------------
      write(*,*) '1'
!
      call startup
      
      write(*,*) '2'
!
!     start of main iteration loop
!     ----------------------------
!
      
      is_End=0
!     if is_End==1, the computation stops
      
      istart=0
!     istart:store the advancing time steps for every new start up
!
!     start a new time step
!     ---------------------
 

 100  is_print=0
!     if is_print==1 output the results

      n=n+1
!     n: accumulated time step of the computation
      istart=istart+1
      
      if(mod(n,1).eq.0) then
          write(*,*) 'now, computing the ',n,'th step'
          write(2,*) 'now, computing the ',n,'th step'
      end if
      if(n.eq.nmax) is_End=1
      
      is_print=0
      if(mod(n,iprint).eq.0) is_print=1
      
!     compute the minimum time step to guarantee CFL condition
      call time_step
!     ============= 
!     
      do nunsloop=1,2
!       nunsloop: number of step in Runge-Kutta time stepping scheme 
!
!     solve Euler equations
!     ------------------
!
      call solver(nunsloop)
!     subroutine solver: solve the ns equation for one Runge-Kutta stage      
      
      end do
!
!
      if(.not.steady1) ttime=ttime+step(ib,jb)
!
!     store variables of previous time steps
!     ---------------------------------------
!

      rms1=0.
      rms2=0.
      rms3=0.
      rms4=0.
      
      
      do i=ib,im-1
      do j=jb,jm-1
      
      tinv=1./step(i,j)
      
      r1=abs(rho(i,j)-rho_m1(i,j))*tinv
      r2=abs(rho_vx(i,j)-Rho_vx_m1(i,j))*tinv
      r3=abs(rho_vy(i,j)-Rho_vy_m1(i,j))*tinv
      r4=abs(rho_Et(i,j)-Rho_Et_m1(i,j))*tinv
      
      if(r1>rms1) then
        rsm1=r1
        imax1=i
        jmax1=j
      end if 
      
      if(r2>rms2) then
        rsm2=r2
        imax2=i
        jmax2=j
      end if 
      
      if(r3>rms3) then
        rsm3=r3
        imax3=i
        jmax3=j
      end if
      
      if(r4>rms4) then
        rsm4=r4
        imax4=i
        jmax4=j
      end if  
      
      end do
      end do
      
      if(mod(n,1).eq.0) then
          write(*,*) 'time', ttime*L_ref/V_ref
          write(*,*) rsm1,imax1,jmax1
          write(*,*) rsm2,imax2,jmax2
          write(*,*) rsm3,imax3,jmax3
          write(*,*) rsm4,imax4,jmax4
        
          write(2,*) rsm1,imax1,jmax1
          write(2,*) rsm2,imax2,jmax2
          write(2,*) rsm3,imax3,jmax3
          write(2,*) rsm4,imax4,jmax4
      end if
      
      do i=ib,im
      do j=jb,jm 
!
      Rho_m1(i,j)=rho(i,j)
      Rho_vx_m1(i,j)=rho_vx(i,j)
      Rho_vy_m1(i,j)=rho_vy(i,j)
      Rho_Et_m1(i,j)=rho_Et(i,j)
!
      end do
      end do

      call output
      if(.not.steady1) call outputuns
      
!     ==============
!
      if(is_End.eq.1) then
      close (2)
      write(*,*) 'Normal Termination'
  
      write(*,*) time2-time1
      stop
      endif
      goto 100
      end
!
!     end of main loop
!     ----------------
!
