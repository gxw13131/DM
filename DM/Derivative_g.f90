subroutine Derivative2rd
    ! compute the Shear strain [SSM2]/ Vorticity[Vort]/Cross Derivative[dKdO]  
    use main
    implicit none
    real*8 :: dudx,dudy,dvdx,dvdy
    real*8 :: dKTdx,dKTdy,dOmegaTdx,dOmegaTdy
    do j=jb,jm-1
      do i=ib,im-1
      !d*dx=dx(*W,*E,*N,*S,i,j)
      dudx=dx(uir,uil,ujr,ujl,i,j)
      dudy=dy(uir,uil,ujr,ujl,i,j)
      dvdx=dx(vir,vil,vjr,vjl,i,j)
      dvdy=dy(vir,vil,vjr,vil,i,j)
      dKTdx=dx(KTir,KTil,KTjr,KTjl,i,j)
      dKTdy=dy(KTir,KTil,KTjr,KTjl,i,j)
      dOmegaTdx=dx(OmegaTir,OmegaTil,OmegaTjr,OmegaTjl,i,j)
      dOmegaTdy=dy(OmegaTir,OmegaTil,OmegaTjr,OmegaTjl,i,j)
      SSM2(i,j)=dudx**2+dvdy**2+0.5*(dudy+dvdx)**2
      Vort(i,j)=0.25*(dudy-dvdx)**2
      dKdO(i,j)=dKTdx*dOmegaTdx+dKTdy*dOmegaTdy
      DIV(i,j)=dudx+dvdy
      end do
      end do
    
    
    end subroutine