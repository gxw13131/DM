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
