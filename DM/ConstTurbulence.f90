module ConstTurbulence
    
real*8 :: Cmu=0.09
real*8 :: beta_=0.09

real*8 :: sigmaK1=1.0/0.85
real*8 :: sigmaK2=1.0
real*8 :: sigmaO1=1.0/0.5
real*8 :: sigmaO2=1/0.856
real*8 :: beta1=0.075
real*8 :: beta2=0.0828
real*8 :: kappa=0.41
real*8 :: a1=0.31
real*8 :: gmt1=0.5532   !beta1/Cmu-kappa**2/sigmaO1/sqrt(Cmu)
real*8 :: gmt2=0.4404   !beta2/Cmu-kappa**2/sigmaO2/sqrt(Cmu)
  
end module ConstTurbulence