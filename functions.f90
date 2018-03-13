module functions
        
  Implicit none

Contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Initial conditions during matter dominance 
!The expressions match those found by Martin
!and Domenico without anisotropic stress in
! PRD 80, 083519 (2009)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initial conditions for matter perturbations
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!****Dark matter density perturbation********

function dm(a,H0,k,B,omegam)
         real*8 ::  a, H0, k, B,dm,omegam
         dm = -2.d0*B*( 1.d0 + (a*k**2)/(3.d0*H0**2*omegam))
end function dm

!********************************************

!****Dark matter velocity perturbation*******


function thetammd(a,H0,k,B,omegam)
     real*8 :: a,H0,k,B,thetammd,omegam
     thetammd= (2.d0*B*k**2*sqrt(a))/(3.d0*H0*Sqrt(omegam))
end function thetammd

!********************************************

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initial conditions for dark energy perturbations
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%%%%%%%%%
!Super-sound horizon scales
!%%%%%%%%%%%%%%%%%%%%%%%%%%


!*****Dark energy density perturbation***********

function dd1(a,H0,k,w,B,omegam) 
      real*8 :: a,H0,k,w,B,dd1,omegam
      dd1 = -2.d0*B*(1.d0+w)*a*k**2*(1.d0/(1.d0-3.d0*w) + 3.d0*H0**2*omegam/(k**2*a) )/(3.d0*H0**2*omegam) 
end function dd1

!************************************************

!***Dark energy velocity perturbation************


function thetademdsltsh(a,H0,k,w,B,omegam)
     real*8 :: a,H0,k,B,thetademdsltsh,w,omegam
     thetademdsltsh = 2.d0*B*k**2*Sqrt(a)/(3.d0*H0*sqrt(omegam)) 
end function thetademdsltsh

!************************************************



!%%%%%%%%%%%%%%%%%%%%%%%%
!Sub-sound horizon scales
!%%%%%%%%%%%%%%%%%%%%%%%%

!******Dark energy density perturbation**********

function dd2(w,cs,B) 
      real*8 :: w,cs,B,dd2
      dd2 =-B*(1.d0 + w)/cs**2
end function dd2

!************************************************

!*******Dark energy velocity perturbation********

!The following is the Domenico's expression divided 
! by (1+w).
function thetademdsstsh(a,H0,k,w,cs,B) 
     real*8 :: a,H0,w,cs,B,thetademdsstsh,k,omegam
     thetademdsstsh = 3.d0*H0*B*sqrt(omegam)*(cs**2-w)/(cs**2*sqrt(a))
end function thetademdsstsh


!************************************************

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Some equations we use along with the equations
!above for the initial conditions
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!**** constant delta_0 in Domenico and Martin's paper ***

function d0(B,k,H0,omegam)
         real*8 :: B,k,H0,omegam,d0
         d0 = -(2.d0*B*k**2)/(3.d0*H0**2*omegam) 
end function d0

!*********************************************************

!******* constant for sub-sound-horizon solutions ********

function d0sub(a,B,k,H0,omegam,cs,w,e,f,g)
         real*8 :: a,B,k,H0,omegam,cs,e,f,g,d0sub,w
         ! This is the constant for the full model 
         !d0sub = (3.d0*(1.0d0 + w)*H0**2*omegam*(6.d0*a*ceff(cs,f)**2*k**2 + &
         !4.d0*f*g**2*H0**2*omegam)*d0(B,k,H0,omegam))/(2.d0*a*cs**2*k**2*(4.d0*a*e*k**2 + &
         !3.d0*H0**2*omegam*(3.d0 - 4.d0*e + 3.d0*w)))
         ! The constant for the model with f_pi = 0 is
         !d0sub = 9.d0*(1.d0+w)*H0**2*omegam*ceff(cs,f)**2*d0(B,k,H0,omegam)/(cs**2*(4.d0*e*a*ceff(cs,f)**2*k**2 +&
         !9.d0*H0**2*omegam*(1.d0+w)))
         !Below the constant for the model with e_pi = 0
         d0sub = d0(B,k,H0,omegam)
end function d0sub

!*********************************************************

!*** constant for super-horizon solutions ****************

function d0sup(a,B,k,H0,omegam,cs,w,e,f,g)
         real*8 :: a,B,k,H0,omegam,cs,w,e,f,g,d0sup
         d0sup = -2.d0*B*a*k**4*(1.d0+w)*(-3.d0+4.d0*f)*(1.d0/(1.d0-3.d0*w) &
         + 3.d0*H0**2*omegam/(k**2*a))/(9.d0*H0**4*omegam**2*(4.d0*e-3.d0*(1.d0+w)))
end function d0sup

function d02sup(B,k,H0,omegam,w,e,f)
         real*8 :: B,k,H0,omegam,w,e,f,d02sup
         d02sup = -2.0d0*B*k**2*(1.0d0+w)*(4.0d0*f-3.d0)/(3.d0*H0**2*omegam*(4.d0*e-3.d0*(1.d0+w)))
end function d02sup

!*********************************************************


! **** Effective sound speed *******************

function ceff(cs,f)
         real*8 :: cs,f,ceff
         ceff = sqrt(cs**2 - 2.d0*f/3.d0)
end function ceff

!***********************************************

!****** The conformal Hubble parameter *********


function cHubble(X,H0,omegam,w) 
          Real*8 :: X,H0,omegam,w,cHubble    
     cHubble= H0*Sqrt(omegam/X + (1.d0 - omegam)*X**(-1.d0 - 3.d0*w))
end function cHubble

!***********************************************


!*** Derivative of the conformal Hubble parameter **


function dcHubble(X,H0,omegam,w) 
     Real*8 :: X,H0,omegam,w,dcHubble
     dcHubble=-(H0**2*(omegam + (1.d0 - omegam)*(1.d0 + 3.d0*w)/X**(3.d0*w)))/(2.d0*X**2*cHubble(X,H0,omegam,w))
end function dcHubble


!***************************************************

!*** Scale factor at which the mode crosses the Hubble horizon in a flat MD universe **** 


function H1(H0,k,speedL,omegam)
     Real*8 :: H0,k,H1,speedL,omegam
     H1=H0**2*omegam/(speedL**2*k**2)
end function H1

!*****************************************************************************************

!*** Scale factor at which the mode crosses the effective sound horizon in a flat MD universe ******


function H2(H0,k,cs,f,omegam) 
     Real*8 :: H0,k,cs,H2,f,omegam
     H2=H0**2*omegam/(ceff(cs,f)**2*k**2)
end function H2

!*****************************************************************************************

!********* Gauge-invariannt comoving density contrast for matter *************************

function Dem(X,H0,omegam,w,k,y1,y3) 
    Real*8 :: X,H0,omegam,w,k,y1,y3,Dem
    Dem = y1 + 3.d0*cHubble(X,H0,omegam,w)*y3/k**2
end function Dem

!*****************************************************************************************

!******** Gauge-invariant comoving density constrast for dark energy *********************

function Dd(X,H0,omegam,w,k,y2,y4) 
    Real*8 :: X,H0,omegam,w,k,y2,y4,Dd
    Dd = y2 + 3.d0*cHubble(X,H0,omegam,w)*(1.d0 + w)*y4/k**2
end function Dd

!*****************************************************************************************

!******** The model of dark energy anisotropic stress ************************************

function sigma(X,H0,omegam,w,k,e,f,g) 
    Real*8 :: X,H0,omegam,w,k,y1,y2,y3,y4,sigma,e,f,g
    sigma = (2.d0*(e*Dem(X,H0,omegam,w,k,y1,y2) + (f + g*cHubble(X,H0,omegam,w)**2/k**2)*Dd(X,H0,omegam,w,k,y2,y4)))/(3.d0*(1.d0+w))
end function sigma

!*****************************************************************************************


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Expressions for the potentials and their 
! initial conditions which are set during 
! matter dominance following solutions 
! found in the paper by Guillermo, Martin 
! and Lukas.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!**************


!function Bf1(k,kp,As,ns,Pi) !Initial potential 
 !    Real*8 :: k,kp,As,Pi,ns,Bf1
  !   Bf1=-(3.d0/5.d0)*Pi*Sqrt(2.0D0)/Sqrt(k**3)*Sqrt(As)*Sqrt((k/kp)**(ns-1.d0))
   !   Bf1=-(3.d0/2.d0)*(71.0d0)**2/k**2
!end function Bf1


!**** Initial condition for the gravitational potential in the MD era ******* 


function Bf2(k,kp,ks,As,ns,Pi) 
     Real*8 :: k,kp,ks,As,ns,Bf2,Pi
     Bf2=-(3.d0*Pi*Sqrt(2.0D0)*Sqrt(As)*Sqrt((k/kp)**(ns-1.d0)))/(5.d0*Sqrt(k**3)*(k**2/(2.d0*ks**2*log(k/ks + 1.d0)) + 1.d0))
end function Bf2

!****************************************************************************

!********************* Potential Phi ****************************************

function phi(X,k,H0,omegam,w,y1,y2,y3,y4)
     Real*8 :: X,k,H0,omegam,w,y1,y2,y3,y4,phi
     phi=-((3.d0*H0**2)/(2.d0*k**2))*(omegam*Dem(X,H0,omegam,w,k,y1,y3)/X + &
          (1 - omegam)*Dd(X,H0,omegam,w,k,y2,y4)/X**(1.d0 + 3.d0*w))
end function phi

!****************************************************************************

!********* Potential Psi ****************************************************

function psi(X,k,H0,omegam,w,y1,y2,y3,y4,e,f,g) 
     Real*8 :: X,k,H0,omegam,w,y1,y2,y3,y4,e,f,g,psi
         psi= -9.d0*H0**2*(1.d0+w)*sigma(X,H0,omegam,w,k,e,f,g)*(1.d0-omegam)/(2.d0*k**2*X**(1.d0+3.d0*w)) &
              + phi(X,k,H0,omegam,w,y1,y2,y3,y4)
end function psi

!****************************************************************************

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Solutions for dark energy perturbations in MD when 
!considering  the model of anisotropic stress defined above
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Super-horizon scales 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!******** Dark energy density perturbation ******************

function dde1(a,H0,k,B,omegam,w,e,f,cs,g)
         real*8 :: a,H0,k,B,omegam,dde1,w,e,f,g,cs
         dde1 = d02sup(B,k,H0,omegam,w,e,f)*(((4.d0*e - 3.d0*(1.d0+w))*( & 
         3.d0*H0**2*omegam))/(k**2*(-3.d0 + 4.d0*f) ) )
end function dde1

!************************************************************

!**** Dark energy velocity perturbation *********************

function Vd1(a,H0,k,B,omegam,w,e,f,cs,g)
         real*8 :: a,H0,k,B,omegam,Vd1,w,e,f,g,cs
         Vd1 = a*d02sup(B,k,H0,omegam,w,e,f)*(4.d0*e-3.d0*(1.d0+w))/(4.d0*f - 3.d0)
end function Vd1

!************************************************************

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sub-sound horizon scales 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!**** Dark energy density perturbation ********************* 

function dd4(a,H0,k,B,omegam,w,e,f,g,cs)
         real*8 :: a,H0,k,B,omegam,w,e,f,g,dd4,cs
         !dd4 = (a*d0sub(1.d-5,B,k,H0,omegam,cs,w,e,f,g)*(4.d0*a*e*k**2 + 3.d0*H0**2*omegam*(3.d0 &
         !- 4.d0*e + 3.d0+w)))/(6.d0*a*ceff(cs,f)**2*k**2 + 4.d0*f*g**2*H0**2*omegam)
         dd4 = a*d0sub(1.d-5,B,k,H0,omegam,cs,w,e,f,g)*(2.d0*e/(3.d0*ceff(cs,f)**2) + &
          3.d0*H0**2*omegam*(1.d0+w)/(2.d0*a*ceff(cs,f)**2*k**2))
end function dd4

!***********************************************************

!****** Dark energy velocity perturbation ******************

function Vd4(a,H0,k,B,omegam,w,e,f,g,cs)
         real*8 :: a,H0,k,B,omegam,w,e,f,g,Vd4,cs
     !    Vd4 = (1.d0/(2.d0*(3.d0*a*ceff(cs,f)**2*k**2 + 2.d0*f*g**2*H0**2*omegam)**2))*(a*&
      !   d0sub(1.d-5,B,k,H0,omegam,cs,w,e,f,g)*(12.d0*a**2*ceff(cs,f)**2*e*k**4*(1.d0 + 3.d0*ceff(cs,f)**2&
       !  +2.d0*f - 3.d0*w) + a*H0**2*k**2*(8.d0*e*f*g**2*(2.d0 + 3.d0*ceff(cs,f)**2 + 2.d0*f - 3.d0*w)&
        ! + 9.d0*ceff(cs,f)**2*(3.d0*ceff(cs,f)**2 + 2.d0*f - 3.d0*w)*(3.d0 - 4.d0*e +3.d0*w)*omegam &
        ! + 6.d0*f*g**2*H0**4*(1.d0 + 3.d0*ceff(cs,f)**2 + 2.d0*f -3.d0*w)*(3.d0 - 4.d0*e + 3.d0*w)*omegam**2)))
         Vd4 = d0sub(1.d-5,B,k,H0,omegam,cs,w,e,f,g)*(2.d0*a*e*(1.d0+3.d0*ceff(cs,f)**2 + 2.d0*f -3.d0*w)/(3.d0*ceff(cs,f)**2) &
         + 3.d0*H0**2*omegam*(3.d0*ceff(cs,f)**2 + 2.d0*f - 3.d0*w)*(1.d0+w)/(2.d0*k**2*ceff(cs,f)**2))
end function Vd4

!***********************************************************

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Solutions for dark matter and dark energy perturbations 
! on super-horizon scales in dark energy dominance, without 
!external anisotropic stress and with w=-1. We fix constant 
! delta_0' in the draft matching solution for matter velocity
! perturbation. We fix it at redshift z~ 0.7.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!*** The constant in the draft ******************************

function d01(a,B,k,H0,omegam,f,g)
         real*8 :: a,B,k,H0,omegam,f,g,d01
         d01 = -(2.d0*a**(1.d0-2.d0*f)*(3.d0+7.d0*f+4.d0*f**2)*k**2*d0(B,k,H0,omegam))/(9.d0*f*H0**2*(1.d0-omegam))
         !(4.d0*H0**2*d0(B,k,H0,omegam)*omegam)/( a**(2.d0*f)*k**2*((36.d0*f**2*(1.d0+4.d0*f)-&
         !6.d0*f*(1.d0+2.d0*f)*g**2+g**4)/(a**2*f*(-2.d0-6.d0*f+8.d0*f**2)*g**4) - &
         !(3.d0*(12.d0*f*(1.d0+5.d0*f+4.d0*f**2)+(13.d0+32.d0*f+16.d0*f**2)*g**2)*H0**2*(1.d0&
         !-omegam))/((3.d0+19.d0*f+32.d0*f**2+16.d0*f**3)*g**2*k**2) - (18.d0*a**2*f*H0**4*(1.d0&
         !-omegam)**2)/((3.d0+7.d0*f+4.d0*f**2)*k**4)     )
end function d01

!************************************************************


!********Dark matter density perturbation********************

function dm3(a,H0,k,B,omegam,f,g) 
         real*8 ::  a, H0, k, B,dm3,omegam,f,g
         dm3 = (d01(0.6d0,B,k,H0,omegam,f,g)*a**(2.d0*f)/2.d0)*( -(18.d0*a**2*f*H0**4*(1.d0&
         -omegam)**2)/((2.d0+8.d0*f**2/3.d0 + 14.d0*f/3.d0)*k**4) - 3.d0*(6.d0*(2.d0+8.d0*f**2 + 10.d0*f)*f + &
         (13.d0+16.d0*f**2+32.d0*f)*g**2)*H0**2*(1.d0 -omegam)/((2.d0+192.d0*f**2/9.d0 + 32.d0*f**3/3.d0 + 38.d0*f/3.d0)*g**2*k**2))         
end function dm3

!***********************************************************


!********Dark matter velocity perturbation******************

function Vm3(a,H0,k,B,omegam,f,g)
         real*8 :: a,H0,k,B,Vm3,omegam,f,g
         Vm3 = -d01(0.6d0,B,k,H0,omegam,f,g)*a**(2.d0*f)*3.d0*f*H0**2*(1.d0-omegam)/(k**2*(2.d0+8.d0*f**2/3.d0 + 14.d0*f/3.d0))
end function Vm3

!***********************************************************

!******** Dark energy density perturbation******************

function dd3(a,H0,k,B,omegam,f,g)
         real*8 :: a,H0,k,B,omegam,dd3,f,g
         dd3 = d01(0.6d0,B,k,H0,omegam,f,g)*a**(2.d0*f)*(3.d0+2.d0*f)*H0**2*(1.d0-omegam)/(k**2*(1.d0+4.d0*f/3.d0))
end function dd3

!***********************************************************

!*** Dark energy velocity perturbation**********************

function Vd3(a,H0,k,B,omegam,f,g)
         real*8 :: a,H0,k,B,omegam,Vd3,f,g
         Vd3 = d01(0.6d0,B,k,H0,omegam,f,g)*a**(-2.d0+2.d0*f)
end function Vd3

!***********************************************************



end module functions
