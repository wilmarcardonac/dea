Module background
  
  use fgsl
  use, intrinsic :: iso_c_binding
  use fiducial

  Implicit none

Contains

  !##################################
  ! CONFORMAL HUBBLE PARAMETER STARTS
  !##################################

  function conformal_Hubble_parameter(a)

    Implicit none

    Real*8 :: a,conformal_Hubble_parameter

    
    If (MG_parametrisation .eq. 'GR_DE') then

       conformal_Hubble_parameter = H0*Sqrt(Omega_m/a + (1.d0 - Omega_m)*a**(-1.d0 - 3.d0*w0_fld))

    Else if (MG_parametrisation .eq. 'HS_Basilakos') then

       conformal_Hubble_parameter = a*H0*Sqrt(1.d0 - Omega_m + Omega_m/a**3 + (2.d0*a**2*b_fR*(-1.d0 + &
            Omega_m)**2*(12.d0*a**7*(-1.d0 + Omega_m)**2 + 3.d0*a**4*(-1.d0 + Omega_m)*Omega_m - &
            6.d0*a*Omega_m**2))/(4.d0*a**3*(-1.d0 + Omega_m) - Omega_m)**3 - (a**5*b_fR**2*(-1.d0 + &
            Omega_m)**3*(-1024.d0*a**19*(-1.d0 + Omega_m)**6 - 9216.d0*a**16*(-1.d0 + Omega_m)**5*Omega_m + &
            22848.d0*a**13*(-1.d0 + Omega_m)**4*Omega_m**2 - 25408.d0*a**10*(-1.d0 + Omega_m)**3*Omega_m**3 + &
            7452.d0*a**7*(-1.d0 + Omega_m)**2*Omega_m**4 + 4656.d0*a**4*(-1.d0 + Omega_m)*Omega_m**5 - &
            37.d0*a*Omega_m**6))/(-4.d0*a**3*(-1.d0 + Omega_m) + Omega_m)**8)

    Else if (MG_parametrisation .eq. 'Starobinsky_Basilakos') then

       conformal_Hubble_parameter = -2.d10

    Else if ( (MG_parametrisation .eq. 'Savvas') .or. (MG_parametrisation .eq. 'GR_LAMBDA') ) then

       conformal_Hubble_parameter = H0*Sqrt( Omega_m/a + (1.d0 - Omega_m)*a**2 )

    Else

       conformal_Hubble_parameter = -1.d10

    End if

  end function conformal_Hubble_parameter

  !################################
  ! CONFORMAL HUBBLE PARAMETER ENDS
  !################################

  !####################################################
  ! DERIVATIVE OF THE CONFORMAL HUBBLE PARAMETER STARTS
  !#################################################### 

  function derivative_conformal_Hubble_parameter(a) 

    Implicit none 

    Real*8 :: a,derivative_conformal_Hubble_parameter

    If (MG_parametrisation .eq. 'GR_DE') then

       derivative_conformal_Hubble_parameter = -(H0**2*(Omega_m + (1.d0 - Omega_m)*(1.d0 + &
            3.d0*w0_fld)/a**(3.d0*w0_fld)))/(2.d0*a**2*conformal_Hubble_parameter(a))

    Else if (MG_parametrisation .eq. 'HS_Basilakos') then

       derivative_conformal_Hubble_parameter = (a*H0*((-3.d0*Omega_m)/a**4 + (6.d0*a**3*b_fR*(-1.d0 + &
            Omega_m)**2*(24*a**5*(-1 + Omega_m)**2 + 3*a**2*(-1 + Omega_m)*Omega_m))/(4*a**3*(-1 + &
            Omega_m) - Omega_m)**3 + (18*a**2*b_fR*(-1 + Omega_m)**2*(4*a**6*(-1 + Omega_m)**2 + &
            a**3*(-1 + Omega_m)*Omega_m - 2*Omega_m**2))/(4*a**3*(-1 + Omega_m) - Omega_m)**3 - &
            (216*a**5*b_fR*(-1 + Omega_m)**3*(4*a**6*(-1 + Omega_m)**2 + a**3*(-1 + Omega_m)*Omega_m - &
            2*Omega_m**2))/(4*a**3*(-1 + Omega_m) - Omega_m)**4 - (a**5*b_fR**2*(-1 + &
            Omega_m)**3*(-19456*a**18*(-1 + Omega_m)**6 - 147456*a**15*(-1 + Omega_m)**5*Omega_m + &
            297024*a**12*(-1 + Omega_m)**4*Omega_m**2 - 254080*a**9*(-1 + Omega_m)**3*Omega_m**3 + &
            52164*a**6*(-1 + Omega_m)**2*Omega_m**4 + 18624*a**3*(-1 + Omega_m)*Omega_m**5 - &
            37*Omega_m**6))/(-4*a**3*(-1 + Omega_m) + Omega_m)**8 - (96*a**7*b_fR**2*(-1 + &
            Omega_m)**4*(-1024*a**19*(-1 + Omega_m)**6 - 9216*a**16*(-1 + Omega_m)**5*Omega_m + &
            22848*a**13*(-1 + Omega_m)**4*Omega_m**2 - 25408*a**10*(-1 + Omega_m)**3*Omega_m**3 + &
            7452*a**7*(-1 + Omega_m)**2*Omega_m**4 + 4656*a**4*(-1 + Omega_m)*Omega_m**5 - &
            37*a*Omega_m**6))/(-4*a**3*(-1 + Omega_m) + Omega_m)**9 - (5*a**4*b_fR**2*(-1 + &
            Omega_m)**3*(-1024*a**19*(-1 + Omega_m)**6 - 9216*a**16*(-1 + Omega_m)**5*Omega_m + &
            22848*a**13*(-1 + Omega_m)**4*Omega_m**2 - 25408*a**10*(-1 + Omega_m)**3*Omega_m**3 + &
            7452*a**7*(-1 + Omega_m)**2*Omega_m**4 + 4656*a**4*(-1 + Omega_m)*Omega_m**5 - &
            37*a*Omega_m**6))/(-4*a**3*(-1 + Omega_m) + Omega_m)**8))/(2.*Sqrt(1 - Omega_m + Omega_m/a**3 + &
            (6*a**3*b_fR*(-1 + Omega_m)**2*(4*a**6*(-1 + Omega_m)**2 + a**3*(-1 + Omega_m)*Omega_m - &
            2*Omega_m**2))/(4*a**3*(-1 + Omega_m) - Omega_m)**3 - (a**5*b_fR**2*(-1 + &
            Omega_m)**3*(-1024*a**19*(-1 + Omega_m)**6 - 9216*a**16*(-1 + Omega_m)**5*Omega_m + &
            22848*a**13*(-1 + Omega_m)**4*Omega_m**2 - 25408*a**10*(-1 + Omega_m)**3*Omega_m**3 + &
            7452*a**7*(-1 + Omega_m)**2*Omega_m**4 + 4656*a**4*(-1 + Omega_m)*Omega_m**5 - &
            37*a*Omega_m**6))/(-4*a**3*(-1 + Omega_m) + Omega_m)**8))+ H0*Sqrt(1 - Omega_m + &
            Omega_m/a**3 + (6*a**3*b_fR*(-1 + Omega_m)**2*(4*a**6*(-1 + Omega_m)**2 + a**3*(-1 + &
            Omega_m)*Omega_m - 2*Omega_m**2))/(4*a**3*(-1 + Omega_m) - Omega_m)**3 - (a**5*b_fR**2*(-1 + &
            Omega_m)**3*(-1024*a**19*(-1 + Omega_m)**6 - 9216*a**16*(-1 + Omega_m)**5*Omega_m + &
            22848*a**13*(-1 + Omega_m)**4*Omega_m**2 - 25408*a**10*(-1 + Omega_m)**3*Omega_m**3 + &
            7452*a**7*(-1 + Omega_m)**2*Omega_m**4 + 4656*a**4*(-1 + Omega_m)*Omega_m**5 - &
            37*a*Omega_m**6))/(-4*a**3*(-1 + Omega_m) + Omega_m)**8)

    Else if (MG_parametrisation .eq. 'Starobinsky_Basilakos') then

       derivative_conformal_Hubble_parameter = -2.d10

    Else if ( (MG_parametrisation .eq. 'Savvas') .or. (MG_parametrisation .eq. 'GR_LAMBDA') ) then

       derivative_conformal_Hubble_parameter = -(H0**2*(Omega_m - 2.d0*(1.d0 - Omega_m)*a**(3.d0)&
            ))/(2.d0*a**2*conformal_Hubble_parameter(a))

    Else

       derivative_conformal_Hubble_parameter = -1.d10

    End if

  end function derivative_conformal_Hubble_parameter

  !##################################################
  ! DERIVATIVE OF THE CONFORMAL HUBBLE PARAMETER ENDS
  !################################################## 

  !###########################################################
  ! SECOND DERIVATIVE OF THE CONFORMAL HUBBLE PARAMETER STARTS
  !###########################################################

  function second_derivative_conformal_Hubble_parameter(a) 

    Implicit none 

    Real*8 :: a,second_derivative_conformal_Hubble_parameter

    If (MG_parametrisation .eq. 'GR_DE') then

       second_derivative_conformal_Hubble_parameter = (3.d0*H0*Sqrt(a**(-1.d0 - &
            3.d0*w0_fld)*(1.d0 + (-1.d0 + a**(3.d0*w0_fld))*Omega_m))*(a**(6.d0*w0_fld)*Omega_m**2 + &
            (-1.d0 + Omega_m)**2*(1.d0 + w0_fld)*(1.d0 + 3.d0*w0_fld) - 2.d0*a**(3.d0*w0_fld)*(-1.d0 &
            + Omega_m)*Omega_m*(1.d0 + w0_fld*(2.d0 + 3.d0*w0_fld))))/(4.d0*a**2*(1.d0 + &
            (-1.d0 + a**(3.d0*w0_fld))*Omega_m)**2)

    Else if (MG_parametrisation .eq. 'HS_Basilakos') then

       second_derivative_conformal_Hubble_parameter = (H0*((-3*Omega_m)/a**4 + (6*a**3*b_fR*(-1 + &
            Omega_m)**2*(24*a**5*(-1 + Omega_m)**2 + 3*a**2*(-1 + Omega_m)*Omega_m))/(4*a**3*(-1 + &
            Omega_m) - Omega_m)**3 + (18*a**2*b_fR*(-1 + Omega_m)**2*(4*a**6*(-1 + Omega_m)**2 + &
            a**3*(-1 + Omega_m)*Omega_m - 2*Omega_m**2))/(4*a**3*(-1 + Omega_m) - Omega_m)**3 - &
            (216*a**5*b_fR*(-1 + Omega_m)**3*(4*a**6*(-1 + Omega_m)**2 + a**3*(-1 + Omega_m)*Omega_m - &
            2*Omega_m**2))/(4*a**3*(-1 + Omega_m) - Omega_m)**4 - (a**5*b_fR**2*(-1 + &
            Omega_m)**3*(-19456*a**18*(-1 + Omega_m)**6 - 147456*a**15*(-1 + Omega_m)**5*Omega_m + &
            297024*a**12*(-1 + Omega_m)**4*Omega_m**2 - 254080*a**9*(-1 + Omega_m)**3*Omega_m**3 + &
            52164*a**6*(-1 + Omega_m)**2*Omega_m**4 + 18624*a**3*(-1 + Omega_m)*Omega_m**5 - &
            37*Omega_m**6))/(-4*a**3*(-1 + Omega_m) + Omega_m)**8 - (96*a**7*b_fR**2*(-1 + &
            Omega_m)**4*(-1024*a**19*(-1 + Omega_m)**6 - 9216*a**16*(-1 + Omega_m)**5*Omega_m + &
            22848*a**13*(-1 + Omega_m)**4*Omega_m**2 - 25408*a**10*(-1 + Omega_m)**3*Omega_m**3 + &
            7452*a**7*(-1 + Omega_m)**2*Omega_m**4 + 4656*a**4*(-1 + Omega_m)*Omega_m**5 - &
            37*a*Omega_m**6))/(-4*a**3*(-1 + Omega_m) + Omega_m)**9 - (5*a**4*b_fR**2*(-1 + &
            Omega_m)**3*(-1024*a**19*(-1 + Omega_m)**6 - 9216*a**16*(-1 + Omega_m)**5*Omega_m + &
            22848*a**13*(-1 + Omega_m)**4*Omega_m**2 - 25408*a**10*(-1 + Omega_m)**3*Omega_m**3 + &
            7452*a**7*(-1 + Omega_m)**2*Omega_m**4 + 4656*a**4*(-1 + Omega_m)*Omega_m**5 - &
            37*a*Omega_m**6))/(-4*a**3*(-1 + Omega_m) + Omega_m)**8))/&
            Sqrt(1 - Omega_m + Omega_m/a**3 + (6*a**3*b_fR*(-1 + Omega_m)**2*(4*a**6*(-1 + &
            Omega_m)**2 + a**3*(-1 + Omega_m)*Omega_m - 2*Omega_m**2))/(4*a**3*(-1 + Omega_m) - Omega_m)**3 - &
            (a**5*b_fR**2*(-1 + Omega_m)**3*(-1024*a**19*(-1 + Omega_m)**6 - 9216*a**16*(-1 + &
            Omega_m)**5*Omega_m + 22848*a**13*(-1 + Omega_m)**4*Omega_m**2 - 25408*a**10*(-1 + &
            Omega_m)**3*Omega_m**3 + 7452*a**7*(-1 + Omega_m)**2*Omega_m**4 + 4656*a**4*(-1 + &
            Omega_m)*Omega_m**5 - 37*a*Omega_m**6))/(-4*a**3*(-1 + Omega_m) + Omega_m)**8) + &
            a*H0*(-((-3*Omega_m)/a**4 + (6*a**3*b_fR*(-1 + Omega_m)**2*(24*a**5*(-1 + Omega_m)**2 + &
            3*a**2*(-1 + Omega_m)*Omega_m))/(4*a**3*(-1 + Omega_m) - Omega_m)**3 + &
            (18*a**2*b_fR*(-1 + Omega_m)**2*(4*a**6*(-1 + Omega_m)**2 + a**3*(-1 + Omega_m)*Omega_m - &
            2*Omega_m**2))/(4*a**3*(-1 + Omega_m) - Omega_m)**3 - (216*a**5*b_fR*(-1 + Omega_m)**3*(4*a**6*(-1 &
            + Omega_m)**2 + a**3*(-1 + Omega_m)*Omega_m - 2*Omega_m**2))/(4*a**3*(-1 + Omega_m) - Omega_m)**4 - &
            (a**5*b_fR**2*(-1 + Omega_m)**3*(-19456*a**18*(-1 + Omega_m)**6 - 147456*a**15*(-1 + &
            Omega_m)**5*Omega_m + 297024*a**12*(-1 + Omega_m)**4*Omega_m**2 - 254080*a**9*(-1 + &
            Omega_m)**3*Omega_m**3 + 52164*a**6*(-1 + Omega_m)**2*Omega_m**4 + 18624*a**3*(-1 + &
            Omega_m)*Omega_m**5 - 37*Omega_m**6))/(-4*a**3*(-1 + Omega_m) + Omega_m)**8&
            - (96*a**7*b_fR**2*(-1 + Omega_m)**4*(-1024*a**19*(-1 + Omega_m)**6 - 9216*a**16*(-1 + &
            Omega_m)**5*Omega_m + 22848*a**13*(-1 + Omega_m)**4*Omega_m**2 - 25408*a**10*(-1 + &
            Omega_m)**3*Omega_m**3 + 7452*a**7*(-1 + Omega_m)**2*Omega_m**4 + 4656*a**4*(-1 + &
            Omega_m)*Omega_m**5 - 37*a*Omega_m**6))/(-4*a**3*(-1 + Omega_m) + Omega_m)**9&
            - (5*a**4*b_fR**2*(-1 + Omega_m)**3*(-1024*a**19*(-1 + Omega_m)**6 - 9216*a**16*(-1 + &
            Omega_m)**5*Omega_m + 22848*a**13*(-1 + Omega_m)**4*Omega_m**2 - 25408*a**10*(-1 + &
            Omega_m)**3*Omega_m**3 + 7452*a**7*(-1 + Omega_m)**2*Omega_m**4 + 4656*a**4*(-1 + &
            Omega_m)*Omega_m**5 - 37*a*Omega_m**6))/(-4*a**3*(-1 + Omega_m) + Omega_m)**8)**2/(4.*(1 - &
            Omega_m + Omega_m/a**3 + (6*a**3*b_fR*(-1 + Omega_m)**2*(4*a**6*(-1 + Omega_m)**2 + a**3*(-1 + &
            Omega_m)*Omega_m - 2*Omega_m**2))/(4*a**3*(-1 + Omega_m) - Omega_m)**3 - (a**5*b_fR**2*(-1 + &
            Omega_m)**3*(-1024*a**19*(-1 + Omega_m)**6 - 9216*a**16*(-1 + Omega_m)**5*Omega_m + 22848*a**13*(-1 + &
            Omega_m)**4*Omega_m**2 - 25408*a**10*(-1 + Omega_m)**3*Omega_m**3 + 7452*a**7*(-1 + &
            Omega_m)**2*Omega_m**4 + 4656*a**4*(-1 + Omega_m)*Omega_m**5 - 37*a*Omega_m**6))/(-4*a**3*(-1 + &
            Omega_m) + Omega_m)**8)**1.5) + ((12*Omega_m)/a**5 + (6*a**3*b_fR*(-1 + Omega_m)**2*(120*a**4*(-1 + &
            Omega_m)**2 + 6*a*(-1 + Omega_m)*Omega_m))/(4*a**3*(-1 + Omega_m) - Omega_m)**3 + &
            (36*a**2*b_fR*(-1 + Omega_m)**2*(24*a**5*(-1 + Omega_m)**2 + 3*a**2*(-1 + &
            Omega_m)*Omega_m))/(4*a**3*(-1 + Omega_m) - Omega_m)**3 - (432*a**5*b_fR*(-1 + &
            Omega_m)**3*(24*a**5*(-1 + Omega_m)**2 + 3*a**2*(-1 + Omega_m)*Omega_m))/(4*a**3*(-1 + &
            Omega_m) - Omega_m)**4 + (36*a*b_fR*(-1 + Omega_m)**2*(4*a**6*(-1 + Omega_m)**2 + a**3*(-1 + &
            Omega_m)*Omega_m - 2*Omega_m**2))/(4*a**3*(-1 + Omega_m) - Omega_m)**3 - &
            (1728*a**4*b_fR*(-1 + Omega_m)**3*(4*a**6*(-1 + Omega_m)**2 + a**3*(-1 + Omega_m)*Omega_m - &
            2*Omega_m**2))/(4*a**3*(-1 + Omega_m) - Omega_m)**4 + (10368*a**7*b_fR*(-1 + Omega_m)**4*(4*a**6*(-1 + &
            Omega_m)**2 + a**3*(-1 + Omega_m)*Omega_m - 2*Omega_m**2))/(4*a**3*(-1 + Omega_m) - Omega_m)**5 - &
            (a**5*b_fR**2*(-1 + Omega_m)**3*(-350208*a**17*(-1 + Omega_m)**6 - 2211840*a**14*(-1 + &
            Omega_m)**5*Omega_m + 3564288*a**11*(-1 + Omega_m)**4*Omega_m**2 - 2286720*a**8*(-1 + &
            Omega_m)**3*Omega_m**3 + 312984*a**5*(-1 + Omega_m)**2*Omega_m**4 + 55872*a**2*(-1 + &
            Omega_m)*Omega_m**5))/(-4*a**3*(-1 + Omega_m) + Omega_m)**8 - (192*a**7*b_fR**2*(-1 + &
            Omega_m)**4*(-19456*a**18*(-1 + Omega_m)**6 - 147456*a**15*(-1 + Omega_m)**5*Omega_m + &
            297024*a**12*(-1 + Omega_m)**4*Omega_m**2 - 254080*a**9*(-1 + Omega_m)**3*Omega_m**3 + &
            52164*a**6*(-1 + Omega_m)**2*Omega_m**4 + 18624*a**3*(-1 + Omega_m)*Omega_m**5 - &
            37*Omega_m**6))/(-4*a**3*(-1 + Omega_m) + Omega_m)**9 - (10*a**4*b_fR**2*(-1 + &
            Omega_m)**3*(-19456*a**18*(-1 + Omega_m)**6 - 147456*a**15*(-1 + Omega_m)**5*Omega_m + &
            297024*a**12*(-1 + Omega_m)**4*Omega_m**2 - 254080*a**9*(-1 + Omega_m)**3*Omega_m**3 + &
            52164*a**6*(-1 + Omega_m)**2*Omega_m**4 + 18624*a**3*(-1 + Omega_m)*Omega_m**5 - &
            37*Omega_m**6))/(-4*a**3*(-1 + Omega_m) + Omega_m)**8 - (10368*a**9*b_fR**2*(-1 + &
            Omega_m)**5*(-1024*a**19*(-1 + Omega_m)**6 - 9216*a**16*(-1 + Omega_m)**5*Omega_m + &
            22848*a**13*(-1 + Omega_m)**4*Omega_m**2 - 25408*a**10*(-1 + Omega_m)**3*Omega_m**3 + &
            7452*a**7*(-1 + Omega_m)**2*Omega_m**4 + 4656*a**4*(-1 + Omega_m)*Omega_m**5 - &
            37*a*Omega_m**6))/(-4*a**3*(-1 + Omega_m) + Omega_m)**10 - (1152*a**6*b_fR**2*(-1 + &
            Omega_m)**4*(-1024*a**19*(-1 + Omega_m)**6 - 9216*a**16*(-1 + Omega_m)**5*Omega_m + &
            22848*a**13*(-1 + Omega_m)**4*Omega_m**2 - 25408*a**10*(-1 + Omega_m)**3*Omega_m**3 + &
            7452*a**7*(-1 + Omega_m)**2*Omega_m**4 + 4656*a**4*(-1 + Omega_m)*Omega_m**5 - &
            37*a*Omega_m**6))/(-4*a**3*(-1 + Omega_m) + Omega_m)**9 - (20*a**3*b_fR**2*(-1 + &
            Omega_m)**3*(-1024*a**19*(-1 + Omega_m)**6 - 9216*a**16*(-1 + Omega_m)**5*Omega_m + &
            22848*a**13*(-1 + Omega_m)**4*Omega_m**2 - 25408*a**10*(-1 + Omega_m)**3*Omega_m**3 + &
            7452*a**7*(-1 + Omega_m)**2*Omega_m**4 + 4656*a**4*(-1 + Omega_m)*Omega_m**5 - &
            37*a*Omega_m**6))/(-4*a**3*(-1 + Omega_m) + Omega_m)**8)/(2.*Sqrt(1 - Omega_m + Omega_m/a**3 + &
            (6*a**3*b_fR*(-1 + Omega_m)**2*(4*a**6*(-1 + Omega_m)**2 + a**3*(-1 + Omega_m)*Omega_m - &
            2*Omega_m**2))/(4*a**3*(-1 + Omega_m) - Omega_m)**3 - (a**5*b_fR**2*(-1 + &
            Omega_m)**3*(-1024*a**19*(-1 + Omega_m)**6 - 9216*a**16*(-1 + Omega_m)**5*Omega_m + &
            22848*a**13*(-1 + Omega_m)**4*Omega_m**2 - 25408*a**10*(-1 + Omega_m)**3*Omega_m**3 + &
            7452*a**7*(-1 + Omega_m)**2*Omega_m**4 + 4656*a**4*(-1 + Omega_m)*Omega_m**5 - &
            37*a*Omega_m**6))/(-4*a**3*(-1 + Omega_m) + Omega_m)**8)))

    Else if (MG_parametrisation .eq. 'Starobinsky_Basilakos') then

       second_derivative_conformal_Hubble_parameter = -2.d10

    Else if ( (MG_parametrisation .eq. 'Savvas') .or. (MG_parametrisation .eq. 'GR_LAMBDA') ) then

       second_derivative_conformal_Hubble_parameter = (3.d0*H0*Sqrt(a**2*(1.d0 + (-1.d0 + &
            a**(-3.d0) )*Omega_m))*(a**(-6.d0)*Omega_m**2 - 4.d0/a**3.d0*(-1.d0 + Omega_m)*Omega_m ))/(4.d0*a**2*(1.d0 + &
            (-1.d0 + a**(-3.d0))*Omega_m)**2)

    Else

       second_derivative_conformal_Hubble_parameter = -1.d10

    End if

  end function second_derivative_conformal_Hubble_parameter

  !#########################################################
  ! SECOND DERIVATIVE OF THE CONFORMAL HUBBLE PARAMETER ENDS
  !#########################################################


  !###########################################################
  ! THIRD DERIVATIVE OF THE CONFORMAL HUBBLE PARAMETER STARTS
  !###########################################################

  function third_derivative_conformal_Hubble_parameter(a) 

    Implicit none 

    Real*8 :: a,third_derivative_conformal_Hubble_parameter

    If (MG_parametrisation .eq. 'GR_DE') then

       third_derivative_conformal_Hubble_parameter = 0.d0 

    Else if (MG_parametrisation .eq. 'HS_Basilakos') then

       third_derivative_conformal_Hubble_parameter = 0.d0

    Else if (MG_parametrisation .eq. 'Starobinsky_Basilakos') then

       third_derivative_conformal_Hubble_parameter = 0.d0

    Else if ( (MG_parametrisation .eq. 'Savvas') .or. (MG_parametrisation .eq. 'GR_LAMBDA') ) then

       third_derivative_conformal_Hubble_parameter = (3*H0*Omega_m*(-32*a**6*(-1 + Omega_m)**2 + &
            10*a**3*(-1 + Omega_m)*Omega_m - 5*Omega_m**2))/(8.*a**6*((a**3 + Omega_m - &
            a**3*Omega_m)/a)**2.5)

    Else

       third_derivative_conformal_Hubble_parameter = 0.d0

    End if

  end function third_derivative_conformal_Hubble_parameter

  !#########################################################
  ! THIRD DERIVATIVE OF THE CONFORMAL HUBBLE PARAMETER ENDS
  !#########################################################

  !##########################
  ! Eq. (3.10) BATTYE'S PAPER
  !##########################

  function g_K(a)

    Implicit none

    Real*8 :: a,g_K,K

    K = wavenumber_k/conformal_Hubble_parameter(a)

    g_K = 1.d0 + K**2/3.d0/E_H(a)

  End function g_K
  
  !##########################
  ! Eq. (2.4) BATTYE'S PAPER
  !##########################

  function E_H(a)

    Implicit none

    Real*8 :: a,E_H

    E_H = 1.d0 -a*derivative_conformal_Hubble_parameter(a)/conformal_Hubble_parameter(a)

  End function E_H

  !###########################
  ! Eq. (4.12a) BATTYE'S PAPER
  !###########################

  function Zeta_de(a)

    Implicit none

    Real*8 :: a,Zeta_de

    Zeta_de = (4.d0*g_K(a) - 1.d0)/3.d0/g_K(a) 

  End function Zeta_de

  !###########################
  ! Eq. (4.12b) BATTYE'S PAPER
  !###########################

  function Zeta_m(a)

    Implicit none

    Real*8 :: a,Zeta_m

    Zeta_m = (g_K(a) - 1.d0)/3.d0/g_K(a) 

  End function Zeta_m

  !####################
  ! RICCI SCALAR STARTS
  !####################

  function ricci_scalar(a)

    Implicit none

    Real*8 :: a,ricci_scalar,H,Hprime

    H = conformal_Hubble_parameter(a) 

    Hprime = derivative_conformal_Hubble_parameter(a)

    ricci_scalar = 6.d0*(H**2 + a*H*Hprime )/a**2

  end function ricci_scalar

  function ricci_scalar_fgsl(a, params) bind(c)
    
    real(c_double), value :: a
    type(c_ptr), value :: params
    real(c_double) :: ricci_scalar_fgsl

    ricci_scalar_fgsl = ricci_scalar(a)

  end function ricci_scalar_fgsl

  !##################
  ! RICCI SCALAR ENDS
  !##################

  !######################################
  ! DERIVATIVE OF THE RICCI SCALAR STARTS
  !######################################

  function ricci_scalar_prime(a)

    real(fgsl_double) :: result, abserr
    integer(fgsl_int) :: status
    type(fgsl_function) :: pwr
    real(fgsl_double) :: ricci_scalar_prime,a

    pwr = fgsl_function_init(ricci_scalar_fgsl, c_null_ptr)

    status = fgsl_deriv_central (pwr, a , 1.E-6_fgsl_double, &
         result, abserr)

    ricci_scalar_prime = result

  end function ricci_scalar_prime

  !####################################
  ! DERIVATIVE OF THE RICCI SCALAR ENDS
  !####################################

  !###############
  ! f(R(a)) STARTS
  !###############

  function fMG(a)

    Implicit none

    Real*8 :: fMG,a,R,x,y
    Real*8,parameter :: b2 = (-7.d0+sqrt(73.d0))/12.d0
    Real*8,parameter :: a2 = 1.d0 + b2
    Real*8,parameter :: b = 1.5d0 + b2
    Real*8,parameter :: c = 13.d0/6.d0 + 2.d0*b2

    R = ricci_scalar(a)

    x = Lambda/(R-3.d0*Lambda)

    y = R/(R-3.d0*Lambda)

    If (MG_parametrisation .eq. 'GR_DE') then

       fMG = 0.d0

    Else if (MG_parametrisation .eq. 'GR_LAMBDA') then

       fMG = -2.d0*Lambda

    Else if (MG_parametrisation .eq. 'HS_Basilakos') then

       fMG = -2.d0*Lambda/(1.d0 + b_fR*Lambda/R)

    Else if (MG_parametrisation .eq. 'Starobinsky_Basilakos') then

       fMG = 1.d10

    Else if (MG_parametrisation .eq. 'Savvas') then

       fMG = -2.d0*Lambda + alpha*H0**2*x**b2*fgsl_sf_hyperg_2f1(b2,b,c,x)

    Else

       fMG = 2.d10

    End if

  end function fMG

  !#############
  ! f(R(a)) ENDS
  !#############

  !#################
  ! f_R(R(a)) STARTS
  !#################

  function fMG_R(a)

    Implicit none

    Real*8 :: fMG_R,a,R,x,y
    Real*8,parameter :: b2 = (-7.d0+sqrt(73.d0))/12.d0
    Real*8,parameter :: a2 = 1.d0 + b2
    Real*8,parameter :: b = 1.5d0 + b2
    Real*8,parameter :: c = 13.d0/6.d0 + 2.d0*b2

    R = ricci_scalar(a)

    x = Lambda/(R-3.d0*Lambda)

    y = R/(R-3.d0*Lambda)

    If (MG_parametrisation .eq. 'GR_DE') then

       fMG_R = 0.d0 

    Else if (MG_parametrisation .eq. 'GR_LAMBDA') then

       fMG_R = 0.d0

    Else if (MG_parametrisation .eq. 'HS_Basilakos') then

       fMG_R = -2.d0*b_fR*Lambda**2/R**2/(1.d0 + b_fR*Lambda/R)**2

    Else if (MG_parametrisation .eq. 'Starobinsky_Basilakos') then

       fMG_R = 0.d0

    Else if (MG_parametrisation .eq. 'Savvas') then

       fMG_R = -b2*H0**2*alpha/R*y**a2*(Lambda/R)**b2*fgsl_sf_hyperg_2f1(a2,b,c,x)

    Else

       fMG_R = 2.d10

    End if

  end function fMG_R

  !###############
  ! f_R(R(a)) ENDS
  !###############

  !#######################
  ! f_R_prime(R(a)) STARTS
  !#######################

  function fMG_R_prime(a)

    Implicit none

    Real*8 :: fMG_R_prime,a,R,x,y,R_prime
    Real*8,parameter :: b2 = (-7.d0+sqrt(73.d0))/12.d0
    Real*8,parameter :: a2 = 1.d0 + b2
    Real*8,parameter :: a3 = 2.d0 + b2
    Real*8,parameter :: b = 1.5d0 + b2
    Real*8,parameter :: c = 13.d0/6.d0 + 2.d0*b2

    R = ricci_scalar(a)

    R_prime = ricci_scalar_prime(a)

    x = Lambda/(R-3.d0*Lambda)

    y = R/(R-3.d0*Lambda)

    If (MG_parametrisation .eq. 'GR_DE') then

       fMG_R_prime = 0.d0

    Else if (MG_parametrisation .eq. 'GR_LAMBDA') then

       fMG_R_prime = 0.d0

    Else if (MG_parametrisation .eq. 'HS_Basilakos') then

       fMG_R_prime = 4.d0*b_fR*Lambda**2*R_prime/(b_fR*Lambda + R)**3

    Else if (MG_parametrisation .eq. 'Starobinsky_Basilakos') then

       fMG_R_prime = 0.d0

    Else if (MG_parametrisation .eq. 'Savvas') then

       fMG_R_prime = b2*a2*H0**2*alpha*(Lambda/R)**b2*(1.d0 + 3.d0*x)**b2*R_prime/(R - 3.d0*Lambda )**2*&
         fgsl_sf_hyperg_2f1(b,a3,c,x)

    Else

       fMG_R_prime = 2.d10

    End if

  end function fMG_R_prime

  function fMG_R_prime_fgsl(a, params) bind(c)
    
    real(c_double), value :: a
    type(c_ptr), value :: params
    real(c_double) :: fMG_R_prime_fgsl

    fMG_R_prime_fgsl = fMG_R_prime(a)

  end function fMG_R_prime_fgsl

  !#####################
  ! f_R_prime(R(a)) ENDS
  !#####################

  !##############################
  ! f_R_double_prime(R(a)) STARTS
  !##############################

  function fMG_R_double_prime(a)

    real(fgsl_double) :: result, abserr
    integer(fgsl_int) :: status
    type(fgsl_function) :: pwr
    real(fgsl_double) :: fMG_R_double_prime,a

    pwr = fgsl_function_init(fMG_R_prime_fgsl, c_null_ptr)

    status = fgsl_deriv_central (pwr, a , 1.E-8_fgsl_double, &
         result, abserr)

    fMG_R_double_prime = result

  end function fMG_R_double_prime

  !############################
  ! f_R_double_prime(R(a)) ENDS
  !############################

  !##################
  ! f_RR(R(a)) STARTS
  !##################

  function fMG_RR(a)

    Implicit none

    Real*8 :: fMG_RR,a,R,x,y
    Real*8,parameter :: b2 = (-7.d0+sqrt(73.d0))/12.d0
    Real*8,parameter :: a2 = 1.d0 + b2
    Real*8,parameter :: b = 1.5d0 + b2
    Real*8,parameter :: c = 13.d0/6.d0 + 2.d0*b2

    R = ricci_scalar(a)

    x = Lambda/(R-3.d0*Lambda)

    y = R/(R-3.d0*Lambda)

    If (MG_parametrisation .eq. 'GR_DE') then

       fMG_RR = 0.d0 

    Else if (MG_parametrisation .eq. 'GR_LAMBDA') then

       fMG_RR = 0.d0

    Else if (MG_parametrisation .eq. 'HS_Basilakos') then

       fMG_RR = 4.d0*b_fR*Lambda**2/(R + b_fR*Lambda)**3

    Else if (MG_parametrisation .eq. 'Starobinsky_Basilakos') then

       fMG_RR = 0.d0

    Else if (MG_parametrisation .eq. 'Savvas') then

       fMG_RR = b2*H0**2*alpha*y**b2*(Lambda/R)**b2/6.d0/(R-4.d0*Lambda)/(R-3.d0*Lambda)**2*(&
            (7.d0+6.d0*b2)*(R-3.d0*Lambda)*fgsl_sf_hyperg_2f1(b2,b,c,x) - &
            (R-6.d0*Lambda)*fgsl_sf_hyperg_2f1(a2,b,c,x))

    Else

       fMG_RR = 2.d10

    End if

  end function fMG_RR

  !################
  ! f_RR(R(a)) ENDS
  !################

  !########################
  ! f_RR_prime(R(a)) STARTS
  !########################

  function fMG_RR_prime(a)

    Implicit none

    Real*8 :: fMG_RR_prime,a,R,x,y,R_prime
    Real*8,parameter :: b2 = (-7.d0+sqrt(73.d0))/12.d0
    Real*8,parameter :: a2 = 1.d0 + b2
    Real*8,parameter :: b = 1.5d0 + b2
    Real*8,parameter :: c = 13.d0/6.d0 + 2.d0*b2

    R = ricci_scalar(a)

    R_prime = ricci_scalar_prime(a)

    x = Lambda/(R-3.d0*Lambda)

    y = R/(R-3.d0*Lambda)

    If (MG_parametrisation .eq. 'GR_DE') then

       fMG_RR_prime = 0.d0 

    Else if (MG_parametrisation .eq. 'GR_LAMBDA') then

       fMG_RR_prime = 0.d0

    Else if (MG_parametrisation .eq. 'HS_Basilakos') then

       fMG_RR_prime = -12.d0*b_fR*Lambda**2*R_prime/(R + b_fR*Lambda)**4

    Else if (MG_parametrisation .eq. 'Starobinsky_Basilakos') then

       fMG_RR_prime = 0.d0

    Else if (MG_parametrisation .eq. 'Savvas') then

       fMG_RR_prime = b2*H0**2*alpha*(Lambda/R)**b2*(1.d0 + 3.d0*x)**b2/36.d0*( &
            ( 7.d0 + 6.d0*b2 )*fgsl_sf_hyperg_2f1(b2,b,c,x)*(36.d0*Lambda - 11.d0*R)/(12.d0*Lambda**2 - &
            7.d0*Lambda*R + R**2)**2 + fgsl_sf_hyperg_2f1(a2,b,c,x)/(3.d0*Lambda - R)**3/(R-4.d0*Lambda)**2*&
            (72.d0*(-2.d0 + b2*(7.d0 + 6.d0*b2))*Lambda**2 - &
            6.d0*(-10.d0 + 7.d0*b2*(7.d0 + 6.d0*b2))*Lambda*R + &
            (-5.d0 + 6.d0*b2*(7.d0 + 6.d0*b2))*R**2))*R_prime

    Else

       fMG_RR_prime = 2.d10

    End if

  end function fMG_RR_prime

  !######################
  ! f_RR_prime(R(a)) ENDS
  !######################

  !###############################
  ! f_RR_double_prime(R(a)) STARTS
  !###############################

  function fMG_RR_double_prime(a)

    Implicit none

    Real*8 :: fMG_RR_double_prime,a

    fMG_RR_double_prime = FR_double_prime(a)

  end function fMG_RR_double_prime

  !#############################
  ! f_RR_double_prime(R(a)) ENDS
  !#############################

  !###########################
  ! DARK ENERGY DENSITY STARTS 
  !###########################

  function dark_energy_density(a) ! ACTUALLY THIS IS THE DARK ENERGY DENSITY MULTIPLIED BY \kappa

    Implicit none

    Real*8 :: a,dark_energy_density,R,H
    
    R = ricci_scalar(a)

    H = conformal_Hubble_parameter(a) 

    If ( (MG_parametrisation .eq. 'GR_DE') .or. (MG_parametrisation .eq. 'Savvas') ) then

       dark_energy_density = 3.d0*H0**2*(1.d0 - Omega_m)*a**(-3.d0 - 3.d0*equation_of_state(a))

    Else if (MG_parametrisation .eq. 'GR_LAMBDA') then

       dark_energy_density = -fMG(a)/2.d0

    Else

       dark_energy_density = (fMG_R(a)*R - fMG(a))/2.d0 - 3.d0*H**2*fMG_R(a)/a**2 - &
            3.d0*H**2*fMG_R_prime(a)/a

    End if

  end function dark_energy_density

  !#########################
  ! DARK ENERGY DENSITY ENDS 
  !#########################

  !#########################
  ! EQUATION OF STATE STARTS
  !#########################

  function equation_of_state(a)

    Implicit none

    Real*8 :: a,equation_of_state,H,H_prime,R

    H = conformal_Hubble_parameter(a)

    H_prime = derivative_conformal_Hubble_parameter(a)

    R = ricci_scalar(a)

    If (MG_parametrisation .eq. 'GR_DE') then

       equation_of_state = w0_fld

    Else if ((MG_parametrisation .eq. 'Savvas') .or. (MG_parametrisation .eq. 'GR_LAMBDA')) then

       equation_of_state = -1.d0 

    Else

       equation_of_state = -1.d0/3.d0 -2.d0/3.d0*(H**2*fMG_R(a) - a**2*fMG(a)/6.d0 - H**2*a*fMG_R_prime(a) - &
            fMG_R_double_prime(a)*H**2*a**2/2.d0 - fMG_R_prime(a)*H_prime*H*a**2/2.d0)/(-H**2*fMG_R(a) - &
            a**2*fMG(a)/6.d0 - H**2*a*fMG_R_prime(a) + a**2*fMG_R(a)*R/6.d0 ) 

    End if

  end function equation_of_state

  !#######################
  ! EQUATION OF STATE ENDS
  !#######################

  !######################################
  ! DERIVATIVE DARK ENERGY DENSITY STARTS
  !######################################

  function derivative_dark_energy_density(a)

    Implicit none

    Real*8 :: a,derivative_dark_energy_density

    If (MG_parametrisation .eq. 'HS_Basilakos') then

       derivative_dark_energy_density = (-18*a**2*b_fR*H0**2*(-1 + Omega_m)**2*Omega_m*(2048*a**21*(24 + 11*b_fR)*(-1 &
            + Omega_m)**7 - 768*a**18*(136 + 77*b_fR)*(-1 + Omega_m)**6*Omega_m + 192 *a**15*(408 + 437*b_fR)*(-1 + &
            Omega_m)**5*Omega_m**2 + 64*a**12*(-420 + 61*b_fR)*(-1 + Omega_m)**4*Omega_m**3 - 24*a**9*(-160 + &
            2561*b_fR)*(-1 + Omega_m)**3*Omega_m**4 - 12*a**6*(-6 + 545*b_fR)*(-1 + Omega_m)**2*Omega_m**5 + a**3*(-78 + &
            37*b_fR)*(-1 + Omega_m)*Omega_m**6 + 6*Omega_m**7))/(4*a**3*(-1 + Omega_m) - Omega_m)**9

    Else

       derivative_dark_energy_density = 0.d0

    End if

  end function derivative_dark_energy_density

  !####################################
  ! DERIVATIVE DARK ENERGY DENSITY ENDS
  !####################################

  !###################################
  !DERIVATIVE EQUATION OF STATE STARTS
  !###################################

  function derivative_equation_of_state(a)

    Implicit none

    Real*8 :: a,derivative_equation_of_state

    If (MG_parametrisation .eq. 'HS_Basilakos') then

       derivative_equation_of_state = 12*a**2*b_fR*(-1 + Omega_m)*Omega_m*((-24*a**3*(a**3*(-1 + Omega_m) - &
            Omega_m)*(-1 + Omega_m))/(-4*a**3*(-1 + Omega_m) + Omega_m)**4 - (48*a**3*(a**3*(-1 + Omega_m) - &
            Omega_m)*(-1 + Omega_m)*(8*a**3*(-1 + Omega_m) + Omega_m))/(-4*a**3*(-1 + Omega_m) + Omega_m)**5 - &
            (3*(a**3*(-1 + Omega_m) - Omega_m)*(8*a**3*(-1 + Omega_m) + Omega_m))/(-4*a**3*(-1 + Omega_m) + Omega_m)**4 - &
            (3*a**3*(-1 + Omega_m)*(8*a**3*(-1 + Omega_m) + Omega_m))/(-4*a**3*(-1 + Omega_m) + Omega_m)**4 - &
            (12*a**6*b_fR*(-1 + Omega_m)**2*(10240*a**15*(-1 + Omega_m)**5 - 16640*a**12*(-1 + Omega_m)**4*Omega_m + &
            12544*a**9*(-1 + Omega_m)**3*Omega_m**2 + 2306*a**6*(-1 + Omega_m)**2*Omega_m**3 - &
            5419*a**3*(-1 + Omega_m)*Omega_m**4 - 277*Omega_m**5))/(4*a**3*(-1 + Omega_m) - Omega_m)**9 -& 
            (a**3*b_fR*(-1 + Omega_m)*(40960*a**18*(-1 + Omega_m)**6 - 79872*a**15*(-1 + Omega_m)**5*Omega_m + &
            75264*a**12*(-1 + Omega_m)**4*Omega_m**2 + 18448*a**9*(-1 + Omega_m)**3*Omega_m**3 - &
            65028*a**6*(-1 + Omega_m)**2*Omega_m**4 - 6648*a**3*(-1 + Omega_m)*Omega_m**5 + 109*Omega_m**6))/(4*a**3*(-1 + &
            Omega_m) - Omega_m)**9 + (18*a**6*b_fR*(-1 + Omega_m)**2*(40960*a**18*(-1 + Omega_m)**6 - 79872*a**15*(-1 + &
            Omega_m)**5*Omega_m + 75264*a**12*(-1 + Omega_m)**4*Omega_m**2 + 18448*a**9*(-1 + Omega_m)**3*Omega_m**3 - &
            65028*a**6*(-1 + Omega_m)**2*Omega_m**4 - 6648*a**3*(-1 + Omega_m)*Omega_m**5 + 109*Omega_m**6))/(-4*a**3*(-1 + &
            Omega_m) + Omega_m)**10)

    Else

       derivative_equation_of_state = 0.d0

    End if

  end function derivative_equation_of_state

  !#################################
  !DERIVATIVE EQUATION OF STATE ENDS
  !#################################

  !###################
  ! F = 1 + f_R STARTS
  !###################

  function F_MG(a)

    Implicit none 

    Real*8 :: F_MG,a

       F_MG = 1.d0 + fMG_R(a)

  end function F_MG

  !#################
  ! F = 1 + f_R ENDS
  !#################

  !##########
  ! F' STARTS
  !##########

  function F_MG_prime(a)

    Implicit none 

    Real*8 :: F_MG_prime,a

    F_MG_prime = fMG_R_prime(a) 

  end function F_MG_prime

  !########
  ! F' ENDS
  !########

  !##########
  !F'' STARTS
  !##########

  function F_MG_double_prime(a)

    Implicit none

    Real*8 :: F_MG_double_prime,a

    F_MG_double_prime = fMG_R_double_prime(a) 

  end function F_MG_double_prime

  !########
  !F'' ENDS
  !########

  !###########
  ! F_R STARTS
  !###########

  function FR(a)

    Implicit none

    Real*8 :: a,FR

    FR = fMG_RR(a)

  end function FR

  !#########
  ! F_R ENDS
  !#########

  !#################
  ! F_R_prime STARTS
  !#################

  function FR_prime(a)

    Implicit none 

    Real*8 :: FR_prime,a

    FR_prime = fMG_RR_prime(a)

  end function FR_prime

  function FR_prime_fgsl(a, params) bind(c)
    
    real(c_double), value :: a
    type(c_ptr), value :: params
    real(c_double) :: FR_prime_fgsl

    FR_prime_fgsl = FR_prime(a)

  end function FR_prime_fgsl

  !###############
  ! F_R_prime ENDS
  !###############

  !########################
  ! F_R_double_prime STARTS
  !########################

  function FR_double_prime(a)

    real(fgsl_double) :: result, abserr
    integer(fgsl_int) :: status
    type(fgsl_function) :: pwr
    real(fgsl_double) :: FR_double_prime,a

    pwr = fgsl_function_init(FR_prime_fgsl, c_null_ptr)

    status = fgsl_deriv_central (pwr, a , 1.E-8_fgsl_double, &
         result, abserr)

    FR_double_prime = result

  end function FR_double_prime

  !######################
  ! F_R_double_prime ENDS
  !######################
  
  !#################
  ! G_eff/G_N STARTS
  !#################

  function Geff_over_GN(a)

    Implicit none

    Real*8 :: Geff_over_GN,a

    Geff_over_GN = (1.d0 + 4.d0*wavenumber_k**2*FR(a)/a**2/F_MG(a))/&
         (F_MG(a) + 3.d0*wavenumber_k**2*FR(a)/a**2)

  end function Geff_over_GN

  !###############
  ! G_eff/G_N ENDS
  !###############

  !#############
  ! Q_eff STARTS
  !#############

  function Qeff(a)

    Implicit none

    Real*8 :: Qeff,a

    Qeff = (1.d0 + 2.d0*wavenumber_k**2*FR(a)/a**2/F_MG(a))/&
         (F_MG(a) + 3.d0*wavenumber_k**2*FR(a)/a**2)

  end function Qeff

  !###########
  ! Q_eff ENDS
  !###########

  !################
  ! \Omega_m STARTS
  !################

  function Omega_Matter(a)

    Implicit none

    Real*8 :: a, Omega_Matter

    Omega_Matter = Omega_m/a**3

  end function Omega_Matter

  !##############
  ! \Omega_m ENDS
  !##############

  !#################
  ! \Omega_DE STARTS
  !#################

  function Omega_DE(a)

    Implicit none

    Real*8 :: a, Omega_DE

    Omega_DE = dark_energy_density(a)/3.d0/H0**2 

  end function Omega_DE

  function Omega_DE_fgsl(x, params) bind(c)
    
    real(c_double), value :: x
    type(c_ptr), value :: params
    real(c_double) :: Omega_DE_fgsl

    Omega_DE_fgsl = Omega_DE(x)

  end function Omega_DE_fgsl

  !###############
  ! \Omega_DE ENDS
  !###############

  !######################
  ! \Omega_m_prime STARTS
  !######################

  function derivative_Omega_Matter(a)

    Implicit none

    Real*8 :: a, derivative_Omega_Matter

    derivative_Omega_Matter = -3.d0*Omega_m/a**4

  end function derivative_Omega_Matter

  !####################
  ! \Omega_m_prime ENDS
  !####################


  !#######################
  ! \Omega_DE_prime STARTS
  !#######################

  function derivative_Omega_DE(x)

    real(fgsl_double) :: result, abserr
    integer(fgsl_int) :: status
    type(fgsl_function) :: pwr
    real(fgsl_double) :: derivative_Omega_DE,x

    If ( ( (MG_parametrisation .eq. 'GR_DE') .or. (MG_parametrisation .eq. 'GR_LAMBDA') ) .or. &
         (MG_parametrisation .eq. 'Savvas') ) then

       derivative_Omega_DE = -3.d0*(1.d0 - Omega_m)*(1.d0 + equation_of_state(x))*&
            x**(-4.d0 - 3.d0*equation_of_state(x))

    Else

       pwr = fgsl_function_init(Omega_DE_fgsl, c_null_ptr)

       status = fgsl_deriv_central (pwr, x , 1.E-7_fgsl_double, &
            result, abserr)

       derivative_Omega_DE = result

    End if

  end function derivative_Omega_DE

  !#####################
  ! \Omega_DE_prime ENDS
  !#####################

  !############################################################
  ! SUBROUTINE TO COMPUTE AND PLOT BACKGROUND QUANTITIES STARTS
  !############################################################

  subroutine compute_background()

    Implicit none

    Integer*4,parameter :: number_points = 1000
    Integer*4 :: index
    Real*8,dimension(number_points) :: scale_factor

    open(UNIT_TEST,file='./output/background_functions.txt')

    write(UNIT_TEST,*) '# a  conformal_Hubble_parameter(a)  derivative_conformal_Hubble_parameter(a)'&
         '  second_derivative_conformal_Hubble_parameter(a)  ricci_scalar(a)  ricci_scalar_prime(a)'&
         '  f(R(a))  f_R(R(a))  f_R_prime(R(a))  f_R_double_prime(R(a))  DE_density(a)  w_DE(a)'&
         '  f_RR(R(a))  f_RR_prime(R(a))  f_RR_double_prime(R(a))  Omega_M(a)  Omega_DE(a) '&
         '  derivative_Omega_DE(a)  f_R-k**2/a**2*(1+3*f_R)*f_RR/F  third_derivative_conformal_Hubble'

    Do index=1,number_points

       scale_factor(index) = 10**(log10(initial_scale_factor) + real(index-1)*(log10(final_scale_factor) - &
            log10(initial_scale_factor))/real(number_points-1))

       write(UNIT_TEST,99) scale_factor(index),&
            conformal_Hubble_parameter(scale_factor(index)),&
            derivative_conformal_Hubble_parameter(scale_factor(index)),&
            second_derivative_conformal_Hubble_parameter(scale_factor(index)),&
            ricci_scalar(scale_factor(index)),&
            ricci_scalar_prime(scale_factor(index)),&
            fMG(scale_factor(index)),&
            fMG_R(scale_factor(index)),&
            FMG_R_prime(scale_factor(index)),&
            fMG_R_double_prime(scale_factor(index)),&
            dark_energy_density(scale_factor(index)),&
            equation_of_state(scale_factor(index)),&
            fMG_RR(scale_factor(index)),&
            fMG_RR_prime(scale_factor(index)),&
            fMG_RR_double_prime(scale_factor(index)),&
            Omega_Matter(scale_factor(index)),&
            Omega_DE(scale_factor(index)),&
            derivative_Omega_DE(scale_factor(index)),&
!            3.d0*scale_factor(index)**2*F_MG(scale_factor(index))*fMG_R(scale_factor(index))
!!$            (3.d0*scale_factor(index)**2*F_MG(scale_factor(index)) - &
!!$            3.d0*scale_factor(index)**2*F_MG(scale_factor(index))**2 + &
            3.d0*wavenumber_k**2*(1.d0 + 3.d0*fMG_R(scale_factor(index)))*fMG_RR(scale_factor(index)),&
            third_derivative_conformal_Hubble_parameter(scale_factor(index))

99     FORMAT(ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,&
            ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10)

    End do

    close(UNIT_TEST)

    If (MG_parametrisation .eq. 'GR_DE') then

       write(UNIT_EXE_FILE,*) 'PLOTTING BACKGROUND'

       call system('cd figures; python plot_background_GR_DE.py')

    Else if (MG_parametrisation .eq. 'GR_LAMBDA') then

       write(UNIT_EXE_FILE,*) 'PLOTTING BACKGROUND'

       call system('cd figures; python plot_background_GR_LAMBDA.py')
   
    Else if (MG_parametrisation .eq. 'Savvas') then

       write(UNIT_EXE_FILE,*) 'PLOTTING BACKGROUND'

       call system('cd figures; python plot_background_Savvas.py')

    Else if (MG_parametrisation .eq. 'HS_Basilakos') then

       write(UNIT_EXE_FILE,*) 'PLOTTING BACKGROUND'

       call system('cd figures; python plot_background.py')

    Else if (MG_parametrisation .eq. 'Starobinsky_Basilakos') then

       write(UNIT_EXE_FILE,*) 'MG_PARAMETRISATION Starobinsky_Basilakos IS NOT YET IMPLEMENTED'

       stop

    Else 

       write(UNIT_EXE_FILE,*) 'UNKNOWN MG_PARAMETRISATION'

       stop

    End if 

  end subroutine compute_background

  !##########################################################
  ! SUBROUTINE TO COMPUTE AND PLOT BACKGROUND QUANTITIES ENDS
  !##########################################################

end Module background
