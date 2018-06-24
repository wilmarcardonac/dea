Module background
  
  use fgsl
  use, intrinsic :: iso_c_binding
  use fiducial

  Implicit none

Contains

  ! CONFORMAL HUBBLE PARAMETER

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

  ! DERIVATIVE OF THE CONFORMAL HUBBLE PARAMETER 

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

  ! SECOND DERIVATIVE OF THE CONFORMAL HUBBLE PARAMETER 

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

  ! RICCI SCALAR 

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

  ! DERIVATIVE OF THE RICCI SCALAR

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

  ! f(R(a))

  function fMG(a)

    Implicit none

    Real*8 :: fMG,a,R,x,y,R0,x0,y0
    Real*8,parameter :: b2 = (-7.d0+sqrt(73.d0))/12.d0
    Real*8,parameter :: a2 = 1.d0 + b2
    Real*8,parameter :: b = 1.5d0 + b2
    Real*8,parameter :: c = 13.d0/6.d0 + 2.d0*b2
    Real*8 :: alpha   ! PARAMETER IN THE f(R) PARAMETRISATION BY SAVVAS. 

    R = ricci_scalar(a)

    R0 = ricci_scalar(1.d0)

    x = Lambda/(R-3.d0*Lambda)

    x0 = Lambda/(R0-3.d0*Lambda)

    y = R/(R-3.d0*Lambda)

    y0 = R0/(R0-3.d0*Lambda)

    alpha = -fR0/(b2*H0**2/R0*y0**a2*(Lambda/R0)**b2*fgsl_sf_hyperg_2f1(a2,b,c,x0)) 

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

  ! f_R(R(a))

  function fMG_R(a)

    Implicit none

    Real*8 :: fMG_R,a,R,x,y,R0,x0,y0
    Real*8,parameter :: b2 = (-7.d0+sqrt(73.d0))/12.d0
    Real*8,parameter :: a2 = 1.d0 + b2
    Real*8,parameter :: b = 1.5d0 + b2
    Real*8,parameter :: c = 13.d0/6.d0 + 2.d0*b2
    Real*8 :: alpha   ! PARAMETER IN THE f(R) PARAMETRISATION BY SAVVAS. 

    R = ricci_scalar(a)

    R0 = ricci_scalar(1.d0)

    x = Lambda/(R-3.d0*Lambda)

    x0 = Lambda/(R0-3.d0*Lambda)

    y = R/(R-3.d0*Lambda)

    y0 = R0/(R0-3.d0*Lambda)

    alpha = -fR0/(b2*H0**2/R0*y0**a2*(Lambda/R0)**b2*fgsl_sf_hyperg_2f1(a2,b,c,x0)) 

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

  ! f_R_prime(R(a))

  function fMG_R_prime(a)

    Implicit none

    Real*8 :: fMG_R_prime,a,R,x,y,R0,x0,y0,R_prime
    Real*8,parameter :: b2 = (-7.d0+sqrt(73.d0))/12.d0
    Real*8,parameter :: a2 = 1.d0 + b2
    Real*8,parameter :: a3 = 2.d0 + b2
    Real*8,parameter :: b = 1.5d0 + b2
    Real*8,parameter :: c = 13.d0/6.d0 + 2.d0*b2
    Real*8 :: alpha   ! PARAMETER IN THE f(R) PARAMETRISATION BY SAVVAS. 

    R = ricci_scalar(a)

    R0 = ricci_scalar(1.d0)

    R_prime = ricci_scalar_prime(a)

    x = Lambda/(R-3.d0*Lambda)

    x0 = Lambda/(R0-3.d0*Lambda)

    y = R/(R-3.d0*Lambda)

    y0 = R0/(R0-3.d0*Lambda)

    alpha = -fR0/(b2*H0**2/R0*y0**a2*(Lambda/R0)**b2*fgsl_sf_hyperg_2f1(a2,b,c,x0)) 

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

  ! f_R_double_prime(R(a))

  function fMG_R_double_prime(a)

    real(fgsl_double) :: result, abserr
    integer(fgsl_int) :: status
    type(fgsl_function) :: pwr
    real(fgsl_double) :: fMG_R_double_prime,a

    pwr = fgsl_function_init(fMG_R_prime_fgsl, c_null_ptr)

    status = fgsl_deriv_central (pwr, a , 1.E-6_fgsl_double, &
         result, abserr)

    fMG_R_double_prime = result

  end function fMG_R_double_prime

  ! DARK ENERGY DENSITY FOR AN EFFECTIVE FLUID FROM F(R) PARAMMETRISATION (IT INCLUDES COSMOLOGICAL CONSTANT)

  function dark_energy_density(a) ! ACTUALLY THIS IS THE DARK ENERGY DENSITY MULTIPLIED BY \kappa

    Implicit none

    Real*8 :: a,dark_energy_density,R,H
    
    R = ricci_scalar(a)

    H = conformal_Hubble_parameter(a) 

    If (MG_parametrisation .eq. 'GR_DE') then

       dark_energy_density = 3.d0*H0**2*(1.d0 - Omega_m)*a**(-3.d0 - 3.d0*w0_fld)

    Else

       dark_energy_density = (fMG_R(a)*R - fMG(a))/2.d0 - 3.d0*H**2*fMG_R(a)/a**2 - &
            3.d0*H**2*fMG_R_prime(a)/a

    End if

  end function dark_energy_density

  ! EQUATION OF STATE FOR AN EFFECTIVE FLUID FROM F(R) PARAMMETRISATION BY BASILAKOS ET AL.

  function equation_of_state(a)

    Implicit none

    Real*8 :: a,equation_of_state,H,H_prime,R

    H = conformal_Hubble_parameter(a)

    H_prime = derivative_conformal_Hubble_parameter(a)

    R = ricci_scalar(a)

    If (MG_parametrisation .eq. 'GR_DE') then

       equation_of_state = w0_fld

    Else

       equation_of_state = -1.d0/3.d0 -2.d0/3.d0*(H**2*fMG_R(a) - a**2*fMG(a)/6.d0 - H**2*a*fMG_R_prime(a) - &
            fMG_R_double_prime(a)*H**2*a**2/2.d0 - fMG_R_prime(a)*H_prime*H*a**2/2.d0)/(-H**2*fMG_R(a) - &
            a**2*fMG(a)/6.d0 - H**2*a*fMG_R_prime(a) + a**2*fMG_R(a)*R/6.d0 ) 

    End if

  end function equation_of_state



  
  ! DERIVATIVE DARK ENERGY DENSITY FOR AN EFFECTIVE FLUID FROM F(R) PARAMMETRISATION BY BASILAKOS ET AL.

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

  ! DARK ENERGY PRESSURE PERTURBATION OVER DARK ENERGY DENSITY AND OTHER FUNCTIONS FOR AN EFFECTIVE FLUID FROM F(R) PARAMETRISATION BY BASILAKOS ET AL. 


  function F_MG(a)

    Implicit none 

    Real*8 :: F_MG,a

    If (MG_parametrisation .eq. 'Savvas') then

       F_MG = F_MG_savvas(a)

    Else

       F_MG = 1 - (2*a**6*b_fR*(-1 + Omega_m)**2*(-4*a**3*(-1 + Omega_m) + &
            Omega_m)**18)/(-16384*a**30*(-64 + 8*b_fR + b_fR**2)*(-1 + Omega_m)**10 - 8192*a**27*(320 - &
            24*b_fR + b_fR**2)*(-1 + Omega_m)**9*Omega_m + 1536*a**24*(1920 - 40*b_fR + 31*b_fR**2)*(-1 + &
            Omega_m)**8*Omega_m**2 + 128*a**21*(-15360 - 1080*b_fR + 43*b_fR**2)*(-1 + Omega_m)**7*Omega_m**3 + &
            64*a**18*(13440 + 2496*b_fR + 3817*b_fR**2)*(-1 + Omega_m)**6*Omega_m**4 - 384*a**15*(672 + &
            193*b_fR + 844*b_fR**2)*(-1 + Omega_m)**5*Omega_m**5 - 8*a**12*(-6720 - 2226*b_fR + 7307*b_fR**2)*(-1 + &
            Omega_m)**4*Omega_m**6 + 10*a**9*(-768 - 222*b_fR + 37*b_fR**2)*(-1 + Omega_m)**3*Omega_m**7 + &
            120 *a**6*(6 + b_fR)*(-1 + Omega_m)**2*Omega_m**8 - a**3*(40 + b_fR)*(-1 + Omega_m)*Omega_m**9 + Omega_m**10)**2
 
    End if

  end function F_MG

  function F_MG_savvas(a)

    Implicit none

    Real*8,parameter :: b2 = (-7.d0+sqrt(73.d0))/12.d0
    Real*8 :: F_MG_savvas,a,x,y
    Real*8,parameter :: a2 = 1.d0 + b2
    Real*8,parameter :: b = 1.5d0 + b2
    Real*8,parameter :: c = 13.d0/6.d0 + 2.d0*b2

    Real*8,parameter :: alpha = 1.d0

    x = Lambda/(ricci_scalar(a)-3.d0*Lambda)

    y = ricci_scalar(a)/(ricci_scalar(a)-3.d0*Lambda)

    F_MG_savvas = 1.d0 - b2*H0**2*alpha/ricci_scalar(a)*y**a2*(Lambda/ricci_scalar(a))**b2*&
         fgsl_sf_hyperg_2f1(a2,b,c,x)

  end function F_MG_savvas

  function F_MG_savvas_prime(a)

    Implicit none

    Real*8,parameter :: b2 = (-7.d0+sqrt(73.d0))/12.d0
    Real*8 :: F_MG_savvas_prime,a,x
    Real*8,parameter :: a2 = 1.d0 + b2
    Real*8,parameter :: b = 1.5d0 + b2
    Real*8,parameter :: c = 13.d0/6.d0 + 2.d0*b2
    Real*8,parameter :: a3 = 2.d0 + b2

    Real*8,parameter :: alpha = 1.d0

    x = Lambda/(ricci_scalar(a)-3.d0*Lambda)

    F_MG_savvas_prime = b2*a2*H0**2*alpha*(Lambda/ricci_scalar(a))**b2*&
         (1.d0 + 3.d0*x)**b2*ricci_scalar_prime(a)/( ricci_scalar(a) - 3.d0*Lambda )**2*&
         fgsl_sf_hyperg_2f1(b,a3,c,x)

  end function F_MG_savvas_prime

  function F_MG_prime(a)

    Implicit none 

    Real*8 :: F_MG_prime,a

    If (MG_parametrisation .eq. 'Savvas') then

       F_MG_prime = F_MG_savvas_prime(a)

    Else

       F_MG_prime = (-8*a**3*b_fR*H0**4*(-1 + Omega_m)**2*(-2*conformal_Hubble_parameter(a)**2 + &
            a**2*derivative_conformal_Hubble_parameter(a)**2 + &
            a*conformal_Hubble_parameter(a)*(a*second_derivative_conformal_Hubble_parameter(a) + &
            derivative_conformal_Hubble_parameter(a))))/(a**2*b_fR*H0**2*(-1 + Omega_m) - &
            2*conformal_Hubble_parameter(a)**2 - &
            2*a*conformal_Hubble_parameter(a)*derivative_conformal_Hubble_parameter(a))**3 

    End if 

  end function F_MG_prime

  function F_MG_prime_fgsl(a, params) bind(c)
    
    real(c_double), value :: a
    type(c_ptr), value :: params
    real(c_double) :: F_MG_prime_fgsl

    F_MG_prime_fgsl = F_MG_prime(a)

  end function F_MG_prime_fgsl

  function F_MG_double_prime(a)

    real(fgsl_double) :: result, abserr
    integer(fgsl_int) :: status
    type(fgsl_function) :: pwr
    real(fgsl_double) :: F_MG_double_prime,a

    pwr = fgsl_function_init(F_MG_prime_fgsl, c_null_ptr)

    status = fgsl_deriv_central (pwr, a , 1.E-8_fgsl_double, &
         result, abserr)

    F_MG_double_prime = result

  end function F_MG_double_prime

  function FR(a)

    Implicit none

    Real*8 :: a,FR

    FR = (4*a**9*b_fR*(4*a**3*(-1 + Omega_m) - Omega_m)**27*(-1 + Omega_m)**2)/(3.*H0**2*(16384*a**30*(-64 + &
         8*b_fR + b_fR**2)*(-1 + Omega_m)**10 + 8192*a**27*(320 - 24*b_fR + b_fR**2)*(-1 + Omega_m)**9*Omega_m - &
         1536*a**24*(1920 - 40*b_fR + 31*b_fR**2)*(-1 + Omega_m)**8*Omega_m**2 - 128*a**21*(-15360 - 1080*b_fR + &
         43*b_fR**2)*(-1 + Omega_m)**7*Omega_m**3 - 64*a**18*(13440 + 2496*b_fR + 3817*b_fR**2)*(-1 + Omega_m)**6*Omega_m**4 + &
         384*a**15*(672 + 193*b_fR + 844*b_fR**2)*(-1 + Omega_m)**5*Omega_m**5 + 8*a**12*(-6720 - 2226*b_fR + 7307*b_fR**2)*(-1 &
         + Omega_m)**4*Omega_m**6 - 10*a**9*(-768 - 222*b_fR + 37*b_fR**2)*(-1 + Omega_m)**3*Omega_m**7 - &
         120*a**6*(6 + b_fR)*(-1 + Omega_m)**2*Omega_m**8 + a**3*(40 + b_fR)*(-1 + Omega_m)*Omega_m**9 - Omega_m**10)**3)

  end function FR

  function FR_prime(a)

    Implicit none 

    Real*8 :: FR_prime,a

    FR_prime = (12*a**8*b_fR*(-1 + Omega_m)**2*Omega_m*(-4*a**3*(-1 + Omega_m) + &
         Omega_m)**26*(16384*a**30*(64 + 24*b_fR + 11*b_fR**2)*(-1 + Omega_m)**10 - 4096*a**27*(640 + 264*b_fR + &
         77*b_fR**2)*(-1 + Omega_m)**9*Omega_m - 30720*a**24*(-96 - 68*b_fR + 13*b_fR**2)*(-1 + Omega_m)**8*Omega_m**2 - &
         256*a**21*(7680 + 6744*b_fR + 15397*b_fR**2)*(-1 + Omega_m)**7*Omega_m**3 + 320*a**18*(2688 + 2136*b_fR + &
         16439*b_fR**2)*(-1 + Omega_m)**6*Omega_m**4 + 192*a**15*(-1344 - 682*b_fR + 14059*b_fR**2)*(-1 + Omega_m)**5*Omega_m**5 + &
         16*a**12*(3360 + 546*b_fR + 10313*b_fR**2)*(-1 + Omega_m)**4*Omega_m**6 - 20*a**9*(384 - 30*b_fR + 37*b_fR**2)*(-1 + &
         Omega_m)**3*Omega_m**7 - 12*a**6*(-60 + 7*b_fR)*(-1 + Omega_m)**2*Omega_m**8 - 40*a**3*(-1 + Omega_m)*Omega_m**9 + &
         Omega_m**10))/(H0**2*(-16384*a**30*(-64 + 8*b_fR + b_fR**2)*(-1 + Omega_m)**10 - 8192*a**27*(320 - 24*b_fR + &
         b_fR**2)*(-1 + Omega_m)**9*Omega_m + 1536*a**24*(1920 - 40*b_fR + 31*b_fR**2)*(-1 + Omega_m)**8*Omega_m**2 + &
         128*a**21*(-15360 - 1080*b_fR + 43*b_fR**2)*(-1 + Omega_m)**7*Omega_m**3 + 64*a**18*(13440 + 2496*b_fR + &
         3817*b_fR**2)*(-1 + Omega_m)**6*Omega_m**4 - 384*a**15*(672 + 193*b_fR + 844*b_fR**2)*(-1 + Omega_m)**5*Omega_m**5 - &
         8*a**12*(-6720 - 2226*b_fR + 7307*b_fR**2)*(-1 + Omega_m)**4*Omega_m**6 + 10*a**9*(-768 - 222*b_fR + &
         37*b_fR**2)*(-1 + Omega_m)**3*Omega_m**7 + 120*a**6*(6 + b_fR)*(-1 + Omega_m)**2*Omega_m**8 - &
         a**3*(40 + b_fR)*(-1 + Omega_m)*Omega_m**9 + Omega_m**10)**4)

  end function FR_prime

  function FR_prime_fgsl(a, params) bind(c)
    
    real(c_double), value :: a
    type(c_ptr), value :: params
    real(c_double) :: FR_prime_fgsl

    FR_prime_fgsl = FR_prime(a)

  end function FR_prime_fgsl

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


  function Geff_over_GN(a)

    Implicit none

    Real*8 :: Geff_over_GN,a

    Geff_over_GN = (1.d0 + 4.d0*wavenumber_k**2*FR(a)/a**2/F_MG(a))/&
         (F_MG(a) + 3.d0*wavenumber_k**2*FR(a)/a**2)

  end function Geff_over_GN

  function Qeff(a)

    Implicit none

    Real*8 :: Qeff,a

    Qeff = (1.d0 + 2.d0*wavenumber_k**2*FR(a)/a**2/F_MG(a))/&
         (F_MG(a) + 3.d0*wavenumber_k**2*FR(a)/a**2)

  end function Qeff

  function Omega_DE(a)

    Implicit none

    Real*8 :: a, Omega_DE

    Omega_DE = 1.d0 - Omega_m + (2*a**2*b_fR*(-1 + Omega_m)**2*(12*a**7*(-1 + Omega_m)**2 + 3*a**4*(-1 + Omega_m)*Omega_m - &
         6*a*Omega_m**2))/(4*a**3*(-1 + Omega_m) - Omega_m)**3 - (a**5*b_fR**2*(-1 + Omega_m)**3*(-1024*a**19*(-1 + &
         Omega_m)**6 - 9216*a**16*(-1 + Omega_m)**5*Omega_m + 22848*a**13*(-1 + Omega_m)**4*Omega_m**2 - 25408*a**10*(-1 + &
         Omega_m)**3*Omega_m**3 + 7452*a**7*(-1 + Omega_m)**2*Omega_m**4 + 4656*a**4*(-1 + Omega_m)*Omega_m**5 - &
         37*a*Omega_m**6))/(-4*a**3*(-1 + Omega_m) + Omega_m)**8

  end function Omega_DE

  function Omega_DE_fgsl(x, params) bind(c)
    
    real(c_double), value :: x
    type(c_ptr), value :: params
    real(c_double) :: Omega_DE_fgsl

    Omega_DE_fgsl = Omega_DE(x)

  end function Omega_DE_fgsl

  function derivative_Omega_DE(x)

    real(fgsl_double) :: result, abserr
    integer(fgsl_int) :: status
    type(fgsl_function) :: pwr
    real(fgsl_double) :: derivative_Omega_DE,x


    pwr = fgsl_function_init(Omega_DE_fgsl, c_null_ptr)

    status = fgsl_deriv_central (pwr, x , 1.E-8_fgsl_double, &
         result, abserr)

    derivative_Omega_DE = result

  end function derivative_Omega_DE

  subroutine compute_background()

    Implicit none

    Integer*4,parameter :: number_points = 1000
    Integer*4 :: index
    Real*8,dimension(number_points) :: scale_factor

    open(UNIT_TEST,file='./output/background_functions.txt')

    write(UNIT_TEST,*) '# a  conformal_Hubble_parameter(a)  derivative_conformal_Hubble_parameter(a)'&
         '  second_derivative_conformal_Hubble_parameter(a)  ricci_scalar(a)  ricci_scalar_prime(a)'&
         '  f(R(a))  f_R(R(a))  f_R_prime(R(a))  f_R_double_prime(R(a))  DE_density(a)  w_DE(a)'
  !w_prime(a)  H(a)/H0  cs2(a)  ceff2(a)  Omega_m(a)  Omega_DE(a)'&
!         '  pressure_perturbation_over_density  dm_th  Vm_th  dde_th  Vde_th  \pi(a)  Geff/GN  Qeff'&
!         '  ca2(a) x'

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
            equation_of_state(scale_factor(index))

99     FORMAT(ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10)!,&
!            ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10)!,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,&
!            ES20.10,ES20.10,ES20.10,ES20.10)


!            derivative_equation_of_state(scale_factor(index)),&

!            sound_speed_squared(scale_factor(index)),&
!            effective_sound_speed_squared(scale_factor(index),wavenumber_k),&
!            Omega_m/scale_factor(index)**3,&
!            Omega_DE(scale_factor(index)),&
!            dark_energy_pressure_perturbation(scale_factor(index)),&
!            matter_density_perturbation(scale_factor(index),wavenumber_k)/delta_0(wavenumber_k),&
!            matter_velocity_perturbation(scale_factor(index),wavenumber_k)/delta_0(wavenumber_k),&
!            dark_energy_density_perturbation(scale_factor(index),&
!            matter_density_perturbation(scale_factor(index),wavenumber_k))/delta_0(wavenumber_k),&
!            dark_energy_velocity_perturbation(scale_factor(index),&
!            matter_density_perturbation(scale_factor(index),wavenumber_k))/delta_0(wavenumber_k),&
!            anisotropic_stress(scale_factor(index),0.d0,0.d0,0.d0,0.d0),&
!            Geff_over_GN(scale_factor(index)),&
!            Qeff(scale_factor(index)),&
!            adiabatic_sound_speed_squared(scale_factor(index)),&
            !F_MG_savvas_prime(scale_factor(index))
!            Lambda/(ricci_scalar(scale_factor(index))-3.d0*Lambda)
            !
!            effective_sound_speed_squared(scale_factor(index),wavenumber_k),speedL,&
!            de_density_perturbation_sub_sound_horizon(wavenumber_k),&
!            (1.d0 + equation_of_state(scale_factor(index)))*&
!            de_velocity_perturbation_sub_sound_horizon(scale_factor(index),wavenumber_k),&
!            dark_energy_density_perturbation_super_sound_horizon(scale_factor(index),wavenumber_k),&
!            (1.d0 + equation_of_state(scale_factor(index)))*&
!            dark_energy_velocity_perturbation_super_sound_horizon(scale_factor(index),wavenumber_k),&
!            anisotropic_stress(scale_factor(index),0.d0,0.d0,0.d0,0.d0)*Omega_m&
!            /( scale_factor(index)*conformal_Hubble_parameter(scale_factor(index))**2/H0**2 - Omega_m )*&
!            matter_density_perturbation(scale_factor(index),wavenumber_k),&
            
!            derivative_anisotropic_stress(scale_factor(index)),&
!           phi(scale_factor(index),wavenumber_k,matter_density_perturbation(scale_factor(index),wavenumber_k),&

!            matter_velocity_perturbation(scale_factor(index),wavenumber_k),&

!!$            psi(scale_factor(index),wavenumber_k,matter_density_perturbation(scale_factor(index),wavenumber_k),&
!!$            dark_energy_density_perturbation(scale_factor(index),&
!!$            matter_density_perturbation(scale_factor(index),wavenumber_k)),&
!!$            matter_velocity_perturbation(scale_factor(index),wavenumber_k),&
!!$            dark_energy_velocity_perturbation(scale_factor(index),&
!!$            matter_density_perturbation(scale_factor(index),wavenumber_k))),&
!!$            -2.d0/3.d0,equation_of_state(scale_factor(index)),-1.d0,&
!!$            wavenumber_k**2*FR(scale_factor(index))/scale_factor(index)**2/F_MG(scale_factor(index)),& !term1
!!$            3.d0*conformal_Hubble_parameter(scale_factor(index))**2*FR_prime(scale_factor(index))/& !term2
!!$            2.d0/scale_factor(index) + 3.d0*conformal_Hubble_parameter(scale_factor(index))*&
!!$            derivative_conformal_Hubble_parameter(scale_factor(index))*FR_prime(scale_factor(index))/2.d0&
!!$            + 3.d0*conformal_Hubble_parameter(scale_factor(index))**2*FR_double_prime(scale_factor(index))/2.d0,&
!!$            scale_factor(index)*FR_prime(scale_factor(index))/F_MG(scale_factor(index))*& !term3
!!$            conformal_Hubble_parameter(scale_factor(index))/2.d0,&
!!$            wavenumber_k**2*conformal_Hubble_parameter(scale_factor(index))/scale_factor(index)*& !term4
!!$            FR_prime(scale_factor(index))/F_MG(scale_factor(index)),&
!!$            2.d0*wavenumber_k**2*FR(scale_factor(index))*F_MG_prime(scale_factor(index))*& !term5
!!$            conformal_Hubble_parameter(scale_factor(index))/F_MG(scale_factor(index))**2/&
!!$            scale_factor(index),&
!!$            1.d0 - F_MG(scale_factor(index)) + wavenumber_k**2*(2.d0 - &  ! denominator
!!$            3.d0*F_MG(scale_factor(index)))*FR(scale_factor(index))/scale_factor(index)**2/&
!!$            F_MG(scale_factor(index)),&
!!$            (effective_sound_speed_squared(scale_factor(index),wavenumber_k) - &     ! fullceff2
!!$            3.d0*scale_factor(index)*conformal_Hubble_parameter(scale_factor(index))**2*&
!!$            derivative_equation_of_state(scale_factor(index))/wavenumber_k**2 + &
!!$            3.d0/wavenumber_k**2*( (1.d0 - 3.d0*equation_of_state(scale_factor(index)))*&
!!$            conformal_Hubble_parameter(scale_factor(index))**2 + scale_factor(index)*&
!!$            conformal_Hubble_parameter(scale_factor(index))*derivative_conformal_Hubble_parameter(&
!!$            scale_factor(index)))*( sound_speed_squared(scale_factor(index)) - &
!!$            equation_of_state(scale_factor(index))  ) + &
!!$            3.d0*scale_factor(index)*conformal_Hubble_parameter(scale_factor(index))**2/&
!!$            wavenumber_k**2*derivative_second_dark_energy_pressure_perturbation(scale_factor(index)))*&
!!$            wavenumber_k**2,&
!!$            3.d0*( scale_factor(index)*conformal_Hubble_parameter(scale_factor(index))**2*& ! term6
!!$            F_MG_prime(scale_factor(index)) + scale_factor(index)**2*&
!!$            conformal_Hubble_parameter(scale_factor(index))*&
!!$            derivative_conformal_Hubble_parameter(scale_factor(index))*F_MG_prime(scale_factor(index)) + &
!!$            scale_factor(index)**2*conformal_Hubble_parameter(scale_factor(index))**2*&
!!$            F_MG_double_prime(scale_factor(index)))/2.d0/wavenumber_k**2,&
!!$            6.d0*FR(scale_factor(index))*( conformal_Hubble_parameter(scale_factor(index))**2*& !term7
!!$            F_MG_prime(scale_factor(index))/scale_factor(index) + &
!!$            conformal_Hubble_parameter(scale_factor(index))*&
!!$            derivative_conformal_Hubble_parameter(scale_factor(index))*F_MG_prime(scale_factor(index)) + &
!!$            conformal_Hubble_parameter(scale_factor(index))**2*F_MG_double_prime(scale_factor(index)) )/&
!!$            F_MG(scale_factor(index)),&
!!$            derivative_equation_of_state(scale_factor(index))

            !conformal_Hubble_parameter(scale_factor(index))**2/2.d0/wavenumber_k**2/FR(scale_factor(index)),&

!            phi(scale_factor(index),wavenumber_k,matter_density_perturbation(scale_factor(index),wavenumber_k),&
!            0.d0,0.d0,0.d0),&
!            psi(scale_factor(index),wavenumber_k,matter_density_perturbation(scale_factor(index),wavenumber_k),&
!            0.d0,0.d0,0.d0)

!            derivative_equation_of_state(scale_factor(index)),adiabatic_sound_speed_squared(scale_factor(index)),&
            !dark_energy_density_perturbation(scale_factor(index),matter_density_perturbation(scale_factor(index),&
            !wavenumber_k))*dark_energy_density(scale_factor(index)),&
            !dark_energy_pressure_perturbation_over_dark_energy_density(scale_factor(index),&

!            anisotropic_stress(scale_factor(index),0.d0,1.d0,0.d0,0.d0),
            !fMG(scale_factor(index))
            !F_MG(scale_factor(index))
            !F_MG_prime(scale_factor(index))
            !FR(scale_factor(index))
            !FR_prime(scale_factor(index))
            !FR_double_prime(scale_factor(index))
            !conformal_Hubble_parameter(scale_factor(index))
            !
            !matter_density_perturbation(scale_factor(index),wavenumber_k),&
            !matter_velocity_perturbation(scale_factor(index),wavenumber_k)
       
            !anisotropic_stress(scale_factor(index),matter_density_perturbation(scale_factor(index),wavenumber_k),&
            !0.d0,0.d0,0.d0)/&
            !dark_energy_velocity_perturbation(scale_factor(index),matter_density_perturbation(scale_factor(index),wavenumber_k))
            !3.d0*H0**2*(1.d0-Omega_m)
            !derivative_dark_energy_density(scale_factor(index))
            !dark_energy_density_perturbation(scale_factor(index),matter_density_perturbation(scale_factor(index),wavenumber_k)),&
            !dark_energy_density_perturbation(scale_factor(1),matter_density_perturbation(scale_factor(1),wavenumber_k))*&
            !(scale_factor(index)/scale_factor(1))**(3.75d0) + &
            !2.5d-45*dark_energy_density_perturbation(scale_factor(1),matter_density_perturbation(scale_factor(1),wavenumber_k))*&
            !(scale_factor(index)/scale_factor(1))**(-11.25d0)

            !dark_energy_density_perturbation_sub_sound_horizon(scale_factor(index),wavenumber_k),&
            !-dark_energy_velocity_perturbation_sub_sound_horizon(scale_factor(index),wavenumber_k)*&
            !conformal_Hubble_parameter(scale_factor(index))

            !FR(scale_factor(index))/F_MG(scale_factor(index))/(1.d0 - F_MG(scale_factor(index))),&
            !scale_factor(index)**2/conformal_Hubble_parameter(scale_factor(index))**2/g_pi**2*f_pi

            !(2.d0 - 3.d0*F_MG(scale_factor(index)))*FR(scale_factor(index))/F_MG(scale_factor(index))/(1.d0 - &
            !F_MG(scale_factor(index))),&
            !scale_factor(index)**2/conformal_Hubble_parameter(scale_factor(index))**2/g_pi**2

            !1.d0,&
            !3.d0*scale_factor(index)*F_MG(scale_factor(index))*FR_prime(scale_factor(index))*&
            !conformal_Hubble_parameter(scale_factor(index))**2/2.d0/wavenumber_k**2/FR(scale_factor(index)),&
            !3.d0*scale_factor(index)**2*F_MG(scale_factor(index))*FR_prime(scale_factor(index))*&
            !conformal_Hubble_parameter(scale_factor(index))*derivative_conformal_Hubble_parameter(scale_factor(index))/&
            !2.d0/wavenumber_k**2/FR(scale_factor(index)),&
            !3.d0*scale_factor(index)**2*F_MG(scale_factor(index))*FR_double_prime(scale_factor(index))*&
            !conformal_Hubble_parameter(scale_factor(index))**2/&
            !2.d0/wavenumber_k**2/FR(scale_factor(index))

            !1.d0,F_MG(scale_factor(index)),wavenumber_k**2*(2.d0-3.d0*F_MG(scale_factor(index)))*FR(scale_factor(index))/&
            !scale_factor(index)**2/F_MG(scale_factor(index)),1.d0

            !sound_speed_squared(scale_factor(index))*3.d0/2.d0*&
            !dark_energy_density_perturbation(scale_factor(index),matter_density_perturbation(scale_factor(index),&
            !wavenumber_k))

    End do

    close(UNIT_TEST)

    call system('cd figures; python plot_background.py')

  end subroutine compute_background

end Module background
