Module functions
  
  use fgsl
  use, intrinsic :: iso_c_binding
  
  Implicit none

Contains

  !##########################################################################################################################################################
  !INITIAL CONDITIONS DURING MATTER DOMINANCE. THE EXPRESSIONS MATCH THOSE FOUND BY MARTIN AND DOMENICO (WITHOUT ANISOTROPIC STRESS) IN PRD 80, 083519 (2009)
  !##########################################################################################################################################################
  
  ! MATTER DENSITY PERTURBATIONS. THE RELATION BETWEEN THE CONSTANT 'B' AND \delta_0 IN DOMENICO'S PAPER IS AS FOLLOWS: \delta_0 = -2 B k^2/3/H0^2/\Omega_m 
  
  function matter_density_perturbation(a,k)

    use fiducial

    Implicit none

    Real*8 ::  a,k,B,matter_density_perturbation

    B = initial_condition_gravitational_potential(k)

    matter_density_perturbation = delta_0(k)*(a + 3.d0*H0**2*Omega_m/k**2 )

  end function matter_density_perturbation

  ! MATTER VELOCITY PERTURBATION 

  function matter_velocity_perturbation(a,k)

    use fiducial

    Implicit none

    Real*8 :: a,k,B,matter_velocity_perturbation

    B = initial_condition_gravitational_potential(k)

    matter_velocity_perturbation = -delta_0(k)*sqrt(a)*H0*Sqrt(Omega_m)

  end function matter_velocity_perturbation

  ! DARK ENERGY DENSITY PERTURBATIONS ON SUPER-SOUND HORIZON SCALES 

  function dark_energy_density_perturbation_super_sound_horizon(a,k) 

    use fiducial

    Implicit none 

    Real*8 :: a,k,B,dark_energy_density_perturbation_super_sound_horizon

    B = initial_condition_gravitational_potential(k)

    If (MG_parametrisation .eq. 'GR') then
       
       dark_energy_density_perturbation_super_sound_horizon = delta_0(k)*(1.d0 + &
            w0_fld)*a*(1.d0/(1.d0 - 3.d0*w0_fld) + 3.d0*H0**2*Omega_m/(k**2*a) ) 

    Else

       dark_energy_density_perturbation_super_sound_horizon = delta_0(k)*(1.d0 + &
            equation_of_state(a))*a*k**2*(1.d0/(1.d0 - 3.d0*equation_of_state(a)) &
            + 3.d0*H0**2*Omega_m/(k**2*a) ) 

    End if

  end function dark_energy_density_perturbation_super_sound_horizon

  ! DARK ENERGY VELOCITY PERTURBATIONS ON SUPER-SOUND HORIZON SCALES. NOTE THIS IS \theta_{DE}; IN DOMENICO'S PAPER THE AUTHORS USED V_{DE} = (1+w)\theta_{DE} 

  function dark_energy_velocity_perturbation_super_sound_horizon(a,k)

    use fiducial

    Implicit none

    Real*8 :: a,k,B,dark_energy_velocity_perturbation_super_sound_horizon

    B = initial_condition_gravitational_potential(k)

    dark_energy_velocity_perturbation_super_sound_horizon = -delta_0(k)*H0*Sqrt(a*Omega_m) 

  end function dark_energy_velocity_perturbation_super_sound_horizon

  ! DARK ENERGY DENSITY PERTURBATIONS ON SUB-SOUND HORIZON SCALES 

  function de_density_perturbation_sub_sound_horizon(k) 

    use fiducial

    Implicit none

    Real*8 :: B,k,de_density_perturbation_sub_sound_horizon

    B = initial_condition_gravitational_potential(k)

    If (MG_parametrisation .eq. 'GR') then

       de_density_perturbation_sub_sound_horizon = -B*(1.d0 + w0_fld)/cs2_fld

    Else

       de_density_perturbation_sub_sound_horizon = -B*(1.d0 + &
            equation_of_state(initial_scale_factor))/sound_speed_squared(initial_scale_factor)

    End if

  end function de_density_perturbation_sub_sound_horizon

  ! DARK ENERGY VELOCITY PERTURBATION ON SUB-SOUND HORIZON SCALES (DOMENICO'S EXPRESSION DIVIDED BY (1+w) ).

  function de_velocity_perturbation_sub_sound_horizon(a,k) 

    use fiducial

    Implicit none 

    Real*8 :: a,B,k,de_velocity_perturbation_sub_sound_horizon

    B = initial_condition_gravitational_potential(k)
    
    If (MG_parametrisation .eq. 'GR') then

       de_velocity_perturbation_sub_sound_horizon = -3.d0*conformal_Hubble_parameter(a)*&
            (cs2_fld - w0_fld)*de_density_perturbation_sub_sound_horizon(k)

    Else

       de_velocity_perturbation_sub_sound_horizon = -3.d0*conformal_Hubble_parameter(a)*&
            (sound_speed_squared(a) - equation_of_state(a))*de_density_perturbation_sub_sound_horizon(k)

    End if

  end function de_velocity_perturbation_sub_sound_horizon

  !###################################################################################
  ! SOME EQUATIONS THAT ARE USED ALONG WITH THE EQUATIONS FOR INITIAL CONDITIONS ABOVE
  !###################################################################################

  ! CONSTANT \delta_0 IN DOMENICO AND MARTIN'S PAPER

  function delta_0(k) 

    use fiducial

    Implicit none

    Real*8 :: B,k,delta_0

    B = initial_condition_gravitational_potential(k)

    delta_0 = -(2.d0*B*k**2)/(3.d0*H0**2*Omega_m) 

  end function delta_0

  ! CONSTANT FOR SUPER-HORIZON SOLUTIONS

  function delta_0_super_horizon(k)

    use fiducial

    Implicit none

    Real*8 :: B,k,delta_0_super_horizon

    B = initial_condition_gravitational_potential(k)

    delta_0_super_horizon = -2.0d0*B*k**2*(1.0d0+w0_fld)*(4.0d0*f_pi-3.d0)/(3.d0*H0**2*Omega_m*&
         (4.d0*e_pi-3.d0*(1.d0+w0_fld)))

  end function delta_0_super_horizon

  ! EFFECTIVE SOUND SPEED 

  function effective_sound_speed_squared(a,k)

    use fiducial

    Implicit none

    Real*8 :: a,k,effective_sound_speed_squared

    If ( MG_parametrisation .eq. 'GR' ) then

       effective_sound_speed_squared = cs2_fld - 2.d0*f_pi/3.d0

    Else

       effective_sound_speed_squared = sound_speed_squared(a) - (2.d0/3.d0)*wavenumber_k**2*FR(a)/&
            (a**2*F_MG(a) - a**2*F_MG(a)**2 + wavenumber_k**2*(2.d0 - 3.d0*F_MG(a))*FR(a))

    End if

  end function effective_sound_speed_squared

  ! DARK ENERGY DENSITY FOR AN EFFECTIVE FLUID FROM F(R) PARAMMETRISATION BY BASILAKOS ET AL.

  function dark_energy_density(a) ! ACTUALLY THIS IS THE DARK ENERGY DENSITY MULTIPLIED BY \kappa

    use fiducial

    Implicit none

    Real*8 :: a,dark_energy_density

    If (MG_parametrisation .eq. 'HS_Basilakos') then

       dark_energy_density = &

!!$-fMG(a)/2.d0 + 3.d0*conformal_Hubble_parameter(a)**2/a**2 - &
!!$            3.d0*conformal_Hubble_parameter(a)**2*F_MG_prime(a)/a + &
!!$            3.d0*F_MG(a)*conformal_Hubble_parameter(a)*derivative_conformal_Hubble_parameter(a)/a

-3*H0**2*(-1 + Omega_m) - (18*a**3*b_fR*H0**2*(-1 + Omega_m)**2*(4*a**6*(-1 &
            + Omega_m)**2 + a**3*(-1 + Omega_m)*Omega_m - 2*Omega_m**2))/(-4*a**3*(-1 + Omega_m) + Omega_m)**3 + &
            (3*a**6*b_fR**2*H0**2*(-1 + Omega_m)**3*(1024*a**18*(-1 + Omega_m)**6 + 9216*a**15*(-1 + &
            Omega_m)**5*Omega_m - 22848*a**12*(-1 + Omega_m)**4*Omega_m**2 + 25408*a**9*(-1 + Omega_m)**3*Omega_m**3 - &
            7452*a**6*(-1 + Omega_m)**2*Omega_m**4 - 4656*a**3*(-1 + Omega_m)*Omega_m**5 + &
            37*Omega_m**6))/(-4*a**3*(-1 + Omega_m) + Omega_m)**8

    Else

       dark_energy_density = 0.d0

    End if

  end function dark_energy_density

  ! DERIVATIVE DARK ENERGY DENSITY FOR AN EFFECTIVE FLUID FROM F(R) PARAMMETRISATION BY BASILAKOS ET AL.

  function derivative_dark_energy_density(a)

    use fiducial

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

  ! EQUATION OF STATE FOR AN EFFECTIVE FLUID FROM F(R) PARAMMETRISATION BY BASILAKOS ET AL.

  function equation_of_state(a)

    use fiducial

    Implicit none

    Real*8 :: a,equation_of_state

    If (MG_parametrisation .eq. 'HS_Basilakos') then

       equation_of_state = -1 - (12*a**3*b_fR*(a**3*(-1 + Omega_m) - Omega_m)*(-1 + Omega_m)*Omega_m*(8*a**3*(-1 + &
            Omega_m) + Omega_m))/(-4*a**3*(-1 + Omega_m) + Omega_m)**4 - (2*a**6*b_fR**2*(-1 + &
            Omega_m)**2*Omega_m*(40960*a**18*(-1 + Omega_m)**6 - 79872*a**15*(-1 + Omega_m)**5*Omega_m + &
            75264*a**12*(-1 + Omega_m)**4*Omega_m**2 + 18448*a**9*(-1 + Omega_m)**3*Omega_m**3 - 65028*a**6*(-1 + &
            Omega_m)**2*Omega_m**4 - 6648*a**3*(-1 + Omega_m)*Omega_m**5 + 109*Omega_m**6))/(4*a**3*(-1 + &
            Omega_m) - Omega_m)**9

    Else

       equation_of_state = w0_fld

    End if

   end function equation_of_state

   function derivative_equation_of_state(a)

     use fiducial

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

  ! ADIABATIC SOUND SPEED FOR AN EFFECTIVE FLUID FROM F(R) PARAMETRISATION BY BASILAKOS ET AL. 

  function adiabatic_sound_speed_squared(a)

    use fiducial

    Implicit none

    Real*8 :: a, adiabatic_sound_speed_squared
    Real*8,parameter :: epsilon = 1.d-3  ! VALUE TAKEN FROM 'PHENOMENOLOGICAL APPROACH TO ...' BY M. KUNZ 

    adiabatic_sound_speed_squared = -1 - (12*a**3*b_fR*(a**3*(-1 + Omega_m) - Omega_m)*(-1 + &
         Omega_m)*Omega_m*(8*a**3*(-1 + Omega_m) + Omega_m))/(-4*a**3*(-1 + Omega_m) + Omega_m)**4 - &
         (2*a**6*b_fR**2*(-1 + Omega_m)**2*Omega_m*(40960*a**18*(-1 + Omega_m)**6 - 79872*a**15*(-1 + &
         Omega_m)**5*Omega_m + 75264*a**12*(-1 + Omega_m)**4*Omega_m**2 + 18448*a**9*(-1 + &
         Omega_m)**3*Omega_m**3 - 65028*a**6*(-1 + Omega_m)**2*Omega_m**4 - 6648*a**3*(-1 + &
         Omega_m)*Omega_m**5 + 109*Omega_m**6))/(4*a**3*(-1 + Omega_m) - Omega_m)**9 - &
         (8*a**6*b_fR**2*(-1 + Omega_m)**2*Omega_m**2*((-6*(a**3*(-1 + Omega_m) - Omega_m)*(8*a**3*(-1 + &
         Omega_m) + Omega_m))/(-4*a**3*(-1 + Omega_m) + Omega_m)**4 - (a**3*b_fR*(-1 + Omega_m)*(40960*a**18*(-1&
         + Omega_m)**6 - 79872*a**15*(-1 + Omega_m)**5*Omega_m + 75264*a**12*(-1 + Omega_m)**4*Omega_m**2 + &
         18448*a**9*(-1 + Omega_m)**3*Omega_m**3 - 65028*a**6*(-1 + Omega_m)**2*Omega_m**4 - 6648*a**3*(-1 + &
         Omega_m)*Omega_m**5 + 109*Omega_m**6))/(4*a**3*(-1 + Omega_m) - Omega_m)**9)*((-24*a**3*(a**3*(-1 + &
         Omega_m) - Omega_m)*(-1 + Omega_m))/(-4*a**3*(-1 + Omega_m) + Omega_m)**4 - (48*a**3*(a**3*(-1 + &
         Omega_m) - Omega_m)*(-1 + Omega_m)*(8*a**3*(-1 + Omega_m) + Omega_m))/(-4*a**3*(-1 + Omega_m) + Omega_m)**5 - &
         (3*(a**3*(-1 + Omega_m) - Omega_m)*(8*a**3*(-1 + Omega_m) + Omega_m))/(-4*a**3*(-1 + Omega_m) + Omega_m)**4 -& 
         (3*a**3*(-1 + Omega_m)*(8*a**3*(-1 + Omega_m) + Omega_m))/(-4*a**3*(-1 + Omega_m) + Omega_m)**4 - &
         (12*a**6*b_fR*(-1 + Omega_m)**2*(10240*a**15*(-1 + Omega_m)**5 - 16640*a**12*(-1 + Omega_m)**4*Omega_m + &
         12544*a**9*(-1 + Omega_m)**3*Omega_m**2 + 2306*a**6*(-1 + Omega_m)**2*Omega_m**3 - 5419*a**3*(-1 + &
         Omega_m)*Omega_m**4 - 277*Omega_m**5))/(4*a**3*(-1 + Omega_m) - Omega_m)**9 - (a**3*b_fR*(-1 + &
         Omega_m)*(40960*a**18*(-1 + Omega_m)**6 - 79872*a**15*(-1 + Omega_m)**5*Omega_m + 75264*a**12*(-1 + &
         Omega_m)**4*Omega_m**2 + 18448*a**9*(-1 + Omega_m)**3*Omega_m**3 - 65028*a**6*(-1 + Omega_m)**2*Omega_m**4 - &
         6648*a**3*(-1 + Omega_m)*Omega_m**5 + 109*Omega_m**6))/(4*a**3*(-1 + Omega_m) - Omega_m)**9 + (18*a**6*b_fR*(-1&
         + Omega_m)**2*(40960*a**18*(-1 + Omega_m)**6 - 79872*a**15*(-1 + Omega_m)**5*Omega_m + 75264*a**12*(-1 +&
         Omega_m)**4*Omega_m**2 + 18448*a**9*(-1 + Omega_m)**3*Omega_m**3 - 65028*a**6*(-1 + Omega_m)**2*Omega_m**4 - &
         6648*a**3*(-1 + Omega_m)*Omega_m**5 + 109*Omega_m**6))/(-4*a**3*(-1 + Omega_m) + Omega_m)**10))/&
         (epsilon + ((12*a**3*b_fR*(a**3*(-1 + Omega_m) - Omega_m)*(-1 + Omega_m)*Omega_m*(8*a**3*(-1 + Omega_m) + &
         Omega_m))/(-4*a**3*(-1 + Omega_m) + Omega_m)**4 + (2*a**6*b_fR**2*(-1 + Omega_m)**2*Omega_m*(40960*a**18*(-1 + &
         Omega_m)**6 - 79872*a**15*(-1 + Omega_m)**5*Omega_m + 75264*a**12*(-1 + Omega_m)**4*Omega_m**2 + &
         18448*a**9*(-1 + Omega_m)**3*Omega_m**3 - 65028*a**6*(-1 + Omega_m)**2*Omega_m**4 - 6648*a**3*(-1 + &
         Omega_m)*Omega_m**5 + 109*Omega_m**6))/(4*a**3*(-1 + Omega_m) - Omega_m)**9)**2)

  end function adiabatic_sound_speed_squared

  ! ANISOTROPIC STRESS FOR AN EFFECTIVE FLUID FROM F(R) PARAMETRISATION BY BASILAKOS ET AL.

  function anisotropic_stress(a,y1,y2,y3,y4)

    use fiducial

    Implicit none

    Real*8 :: a,anisotropic_stress,y1,y2,y3,y4

    If (MG_parametrisation .eq. 'HS_Basilakos') then

       anisotropic_stress = wavenumber_k**2*FR(a)/(a**2*F_MG(a)**2 + &
            3.d0*wavenumber_k**2*FR(a)*F_MG(a) )

!!$wavenumber_k**2*FR(a)/(a**2*F_MG(a) - &
!!$         a**2*F_MG(a)**2 + wavenumber_k**2*(2.d0 - 3.d0*F_MG(a))*FR(a))

    Else if (MG_parametrisation .eq. 'Starobinsky_Basilakos') then

       anisotropic_stress = -2.d10

    Else if (MG_parametrisation .eq. 'GR') then

       anisotropic_stress = e_pi*Delta_matter(a,wavenumber_k,y1,y3) + (f_pi/(1.d0 + &
            g_pi**2*conformal_Hubble_parameter(a)**2/wavenumber_k**2))*Delta_dark_energy(a,wavenumber_k,y2,y4)

    Else

       anisotropic_stress = -3.d10

    End if

  end function anisotropic_stress

  function anisotropic_stress_fr(a, params) bind(c)
    
    real(c_double), value :: a
    type(c_ptr), value :: params
    real(c_double) :: anisotropic_stress_fr

    anisotropic_stress_fr = anisotropic_stress(a,0.d0,0.d0,0.d0,0.d0)

  end function anisotropic_stress_fr

  function derivative_anisotropic_stress(a)

    use fiducial

    real(fgsl_double) :: result, abserr
    integer(fgsl_int) :: status
    type(fgsl_function) :: pwr
    real(fgsl_double) :: derivative_anisotropic_stress,a

    If (MG_parametrisation .eq. 'HS_Basilakos') then

       pwr = fgsl_function_init(anisotropic_stress_fr, c_null_ptr)

       status = fgsl_deriv_central (pwr, a , 1.E-8_fgsl_double, &
         result, abserr)

       derivative_anisotropic_stress = result

    Else if (MG_parametrisation .eq. 'Starobinsky_Basilakos') then

       derivative_anisotropic_stress = -2.d10

    Else if (MG_parametrisation .eq. 'GR') then

       derivative_anisotropic_stress = -3.d10 

    Else

       derivative_anisotropic_stress = -4.d10

    End if

  end function derivative_anisotropic_stress

  ! SOUND SPEED SQUARED FOR AN EFFECTIVE FLUID FROM F(R) PARAMETRISATION BY BASILAKOS ET AL. 

  function sound_speed_squared(a)

    use fiducial

    Implicit none

    Real*8 :: a,sound_speed_squared

    If ( a .lt. switch_off_pressure_perturbation_terms) then

       sound_speed_squared = (2.d0*wavenumber_k**2*FR(a))/(3.d0*a**2*F_MG(a) - &
            3.d0*a**2*F_MG(a)**2 + 3.d0*wavenumber_k**2*(2.d0 - 3.d0*F_MG(a))*FR(a)) &

            +(conformal_Hubble_parameter(a)**2*FR_prime(a)/a + conformal_Hubble_parameter(a)*&
            derivative_conformal_Hubble_parameter(a)*FR_prime(a)&
            + conformal_Hubble_parameter(a)**2*FR_double_prime(a))/&
            (1.d0 - F_MG(a) + wavenumber_k**2*(2.d0 - 3.d0*F_MG(a))*FR(a)/a**2/F_MG(a) ) &

            + ( a*conformal_Hubble_parameter(a)**2*& 
            F_MG_prime(a) + a**2*conformal_Hubble_parameter(a)*&
            derivative_conformal_Hubble_parameter(a)*F_MG_prime(a) + a**2*conformal_Hubble_parameter(a)**2*&
            F_MG_double_prime(a))/( wavenumber_k**2 - F_MG(a)*wavenumber_k**2 + wavenumber_k**4*(2.d0 - &
            3.d0*F_MG(a))*FR(a)/a**2/F_MG(a)  )

    Else

       sound_speed_squared = (2.d0*wavenumber_k**2*FR(a))/(3.d0*a**2*F_MG(a) - &
            3.d0*a**2*F_MG(a)**2 + 3.d0*wavenumber_k**2*(2.d0 - 3.d0*F_MG(a))*FR(a)) 

    End if

  end function sound_speed_squared

  ! DARK ENERGY PRESSURE PERTURBATION OVER DARK ENERGY DENSITY AND OTHER FUNCTIONS FOR AN EFFECTIVE FLUID FROM F(R) PARAMETRISATION BY BASILAKOS ET AL. 

  function fMG(a)

    use fiducial

    Implicit none

    Real*8 :: fMG,a

    fMG = 6*conformal_Hubble_parameter(a)*(conformal_Hubble_parameter(a)/a**2 + derivative_conformal_Hubble_parameter(a)/a - &
         (2*H0**2*(-1 + Omega_m)*(conformal_Hubble_parameter(a) + &
         a*derivative_conformal_Hubble_parameter(a)))/(a**2*b_fR*H0**2*(-1 + Omega_m) - &
         2*conformal_Hubble_parameter(a)**2 - 2*a*conformal_Hubble_parameter(a)*derivative_conformal_Hubble_parameter(a)))

  end function fMG

  function F_MG(a)

    use fiducial

    Implicit none 

    Real*8 :: F_MG,a

    F_MG = 1 - (2*a**6*b_fR*(-1 + Omega_m)**2*(-4*a**3*(-1 + Omega_m) + &
         Omega_m)**18)/(-16384*a**30*(-64 + 8*b_fR + b_fR**2)*(-1 + Omega_m)**10 - 8192*a**27*(320 - &
         24*b_fR + b_fR**2)*(-1 + Omega_m)**9*Omega_m + 1536*a**24*(1920 - 40*b_fR + 31*b_fR**2)*(-1 + &
         Omega_m)**8*Omega_m**2 + 128*a**21*(-15360 - 1080*b_fR + 43*b_fR**2)*(-1 + Omega_m)**7*Omega_m**3 + &
         64*a**18*(13440 + 2496*b_fR + 3817*b_fR**2)*(-1 + Omega_m)**6*Omega_m**4 - 384*a**15*(672 + &
         193*b_fR + 844*b_fR**2)*(-1 + Omega_m)**5*Omega_m**5 - 8*a**12*(-6720 - 2226*b_fR + 7307*b_fR**2)*(-1 + &
         Omega_m)**4*Omega_m**6 + 10*a**9*(-768 - 222*b_fR + 37*b_fR**2)*(-1 + Omega_m)**3*Omega_m**7 + &
         120 *a**6*(6 + b_fR)*(-1 + Omega_m)**2*Omega_m**8 - a**3*(40 + b_fR)*(-1 + Omega_m)*Omega_m**9 + Omega_m**10)**2
 
  end function F_MG

  function F_MG_prime(a)

    use fiducial

    Implicit none 

    Real*8 :: F_MG_prime,a

    F_MG_prime = (-8*a**3*b_fR*H0**4*(-1 + Omega_m)**2*(-2*conformal_Hubble_parameter(a)**2 + &
         a**2*derivative_conformal_Hubble_parameter(a)**2 + &
         a*conformal_Hubble_parameter(a)*(a*second_derivative_conformal_Hubble_parameter(a) + &
         derivative_conformal_Hubble_parameter(a))))/(a**2*b_fR*H0**2*(-1 + Omega_m) - &
         2*conformal_Hubble_parameter(a)**2 - &
         2*a*conformal_Hubble_parameter(a)*derivative_conformal_Hubble_parameter(a))**3 
 
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

    use fiducial 

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

    use fiducial 

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

  function dark_energy_pressure_perturbation(a)

    use fiducial

    Implicit none 

    Real*8 :: dark_energy_pressure_perturbation,a

    If ( a .lt. switch_off_pressure_perturbation_terms) then

       dark_energy_pressure_perturbation = 2.d0*wavenumber_k**2*FR(a)/(3.d0*a**2*F_MG(a)**2 + &
            9.d0*wavenumber_k**2*FR(a)*F_MG(a) ) + &

            ((a*conformal_Hubble_parameter(a)**2 + a**2*conformal_Hubble_parameter(a)*&
            derivative_conformal_Hubble_parameter(a))*FR_prime(a) + &
            a**2*conformal_Hubble_parameter(a)**2*FR_double_prime(a) )/(a**2*F_MG(a)  &
            + 3.d0*wavenumber_k**2*FR(a)) + &

            ( a*conformal_Hubble_parameter(a)**2*F_MG_prime(a) + a**2*conformal_Hubble_parameter(a)*&
            derivative_conformal_Hubble_parameter(a)*F_MG_prime(a) + a**2*conformal_Hubble_parameter(a)**2*&
            F_MG_double_prime(a))/(F_MG(a)*wavenumber_k**2 + 3.d0*wavenumber_k**4*FR(a)/a**2  )

    Else

       dark_energy_pressure_perturbation = 2.d0*wavenumber_k**2*FR(a)/(3.d0*a**2*F_MG(a)**2 + &
            9.d0*wavenumber_k**2*FR(a)*F_MG(a) ) 

    End if

  end function dark_energy_pressure_perturbation

  function second_dark_energy_pressure_perturbation(a,params) bind(c)

    !use fiducial

    real(c_double),value :: a
    real(c_double) :: second_dark_energy_pressure_perturbation
    type(c_ptr),value :: params

    second_dark_energy_pressure_perturbation = sound_speed_squared(a) 

  end function second_dark_energy_pressure_perturbation

  function derivative_second_dark_energy_pressure_perturbation(a)

    real(fgsl_double) :: result, abserr
    integer(fgsl_int) :: status
    type(fgsl_function) :: pwr
    real(fgsl_double) :: derivative_second_dark_energy_pressure_perturbation,a

    pwr = fgsl_function_init(second_dark_energy_pressure_perturbation, c_null_ptr)

    status = fgsl_deriv_central (pwr, a , 1.E-8_fgsl_double, &
         result, abserr)

    derivative_second_dark_energy_pressure_perturbation = result

  end function derivative_second_dark_energy_pressure_perturbation

  ! CONFORMAL HUBBLE PARAMETER

  function conformal_Hubble_parameter(a)

    use fiducial

    Implicit none

    Real*8 :: a,conformal_Hubble_parameter
    Real*8,parameter :: om = 3.d-1

    If (MG_parametrisation .eq. 'HS_Basilakos') then

       conformal_Hubble_parameter = a*H0*Sqrt(1.d0 - Omega_m + Omega_m/a**3 + (2.d0*a**2*b_fR*(-1.d0 + &
            Omega_m)**2*(12.d0*a**7*(-1.d0 + Omega_m)**2 + 3.d0*a**4*(-1.d0 + Omega_m)*Omega_m - &
            6.d0*a*Omega_m**2))/(4.d0*a**3*(-1.d0 + Omega_m) - Omega_m)**3 - (a**5*b_fR**2*(-1.d0 + &
            Omega_m)**3*(-1024.d0*a**19*(-1.d0 + Omega_m)**6 - 9216.d0*a**16*(-1.d0 + Omega_m)**5*Omega_m + &
            22848.d0*a**13*(-1.d0 + Omega_m)**4*Omega_m**2 - 25408.d0*a**10*(-1.d0 + Omega_m)**3*Omega_m**3 + &
            7452.d0*a**7*(-1.d0 + Omega_m)**2*Omega_m**4 + 4656.d0*a**4*(-1.d0 + Omega_m)*Omega_m**5 - &
            37.d0*a*Omega_m**6))/(-4.d0*a**3*(-1.d0 + Omega_m) + Omega_m)**8)

    Else if (MG_parametrisation .eq. 'Starobinsky_Basilakos') then

       conformal_Hubble_parameter = -2.d0

    Else if (MG_parametrisation .eq. 'GR') then

       conformal_Hubble_parameter = H0*Sqrt(Omega_m/a + (1.d0 - Omega_m)*a**(-1.d0 - 3.d0*w0_fld))

    Else

       conformal_Hubble_parameter = -1.d10

    End if

  end function conformal_Hubble_parameter

  ! DERIVATIVE OF THE CONFORMAL HUBBLE PARAMETER 

  function derivative_conformal_Hubble_parameter(a) 

    use fiducial

    Implicit none 

    Real*8 :: a,derivative_conformal_Hubble_parameter

    If (MG_parametrisation .eq. 'HS_Basilakos') then

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

    Else if (MG_parametrisation .eq. 'GR') then

       derivative_conformal_Hubble_parameter = -(H0**2*(Omega_m + (1.d0 - Omega_m)*(1.d0 + &
            3.d0*w0_fld)/a**(3.d0*w0_fld)))/(2.d0*a**2*conformal_Hubble_parameter(a))

    Else

       derivative_conformal_Hubble_parameter = -1.d10

    End if

  end function derivative_conformal_Hubble_parameter

  function second_derivative_conformal_Hubble_parameter(a) 

    use fiducial

    Implicit none 

    Real*8 :: a,second_derivative_conformal_Hubble_parameter

    If (MG_parametrisation .eq. 'HS_Basilakos') then

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

    Else if (MG_parametrisation .eq. 'GR') then

       second_derivative_conformal_Hubble_parameter = -3.d10
    Else

       second_derivative_conformal_Hubble_parameter = -1.d10

    End if

  end function second_derivative_conformal_Hubble_parameter

  ! SCALE FACTOR AT WHICH THE MODE CROSSES THE HUBBLE HORIZON IN A FLAT MATTER DOMINATED UNIVERSE

  function scale_factor_horizon_crossing(k)

    use fiducial

    Implicit none

    Real*8 :: k,scale_factor_horizon_crossing

    scale_factor_horizon_crossing = H0**2*Omega_m/(speedL**2*k**2)

  end function scale_factor_horizon_crossing

  ! SCALE FACTOR AT WHICH THE MODE CROSSES THE EFFECTIVE SOUND HORIZON IN A FLAT MATTER DOMINATED UNIVERSE

  function scale_factor_effective_sound_horizon(a,k) 

    use fiducial

    Implicit none

    Real*8 :: a,k,scale_factor_effective_sound_horizon

    scale_factor_effective_sound_horizon = H0**2*Omega_m/(effective_sound_speed_squared(a,k)*k**2)

  end function scale_factor_effective_sound_horizon

  ! GAUGE-INVARIANT COMOVING DENSITY CONTRAST (MATTER) 

  function Delta_matter(a,k,y1,y3) 

    use fiducial

    Implicit none

    Real*8 :: a,k,y1,y3,Delta_matter

    If (MG_parametrisation .eq. 'GR') then

       Delta_matter = y1 + 3.d0*conformal_Hubble_parameter(a)*y3/k**2

    Else

       Delta_matter = y1 + 3.d0*conformal_Hubble_parameter(a)*y3/k**2

    End if

  end function Delta_matter

  ! GAUGE-INVARIANT COMOVING DENSITY CONTRAST (DARK ENERGY) 

  function Delta_dark_energy(a,k,y2,y4) 

    use fiducial

    Implicit none

    Real*8 :: a,k,y2,y4,Delta_dark_energy

    If (MG_parametrisation .eq. 'GR') then

       Delta_dark_energy = y2 + 3.d0*conformal_Hubble_parameter(a)*(1.d0 + w0_fld)*y4/k**2

    Else

       Delta_dark_energy = y2 + 3.d0*conformal_Hubble_parameter(a)*y4/k**2 ! CAREFUL WITH DIFFERENCES IN y4

    End if

  end function Delta_dark_energy

  ! DARK ENERGY ANISOTROPIC STRESS 

  function sigma(a,y1,y2,y3,y4)

    use fiducial

    Implicit none

    Real*8 :: a,y1,y2,y3,y4,sigma

    If (MG_parametrisation .eq. 'GR') then

       sigma = 2.d0*anisotropic_stress(a,y1,y2,y3,y4)/(3.d0*(1.d0+w0_fld))

    Else

       sigma = 2.d0*anisotropic_stress(a,y1,y2,y3,y4)/(3.d0*(1.d0 + &
            equation_of_state(a) ))

    End if

  end function sigma

  function Geff_over_GN(a)

    use fiducial 

    Implicit none

    Real*8 :: Geff_over_GN,a

    Geff_over_GN = (1.d0 + 4.d0*wavenumber_k**2*FR(a)/a**2/F_MG(a))/&
         (F_MG(a) + 3.d0*wavenumber_k**2*FR(a)/a**2)

  end function Geff_over_GN

  function Qeff(a)

    use fiducial

    Implicit none

    Real*8 :: Qeff,a

    Qeff = (1.d0 + 2.d0*wavenumber_k**2*FR(a)/a**2/F_MG(a))/&
         (F_MG(a) + 3.d0*wavenumber_k**2*FR(a)/a**2)

  end function Qeff

  !####################################################################################################################################################
  ! EXPRESSIONS FOR THE POTENTIALS AND THEIR INITIAL CONDITIONS (SET IN MATTER DOMINANCE) FOLLOW SOLUTIONS IN THE PAPER BY GUILLERMO, MARTIN, AND LUKAS
  !####################################################################################################################################################

  ! INITIAL CONDITION FOR THE GRAVITATIONAL POTENTIAL IN MATTER DOMINATED ERA. DIMENSIONLESS

  function initial_condition_gravitational_potential(k) 

    use fiducial

    Implicit none

    Real*8 :: k,initial_condition_gravitational_potential

    If ( wavenumber_k/ks .lt. lower_limit_ks ) then

       initial_condition_gravitational_potential = -3.d0*Sqrt(2.d0*Pi**2*H0**3*&
            primordial_dimensionless_power_spectrum(k)/k**3)/5.d0

    Else if ( wavenumber_k/ks .gt. lower_limit_ks ) then

       initial_condition_gravitational_potential = -3.d0*Sqrt(2.d0*Pi**2*H0**3*&
            primordial_dimensionless_power_spectrum(k)/k**3)/5.d0*3.d0*&
            (ks/k)**2*log(k/ks)
    
    Else

       initial_condition_gravitational_potential = 1.d10

    End if

!    initial_condition_gravitational_potential = -3.d0/2.d0*delta_0(k)*H0**2*Omega_m/k**2

  end function initial_condition_gravitational_potential

  function primordial_dimensionless_power_spectrum(k)

    use fiducial

    Implicit none

    Real*8 :: k, primordial_dimensionless_power_spectrum

    primordial_dimensionless_power_spectrum = A_s*(k/kp)**(n_s- 1.d0)

  end function primordial_dimensionless_power_spectrum

  ! POTENTIAL \Phi

  function phi(a,k,y1,y2,y3,y4)

    use fiducial

    Implicit none

    Real*8 :: a,k,y1,y2,y3,y4,phi

    If (MG_parametrisation .eq. 'GR') then

       phi = -((3.d0*H0**2)/(2.d0*k**2))*(Omega_m*Delta_matter(a,k,y1,y3)/a + &
            (1 - Omega_m)*Delta_dark_energy(a,k,y2,y4)/a**(1.d0 + 3.d0*w0_fld))

    Else

       phi = -(3.d0*H0**2/2.d0/k**2)*(Omega_m*Delta_matter(a,k,y1,y3)/a + &
            Omega_DE(a)*a**2*Delta_dark_energy(a,k,y2,y4) )

    End if

  end function phi

  ! POTENTIAL \Psi 

  function psi(a,k,y1,y2,y3,y4) 

    use fiducial

    Implicit none

    Real*8 :: a,k,y1,y2,y3,y4,psi

    If (MG_parametrisation .eq. 'GR') then

       psi = -9.d0*H0**2*(1.d0+w0_fld)*sigma(a,y1,y2,y3,y4)*(1.d0-&
            Omega_m)/(2.d0*k**2*a**(1.d0+3.d0*w0_fld)) + phi(a,k,y1,y2,y3,y4)

    Else

       psi = phi(a,k,y1,y2,y3,y4) - 3.d0*H0**2/k**2*anisotropic_stress(a,y1,y2,y3,y4)*&
            Omega_DE(a)*a**2

    End if

  end function psi

  function Omega_DE(a)

    use fiducial

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

  !################################################################################################################
  ! SOLUTIONS FOR DARK ENERGY PERTURBATIONS IN MATTER DOMINATED REGIME WHEN CONSIDERING THE DARK ENERGY ANISOTROPIC
  ! STRESS MODEL DEFINED ABOVE
  !################################################################################################################

  ! DARK ENERGY DENSITY PERTURBATIONS ON SUPER-HORIZON SCALES 

  function dark_energy_density_perturbation_super_horizon(k)

    use fiducial

    Implicit none

    Real*8 :: k,dark_energy_density_perturbation_super_horizon

    dark_energy_density_perturbation_super_horizon = delta_0_super_horizon(k)*(((4.d0*e_pi - &
         3.d0*(1.d0 + w0_fld))*( 3.d0*H0**2*Omega_m))/(k**2*(-3.d0 + 4.d0*f_pi) ) )

  end function dark_energy_density_perturbation_super_horizon

  ! DARK ENERGY VELOCITY PERTURBATIONS ON SUPER-HORIZON SCALES. NOTE THE VARIABLE USED HERE IS V_{de} = -(1+w)*\theta/conformal_H 

  function dark_energy_velocity_perturbation_super_horizon(a,k)

    use fiducial

    Implicit none 

    Real*8 :: a,k,dark_energy_velocity_perturbation_super_horizon

    dark_energy_velocity_perturbation_super_horizon = a*delta_0_super_horizon(k)*(4.d0*e_pi-&
         3.d0*(1.d0+w0_fld))/(4.d0*f_pi - 3.d0)

  end function dark_energy_velocity_perturbation_super_horizon

  ! DARK ENERGY DENSITY PERTURBATIONS ON SUB-SOUND HORIZON SCALES

  function dark_energy_density_perturbation_sub_sound_horizon(a,k)

    use fiducial

    Implicit none

    Real*8 :: a,k,dark_energy_density_perturbation_sub_sound_horizon

    dark_energy_density_perturbation_sub_sound_horizon = a*delta_0(k)*(2.d0*e_pi/(3.d0*&
         effective_sound_speed_squared(a,k)) + 3.d0*H0**2*Omega_m*(1.d0 + &
         equation_of_state(a))/(2.d0*a*effective_sound_speed_squared(a,k)*k**2))

  end function dark_energy_density_perturbation_sub_sound_horizon

  ! DARK ENERGY VELOCITY PERTURBATIONS ON SUB-SOUND HORIZON SCALES

  function dark_energy_velocity_perturbation_sub_sound_horizon(a,k)

    use fiducial

    Implicit none

    Real*8 :: a,k,dark_energy_velocity_perturbation_sub_sound_horizon

    dark_energy_velocity_perturbation_sub_sound_horizon = delta_0(k)*(2.d0*a*e_pi*(1.d0 + &
         3.d0*effective_sound_speed_squared(a,k) + 2.d0*f_pi - &
         3.d0*equation_of_state(a))/(3.d0*effective_sound_speed_squared(a,k)) + &
         3.d0*H0**2*Omega_m*(3.d0*effective_sound_speed_squared(a,k) + &
         2.d0*f_pi - 3.d0*equation_of_state(a))*(1.d0 + equation_of_state(a))/&
         (2.d0*k**2*effective_sound_speed_squared(a,k)))

  end function dark_energy_velocity_perturbation_sub_sound_horizon

  !##########################################################################################################################
  ! SOLUTIONS FOR BOTH DARK MATTER AND DARK ENERGY PERTURBATIONS ON SUPER-HORIZON SCALES IN DARK ENERGY DOMINANCE. THE DARK 
  ! ENERGY ANISOTROPIC STRESS MODEL IS NOT EXTERNALLY SOURCED AND w=-1. THE CONSTANT \delta_0' IN THE DRAFT IS FIXED TO MATCH 
  ! THE SOLUTIONS FOR MATTER DENSITY PERTURBATIONS; IT IS FIXED ET REDSHIFT z~ 0.7 
  !##########################################################################################################################

  ! THE CONSTANT IN THE DRAFT

  function d01(a,k)

    use fiducial

    Implicit none

    Real*8 :: a,B,k,d01

    B = initial_condition_gravitational_potential(k)

    d01 = -(2.d0*a**(1.d0-2.d0*f_pi)*(3.d0+7.d0*f_pi+4.d0*f_pi**2)*k**2*delta_0(k))/(9.d0*f_pi*H0**2*(1.d0-Omega_m))

  end function d01

  ! MATTER DENSITY PERTURBATION

  function dm3(a,k) 

    use fiducial

    Implicit none 

    Real*8 ::  a,k,dm3

    dm3 = (d01(0.6d0,k)*a**(2.d0*f_pi)/2.d0)*( -(18.d0*a**2*f_pi*H0**4*(1.d0&
         -Omega_m)**2)/((2.d0+8.d0*f_pi**2/3.d0 + 14.d0*f_pi/3.d0)*k**4) - 3.d0*(6.d0*(2.d0+8.d0*f_pi**2 + 10.d0*f_pi)*f_pi + &
         (13.d0+16.d0*f_pi**2+32.d0*f_pi)*g_pi**2)*H0**2*(1.d0 - Omega_m)/((2.d0+192.d0*f_pi**2/9.d0 + 32.d0*f_pi**3/3.d0 + &
         38.d0*f_pi/3.d0)*g_pi**2*k**2))         

  end function dm3

  ! MATTER VELOCITY PERTURBATION

  function Vm3(a,k)

    use fiducial

    Implicit none

    Real*8 :: a,k,Vm3

    Vm3 = -d01(0.6d0,k)*a**(2.d0*f_pi)*3.d0*f_pi*H0**2*(1.d0-Omega_m)/(k**2*(2.d0+8.d0*f_pi**2/3.d0 + 14.d0*f_pi/3.d0))

  end function Vm3

  ! DARK ENERGY DENSITY PERTURBATION 

  function dd3(a,k)

    use fiducial

    Implicit none 

    Real*8 :: a,k,dd3

    dd3 = d01(0.6d0,k)*a**(2.d0*f_pi)*(3.d0+2.d0*f_pi)*H0**2*(1.d0-Omega_m)/(k**2*(1.d0+4.d0*f_pi/3.d0))

  end function dd3

  ! DARK ENERGY VELOCITY PERTURBATION 

  function Vd3(a,k)

    use fiducial

    Implicit none 

    Real*8 :: a,k,Vd3

    Vd3 = d01(0.6d0,k)*a**(-2.d0+2.d0*f_pi)

  end function Vd3

  !###############################################################################################################
  ! SOLUTIONS FOR BOTH DARK ENERGY DENSITY AND BOTH DARK ENERGY VELOCITY PERTURBATIONS IN THE f(R) PARAMETRISATION
  ! BY BASILAKOS ET AL. IN MATTER DOMINATED REGIME
  !###############################################################################################################

  function dark_energy_density_perturbation(a,y1)

    use fiducial

    Implicit none

    Real*8 :: dark_energy_density_perturbation,a,y1

    dark_energy_density_perturbation = &

!!$-delta_0(wavenumber_k)*b_fR*wavenumber_k**2*&
!!$         (1.d0-Omega_m)*a**5/3.d0/Omega_m**2/H0**2 - 52.d0*a**4*b_fR*delta_0(wavenumber_k)*&
!!$         (1.d0-Omega_m)/35.d0/Omega_m - 81.d0*a**2*b_fR*delta_0(wavenumber_k)*(1.d0-Omega_m)*&
!!$         (2.d0*Omega_m - 9.d0*a**3*b_fR*(1.d0-Omega_m))*H0**4/5.d0/wavenumber_k**4

         ((1.d0 - F_MG(a) + &
         wavenumber_k**2*(2.d0 - 3.d0*F_MG(a))*FR(a)/(a**2*F_MG(a)))/(1.d0 + &
         3.d0*wavenumber_k**2*FR(a)/(a**2*F_MG(a))))*Omega_m*y1/F_MG(a)/(a**3*Omega_DE(a))

  end function dark_energy_density_perturbation

  function dark_energy_velocity_perturbation(a,y1)

    use fiducial

    Implicit none

    Real*8 :: dark_energy_velocity_perturbation,a,y1

    dark_energy_velocity_perturbation = &

!!$12.d0*delta_0(wavenumber_k)*(1.d0-Omega_m)*b_fR*a**(7.d0/2.d0)*H0/&
!!$         5.d0/Sqrt(Omega_m) + H0**3*108.d0*a**(5.d0/2.d0)*b_fR*delta_0(wavenumber_k)*(1.d0-Omega_m)*Sqrt(Omega_m)/&
!!$         13.d0/wavenumber_k**2 

         (Omega_m*y1)/(a**3*Omega_DE(a))*(conformal_Hubble_parameter(a)*wavenumber_k**2*&
         FR_prime(a)/F_MG(a)/a + a*conformal_Hubble_parameter(a)*F_MG_prime(a)/F_MG(a)&
         )/(1.d0 + 3.d0*wavenumber_k**2*FR(a)/a**2/F_MG(a))

  end function dark_energy_velocity_perturbation
  

  function derivative_matter_perturbation(X,y1,y2,y3,y4,y5,y6)

    use fiducial

    Implicit none

    Real*8 :: derivative_matter_perturbation,X,y1,y2,y3,y4,y5,y6

    If (MG_parametrisation .eq. 'GR') then

       derivative_matter_perturbation = (X**(-3.d0 - 6.d0*w0_fld)*(9.d0*H0**3*wavenumber_k**4*&
            X**(1.d0 + 3.d0*w0_fld)*(-2.d0*(-1.d0 + &
          Omega_m)*e_pi + Omega_m*X**(3.d0*w0_fld))*(1.d0 + Omega_m*(-1.d0&
          + X**(3.d0*w0_fld)))*y1 - 9.d0*H0**3*(-1.d0 + Omega_m)*wavenumber_k**2*(1.d0 + &
          Omega_m*(-1.d0 + &
          X**(3.d0*w0_fld)))*(wavenumber_k**2*(1.d0 + 2.d0*f_pi)*X**(1.d0 + 3.d0*w0_fld)&
          + 2.d0*H0**2*g_pi*(1.d0 + Omega_m*(-1.d0 + X**(3.d0*w0_fld))))*y2 + Sqrt(X**(-1.d0 - &
          3.d0*w0_fld)*(1.d0 + Omega_m*(-1.d0 + X**(3.d0*w0_fld))))*(wavenumber_k**2*X**(1.d0&
          + 3.d0*w0_fld)*(9.d0*H0**2*Omega_m*wavenumber_k**2*X**(1.d0 + 6.d0*w0_fld) - &
          2.d0*wavenumber_k**4*X**(2.d0 + &
          6.d0*w0_fld) + 27.d0*H0**4*(2.d0*e_pi - 2.d0*Omega_m*e_pi +&
          Omega_m*X**(3.d0*w0_fld))*(1.d0 + Omega_m*(-1.d0 + X**(3.d0*w0_fld))))*y3 - &
          9.d0*H0**2*(1.d0 + &
          w0_fld)*(-1.d0 + Omega_m)*(wavenumber_k**4*X**(2.d0 + 6.d0*w0_fld) +&
          3.d0*H0**2*wavenumber_k**2*(1.d0 + 2.d0*f_pi)*X**(1.d0 + 3.d0*w0_fld)*(1.d0 + &
          Omega_m*(-1.d0 + X**(3.d0*w0_fld))) + &
          6.d0*H0**4*g_pi*(1.d0 + Omega_m*(-1.d0 +&
          X**(3.d0*w0_fld)))**2)*y4)))/(2.d0*H0*wavenumber_k**6*(1.d0 + &
          Omega_m*(-1.d0 + X**(3.d0*w0_fld))))

    Else if (MG_parametrisation .eq. 'HS_Basilakos') then

       derivative_matter_perturbation = -(3*H0**2*Omega_m*y1 + (-3*H0**2*Omega_m + &
            3*X*conformal_Hubble_parameter(X)**2)*y2 + &
            2*X*(wavenumber_k**2*y5 + conformal_Hubble_parameter(X)*(y3 + 3*conformal_Hubble_parameter(X)*&
            y6)))/(2.*X**2*conformal_Hubble_parameter(X)**2)

    End if

  end function derivative_matter_perturbation

!!$  function f_sigma_8(a,y1,y2,y3,y5,y6)
!!$
!!$    use fiducial
!!$
!!$    Implicit none 
!!$
!!$    Real*8 :: f_sigma_8,a,y1,y2,y3,y5,y6
!!$    Real*8,parameter :: sigma_8 = 8.d-1
!!$
!!$    f_sigma_8 = sigma_8*a*derivative_matter_perturbation(a,y1,y2,y3,y5,y6)
!!$
!!$  end function f_sigma_8

  !####################
  ! TESTING SUBROUTINES
  !####################

  subroutine test_approximations()

    use fiducial

    Implicit none

    Integer*4,parameter :: number_points = 10000
    Integer*4 :: index
    Real*8,dimension(number_points) :: scale_factor

    open(UNIT_TEST,file='./output/test.txt')

    write(UNIT_TEST,*) '# a  \rho_DE  \rho_LAMBDA  \rho_DE_prime  W_DE  W_LAMBDA  W_DE_prime  ca2  cs2  c  ceff2'&
         ' \delta_m  V_m  \delta_de  V_de  \pi_de/\delta_de  '&
         !\delta P_de  \delta \rho_de'&
         !  \delta P_de/\rho_de'&
         ' \sigma F FR H'

    Do index=1,number_points

       scale_factor(index) = 10**(log10(initial_scale_factor) + real(index-1)*(log10(final_scale_factor) - &
            log10(initial_scale_factor))/real(number_points-1))

       write(UNIT_TEST,99) scale_factor(index),dark_energy_density(scale_factor(index)),3.d0*H0**2*(1.d0-Omega_m),&
            derivative_dark_energy_density(scale_factor(index)),equation_of_state(scale_factor(index)),-1.d0,&
            derivative_equation_of_state(scale_factor(index)),adiabatic_sound_speed_squared(scale_factor(index)),&
            sound_speed_squared(scale_factor(index)),speedL,effective_sound_speed_squared(scale_factor(index),wavenumber_k),&
            matter_density_perturbation(scale_factor(index),wavenumber_k),&
            matter_velocity_perturbation(scale_factor(index),wavenumber_k),&
            dark_energy_density_perturbation(scale_factor(index),matter_density_perturbation(scale_factor(index),wavenumber_k)),&
            dark_energy_velocity_perturbation(scale_factor(index),matter_density_perturbation(scale_factor(index),wavenumber_k)),&
            anisotropic_stress(scale_factor(index),0.d0,&
            1.d0,0.d0,0.d0),&
            !dark_energy_pressure_perturbation_over_dark_energy_density(scale_factor(index),&
            !matter_density_perturbation(scale_factor(index),wavenumber_k))*dark_energy_density(scale_factor(index)),&
            !dark_energy_density_perturbation(scale_factor(index),matter_density_perturbation(scale_factor(index),&
            !wavenumber_k))*dark_energy_density(scale_factor(index)),&
            !dark_energy_pressure_perturbation_over_dark_energy_density(scale_factor(index),&
            !matter_density_perturbation(scale_factor(index),wavenumber_k)),&
            sigma(scale_factor(index),0.d0,&
            1.d0,0.d0,0.d0),F_MG(scale_factor(index)),&
            FR(scale_factor(index)),conformal_Hubble_parameter(scale_factor(index))/scale_factor(index)

99     FORMAT(E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,&
            !E20.10,E20.10,&
            !E20.10,
            E20.10,E20.10,E20.10,E20.10)

    End do

    close(UNIT_TEST)

  end subroutine test_approximations

  subroutine test_function()

    use fiducial

    Implicit none

    Integer*4,parameter :: number_points = 1000
    Integer*4 :: index
    Real*8,dimension(number_points) :: scale_factor

    open(UNIT_TEST,file='./output/functions.txt')

    write(UNIT_TEST,*) '# a  w(a)  w_prime(a)  H(a)/H0  cs2(a)  ceff2(a)  Omega_m(a)  Omega_DE(a)'&
         '  pressure_perturbation_over_density  dm_th  Vm_th  dde_th  Vde_th  \pi(a)  Geff/GN  Qeff'

    Do index=1,number_points

       scale_factor(index) = 10**(log10(initial_scale_factor) + real(index-1)*(log10(final_scale_factor) - &
            log10(initial_scale_factor))/real(number_points-1))

       write(UNIT_TEST,99) scale_factor(index),&
            equation_of_state(scale_factor(index)),&
            derivative_equation_of_state(scale_factor(index)),&
            conformal_Hubble_parameter(scale_factor(index))/scale_factor(index)/H0,&
            sound_speed_squared(scale_factor(index)),&
            effective_sound_speed_squared(scale_factor(index),wavenumber_k),&
            Omega_m/scale_factor(index)**3,&
            Omega_DE(scale_factor(index)),&
            dark_energy_pressure_perturbation(scale_factor(index)),&
            matter_density_perturbation(scale_factor(index),wavenumber_k)/delta_0(wavenumber_k),&
            matter_velocity_perturbation(scale_factor(index),wavenumber_k)/delta_0(wavenumber_k),&
            dark_energy_density_perturbation(scale_factor(index),&
            matter_density_perturbation(scale_factor(index),wavenumber_k))/delta_0(wavenumber_k),&
            dark_energy_velocity_perturbation(scale_factor(index),&
            matter_density_perturbation(scale_factor(index),wavenumber_k))/delta_0(wavenumber_k),&
            anisotropic_stress(scale_factor(index),0.d0,0.d0,0.d0,0.d0),&
            Geff_over_GN(scale_factor(index)),&
            Qeff(scale_factor(index))
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
            !dark_energy_density(scale_factor(index))
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

99     FORMAT(ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,&
            ES20.10,ES20.10,ES20.10,ES20.10)!,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,&
!            ES20.10,ES20.10,ES20.10,ES20.10)

    End do

    close(UNIT_TEST)

  end subroutine test_function

  function ft(x, params) bind(c)
    
    !use fgls
    !use, intrinsic :: iso_c_binding
    !Implicit none
    real(c_double), value :: x
    type(c_ptr), value :: params
    real(c_double) :: ft

    ft = F_MG_prime(x)!conformal_Hubble_parameter(x)!x**1.5_c_double

  end function ft


  function derivative_f(x)

    !use fgls
    !use, intrinsic :: iso_c_binding
    !implicit none
    real(fgsl_double) :: result, abserr
    integer(fgsl_int) :: status
    type(fgsl_function) :: pwr
    real(fgsl_double) :: derivative_f,x


    pwr = fgsl_function_init(ft, c_null_ptr)

    status = fgsl_deriv_central (pwr, x , 1.E-8_fgsl_double, &
         result, abserr)

    derivative_f = result

  end function derivative_f

End module functions
