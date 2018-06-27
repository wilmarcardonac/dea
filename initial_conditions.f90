Module initial_conditions
  
  use fgsl
  use, intrinsic :: iso_c_binding
  use background
  use fiducial

  Implicit none

Contains

  !##########################################################################################################################################################
  !INITIAL CONDITIONS DURING MATTER DOMINANCE. THE EXPRESSIONS MATCH THOSE FOUND BY MARTIN AND DOMENICO (WITHOUT ANISOTROPIC STRESS) IN PRD 80, 083519 (2009)
  !##########################################################################################################################################################
  
  ! MATTER DENSITY PERTURBATIONS. THE RELATION BETWEEN THE CONSTANT 'B' AND \delta_0 IN DOMENICO'S PAPER IS AS FOLLOWS: \delta_0 = -2 B k^2/3/H0^2/\Omega_m 
  
  function matter_density_perturbation(a,k)

    Implicit none

    Real*8 ::  a,k,matter_density_perturbation

    matter_density_perturbation = delta_0(k)*(a + 3.d0*H0**2*Omega_m/k**2 )

  end function matter_density_perturbation

  ! MATTER VELOCITY PERTURBATION 

  function matter_velocity_perturbation(a,k)

    Implicit none

    Real*8 :: a,k,matter_velocity_perturbation

    matter_velocity_perturbation = -delta_0(k)*sqrt(a)*H0*Sqrt(Omega_m)

  end function matter_velocity_perturbation

  ! DARK ENERGY DENSITY PERTURBATIONS ON SUPER-SOUND HORIZON SCALES 

  function dark_energy_density_perturbation_super_sound_horizon(a,k) 

    Implicit none 

    Real*8 :: a,k,dark_energy_density_perturbation_super_sound_horizon

    dark_energy_density_perturbation_super_sound_horizon = delta_0(k)*(1.d0 + &
         equation_of_state(a))*a*(1.d0/(1.d0 - 3.d0*equation_of_state(a)) &
         + 3.d0*H0**2*Omega_m/k**2/a ) 

  end function dark_energy_density_perturbation_super_sound_horizon

  ! DARK ENERGY VELOCITY PERTURBATIONS ON SUPER-SOUND HORIZON SCALES. NOTE THIS IS \theta_{DE}; IN DOMENICO'S PAPER THE AUTHORS USED V_{DE} = (1+w)\theta_{DE} 

  function dark_energy_velocity_perturbation_super_sound_horizon(a,k)

    Implicit none

    Real*8 :: a,k,dark_energy_velocity_perturbation_super_sound_horizon

    dark_energy_velocity_perturbation_super_sound_horizon = -delta_0(k)*H0*Sqrt(a*Omega_m) 

  end function dark_energy_velocity_perturbation_super_sound_horizon

  ! DARK ENERGY DENSITY PERTURBATIONS ON SUB-SOUND HORIZON SCALES 

  function de_density_perturbation_sub_sound_horizon(k) 

    Implicit none

    Real*8 :: B,k,de_density_perturbation_sub_sound_horizon

    B = initial_condition_gravitational_potential(k)

    de_density_perturbation_sub_sound_horizon = -B*(1.d0 + &
         equation_of_state(initial_scale_factor))/sound_speed_squared(initial_scale_factor)

  end function de_density_perturbation_sub_sound_horizon

  ! DARK ENERGY VELOCITY PERTURBATION ON SUB-SOUND HORIZON SCALES (DOMENICO'S EXPRESSION DIVIDED BY (1+w) ).

  function de_velocity_perturbation_sub_sound_horizon(a,k) 

    Implicit none 

    Real*8 :: a,k,de_velocity_perturbation_sub_sound_horizon

    de_velocity_perturbation_sub_sound_horizon = -3.d0*conformal_Hubble_parameter(a)*&
         (sound_speed_squared(a) - equation_of_state(a))*de_density_perturbation_sub_sound_horizon(k)

  end function de_velocity_perturbation_sub_sound_horizon

  ! CONSTANT \delta_0 IN DOMENICO AND MARTIN'S PAPER

  function delta_0(k) 

    Implicit none

    Real*8 :: B,k,delta_0

    B = initial_condition_gravitational_potential(k)

    delta_0 = -(2.d0*B*k**2)/(3.d0*H0**2*Omega_m) 

  end function delta_0

  ! CONSTANT FOR SUPER-HORIZON SOLUTIONS

  function delta_0_super_horizon(k)

    Implicit none

    Real*8 :: B,k,delta_0_super_horizon

    B = initial_condition_gravitational_potential(k)

    delta_0_super_horizon = -2.0d0*B*k**2*(1.0d0+w0_fld)*(4.0d0*f_pi-3.d0)/(3.d0*H0**2*Omega_m*&
         (4.d0*e_pi-3.d0*(1.d0+w0_fld)))

  end function delta_0_super_horizon

  ! INITIAL CONDITION FOR THE GRAVITATIONAL POTENTIAL IN MATTER DOMINATED ERA. DIMENSIONLESS

  function initial_condition_gravitational_potential(k) 

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

    Implicit none

    Real*8 :: k, primordial_dimensionless_power_spectrum

    primordial_dimensionless_power_spectrum = A_s*(k/kp)**(n_s- 1.d0)

  end function primordial_dimensionless_power_spectrum

  ! SOUND SPEED SQUARED FOR AN EFFECTIVE FLUID FROM F(R) PARAMETRISATION BY BASILAKOS ET AL. 

  function sound_speed_squared(a)

    Implicit none

    Real*8 :: a,sound_speed_squared

    If (MG_parametrisation .eq. 'GR_DE') then

       sound_speed_squared = cs2_fld

    Else if (MG_parametrisation .eq. 'GR_LAMBDA') then

       sound_speed_squared = 0.d0

    Else if ( ( (MG_parametrisation .eq. 'Starobinsky_Basilakos') .or. &
         (MG_parametrisation .eq. 'HS_Basilakos') ) .or. (MG_parametrisation .eq. 'Savvas') ) then
       
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

    End if

!!$    If ( a .lt. switch_off_pressure_perturbation_terms) then
!!$
!!$
!!$    Else
!!$
!!$       sound_speed_squared = (2.d0*wavenumber_k**2*FR(a))/(3.d0*a**2*F_MG(a) - &
!!$            3.d0*a**2*F_MG(a)**2 + 3.d0*wavenumber_k**2*(2.d0 - 3.d0*F_MG(a))*FR(a)) 
!!$
!!$    End if

  end function sound_speed_squared

end Module initial_conditions
