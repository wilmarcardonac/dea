Module perturbations
  
  use fgsl
  use, intrinsic :: iso_c_binding
  use background
  use fiducial
  use initial_conditions

  Implicit none

Contains

  !##################################
  ! SETTING INITIAL CONDITIONS STARTS
  !##################################

  subroutine set_initial_conditions()

    Implicit none

    Logical :: mode_is_subhorizon,mode_is_subsoundhorizon 

    write(UNIT_EXE_FILE,*) 'THE CONFORMAL HUBBLE PARAMETER H(a) AT INITIAL SCALE FACTOR ',&
         initial_scale_factor, ' IS : ', conformal_Hubble_parameter(initial_scale_factor), ' Mpc^{-1}'

    write(UNIT_EXE_FILE,*) 'THE WAVENUMBER K CORRESPONDING TO THE HORIZON AT INITIAL SCALE FACTOR ', &
         initial_scale_factor, ' IS : ', conformal_Hubble_parameter(initial_scale_factor)/speedL, ' Mpc^{-1}'

    write(UNIT_EXE_FILE,*) 'THE CONFORMAL HUBBLE PARAMETER H(a) AT FINAL SCALE FACTOR ',&
         final_scale_factor, ' IS : ', conformal_Hubble_parameter(final_scale_factor),' Mpc^{-1}'

    write(UNIT_EXE_FILE,*) 'THE WAVENUMBER K CORRESPONDING TO THE HORIZON AT FINAL SCALE FACTOR ', &
         final_scale_factor, ' IS : ', conformal_Hubble_parameter(final_scale_factor)/speedL, ' Mpc^{-1}'

    If (wavenumber_k .lt. conformal_Hubble_parameter(initial_scale_factor)/speedL) then

       write(UNIT_EXE_FILE,*) 'CURRENT MODE STARTS BEYOND THE HORIZON (SUPER-HORIZON)' 

       mode_is_subhorizon = .false.

       mode_is_subsoundhorizon = .false. 

       If (wavenumber_k .lt. conformal_Hubble_parameter(final_scale_factor)/speedL) then

          write(UNIT_EXE_FILE,*) 'CURRENT MODE ENDS BEYOND THE HORIZON (SUPER-HORIZON)'

       Else

          write(UNIT_EXE_FILE,*) 'CURRENT MODE ENDS INSIDE THE HORIZON (SUB-HORIZON)'

          If (MG_parametrisation .eq. 'GR_DE') then 

             write(UNIT_EXE_FILE,*) 'THE MODE CROSSES THE HORIZON (IN MATTER DOMINANCE) AT SCALE FACTOR : ', &
                  scale_factor_horizon_crossing(wavenumber_k) 

             If ( (wavenumber_k/conformal_Hubble_parameter(final_scale_factor))**2 .lt. &
                  (1.d0/effective_sound_speed_squared(final_scale_factor,wavenumber_k)) ) then 

                write(UNIT_EXE_FILE,*) 'THE MODE ENDS BEYOND THE EFFECTIVE SOUND HORIZON (SUPER-SOUND-HORIZON)'

             Else

                write(UNIT_EXE_FILE,*) 'THE MODE CROSSES THE EFFECTIVE SOUND HORIZON (IN MATTER DOMINANCE) AT SCALE FACTOR : ', &
                     scale_factor_effective_sound_horizon(initial_scale_factor,wavenumber_k)

             End if

          Else if (MG_parametrisation .eq. 'GR_LAMBDA') then

             write(UNIT_EXE_FILE,*) 'THE MODE CROSSES THE HORIZON (IN MATTER DOMINANCE) AT SCALE FACTOR : ', &
                  scale_factor_horizon_crossing(wavenumber_k) 

          Else if ( (MG_parametrisation .eq. 'Savvas') .and. (approach .eq. 'GI') ) then

             write(UNIT_EXE_FILE,*) 'THE MODE CROSSES THE HORIZON (IN MATTER DOMINANCE) AT SCALE FACTOR : ', &
                  scale_factor_horizon_crossing(wavenumber_k) 

          Else if (((MG_parametrisation .eq. 'HS_Basilakos') .or. (MG_parametrisation .eq. 'Starobinsky_Basilakos')) .or.&
               (MG_parametrisation .eq. 'Savvas')) then

             write(UNIT_EXE_FILE,*) 'THE MODE CROSSES THE HORIZON (IN MATTER DOMINANCE) AT SCALE FACTOR : ', &
                  scale_factor_horizon_crossing(wavenumber_k) 

             write(UNIT_EXE_FILE,*) 'CURRENT MG_PARAMETRISATION ASSUMES SUB-HORIZON APPROXIMATION '&
                  'AND CURRENT MODE VIOLATES THIS CONDITION. CODE STOPPED!'

             stop

          End if

       End if

    Else

       write(UNIT_EXE_FILE,*) 'CURRENT MODE STARTS INSIDE THE HORIZON (SUB-HORIZON)'

       mode_is_subhorizon = .true.

       If (MG_parametrisation .eq. 'GR_DE') then

          If ( (wavenumber_k/conformal_Hubble_parameter(initial_scale_factor))**2 .lt. &
               (1.d0/effective_sound_speed_squared(initial_scale_factor,wavenumber_k)) ) then 

             write(UNIT_EXE_FILE,*) 'CURRENT MODE STARTS BEYOND THE SOUND HORIZON (SUPER-SOUND-HORIZON)'

             mode_is_subsoundhorizon = .false.

          Else

             write(UNIT_EXE_FILE,*) 'CURRENT MODE STARTS INSIDE THE SOUND HORIZON (SUB-SOUND-HORIZON)'

             mode_is_subsoundhorizon = .true.

          End if

       Else if (MG_parametrisation .eq. 'GR_LAMBDA') then

          continue

       Else if (MG_parametrisation .eq. 'HS_BASILAKOS') then

          continue

       Else if (MG_parametrisation .eq. 'Savvas') then

          continue

       End if

       If (wavenumber_k .lt. conformal_Hubble_parameter(final_scale_factor)/speedL) then

          write(UNIT_EXE_FILE,*) 'CURRENT MODE ENDS BEYOND THE HORIZON (SUPER-HORIZON)'

       Else

          write(UNIT_EXE_FILE,*) 'CURRENT MODE ENDS INSIDE THE HORIZON (SUB-HORIZON)'

          If (MG_parametrisation .eq. 'GR_DE') then 

             write(UNIT_EXE_FILE,*) 'THE MODE CROSSES THE HORIZON (IN MATTER DOMINANCE) AT SCALE FACTOR : ', &
                  scale_factor_horizon_crossing(wavenumber_k) 

             If ( (wavenumber_k/conformal_Hubble_parameter(final_scale_factor))**2 .lt. &
                  (1.d0/effective_sound_speed_squared(final_scale_factor,wavenumber_k)) ) then 

                write(UNIT_EXE_FILE,*) 'THE MODE ENDS BEYOND THE EFFECTIVE SOUND HORIZON (SUPER-SOUND-HORIZON)'

             Else

                write(UNIT_EXE_FILE,*) 'THE MODE CROSSES THE EFFECTIVE SOUND HORIZON (IN MATTER DOMINANCE) AT SCALE FACTOR : ', &
                     scale_factor_effective_sound_horizon(initial_scale_factor,wavenumber_k)

             End if

             write(UNIT_EXE_FILE,*) 'THE MODE CROSSES THE EFFECTIVE SOUND HORIZON (IN MATTER DOMINANCE) AT SCALE FACTOR : ', &
                  scale_factor_effective_sound_horizon(initial_scale_factor,wavenumber_k)

          Else if (MG_parametrisation .eq. 'GR_LAMBDA') then

             write(UNIT_EXE_FILE,*) 'THE MODE CROSSES THE HORIZON (IN MATTER DOMINANCE) AT SCALE FACTOR : ', &
                  scale_factor_horizon_crossing(wavenumber_k) 

          Else if (MG_parametrisation .eq. 'HS_BASILAKOS') then

             write(UNIT_EXE_FILE,*) 'THE MODE CROSSES THE HORIZON (IN MATTER DOMINANCE) AT SCALE FACTOR : ', &
                  scale_factor_horizon_crossing(wavenumber_k) 

          Else if (MG_parametrisation .eq. 'Savvas') then

             write(UNIT_EXE_FILE,*) 'THE MODE CROSSES THE HORIZON (IN MATTER DOMINANCE) AT SCALE FACTOR : ', &
                  scale_factor_horizon_crossing(wavenumber_k) 

          End if

       End if

    End if

!!$    If ( wavenumber_k/ks .lt. lower_limit_ks ) then
!!$
!!$       write(UNIT_EXE_FILE,*) 'MODE WAVELENGTH IS GREATER THAN SOUND HORIZON AT MATTER-RADIATION EQUALITY'
!!$
!!$    Else if ( wavenumber_k/ks .gt. upper_limit_ks ) then
!!$
!!$       write(UNIT_EXE_FILE,*) 'MODE WAVELENGTH IS SMALLER THAN SOUND HORIZON AT MATTER-RADIATION EQUALITY'
!!$
!!$    Else
!!$
!!$       write(UNIT_EXE_FILE,*) 'MODE DOES NOT SATISFY CONDITION FOR INITIAL POTENTIAL IN 1112.4837'
!!$
!!$       stop
!!$
!!$    End if

    If (approach .eq. 'GI') then 

       write(UNIT_EXE_FILE,*) 'EQUATIONS FOR PERTURBATIONS ARE WRITTEN IN A GAUGE INVARIANT FORMALISM. INITIAL CONDITIONS'&
            'ARE SET DURING MATTER DOMINATION. WE FOLLOW 1508.04569'

    Else

       write(UNIT_EXE_FILE,*) 'EQUATIONS FOR PERTURBATIONS ARE WRITTEN IN THE NEWTONIAN GAUGE. INITIAL CONDITIONS'&
            'ARE SET DURING MATTER DOMINATION WHERE THE TWO GRAVITATIONAL POTENTIALS ARE EQUAL'

    End if

    write(UNIT_EXE_FILE,*) 'INITIAL VALUE FOR POTENTIAL IS: ',initial_condition_gravitational_potential(wavenumber_k)

    write(UNIT_EXE_FILE,*) 'INITIAL CONDITIONS ARE SET AT a : ', initial_scale_factor

    write(UNIT_EXE_FILE,*) 'SYSTEM OF DIFFERENTIAL EQUATIONS IS WRITTEN IN THE CODE AS FOLLOWS: '

    If ( (MG_parametrisation .eq. 'GR_DE') .and. (dimension_system_ode .eq. 4) ) then

       write(UNIT_EXE_FILE,*) 'Y(1) IS \delta_m, MATTER DENSITY PERTURBATION '

       initial_conditions_system(1) = matter_density_perturbation(initial_scale_factor,wavenumber_k) ! MATTER DENSITY PERTURBATIONS

       write(UNIT_EXE_FILE,*) 'Y(2) IS \delta_de, DARK ENERGY DENSITY PERTURBATION '

       If (mode_is_subhorizon) then 

          If (mode_is_subsoundhorizon) then

             initial_conditions_system(2) = dark_energy_density_perturbation_sub_sound_horizon(initial_scale_factor,wavenumber_k) ! DE VELOCITY PERTURBATIONS ON SUB-SOUND SCALES 
             !de_density_perturbation_sub_sound_horizon(wavenumber_k) ! DARK ENERGY DENSITY PERTURBATIONS ON SUB-SOUND HORIZON SCALES

          Else

             If ( (e_pi .eq. 0.d0) .and. (f_pi .eq. 0.d0) ) then

                initial_conditions_system(2) = dark_energy_density_perturbation_super_sound_horizon(initial_scale_factor,&
                     wavenumber_k) ! DE DENSITY PERTURBATIONS ON SUPER-SOUND SCALES

             Else

                write(UNIT_EXE_FILE,*) 'INITIAL CONDITIONS FOR MODE SUPER-SOUND HORIZON BUT SUB-HORIZON IN A MODEL '&
                     'WITH NON-VANISHING ANISOTROPIC STRESS ARE NOT YEY IMPLEMENTED. CODE STOPPED'

                stop

             End if

          End if

       Else 

          initial_conditions_system(2) = dark_energy_density_perturbation_super_horizon(wavenumber_k) ! DE DENSITY PERTURBATIONS ON SUPER-HORIZON SCALES
          !dark_energy_density_perturbation_super_sound_horizon(initial_scale_factor,wavenumber_k) ! DE DENSITY PERTURBATIONS ON SUPER-SOUND SCALES

       End if

       write(UNIT_EXE_FILE,*) 'Y(3) IS \theta_m, MATTER VELOCITY PERTURBATION '

       initial_conditions_system(3) = matter_velocity_perturbation(initial_scale_factor,wavenumber_k) ! MATTER VELOCITY PERTURBATIONS

       write(UNIT_EXE_FILE,*) 'Y(4) IS \theta_de, DARK ENERGY VELOCITY PERTURBATION '

       If (mode_is_subhorizon) then 

          If (mode_is_subsoundhorizon) then

             initial_conditions_system(4) = -dark_energy_velocity_perturbation_sub_sound_horizon(initial_scale_factor,&
                  wavenumber_k)*H0*Sqrt(Omega_m/initial_scale_factor)/(1.d0 + equation_of_state(initial_scale_factor)) ! DE VELOCITY PERTURBATIONS ON SUB-SOUND HORIZON SCALES

             !de_velocity_perturbation_sub_sound_horizon(initial_scale_factor,wavenumber_k) ! DE VELOCITY PERTURBATIONS ON SUB-SOUND HORIZON SCALES

          Else

             If ((e_pi .eq. 0.d0) .and. (f_pi .eq. 0.d0)) then

                initial_conditions_system(4) = dark_energy_velocity_perturbation_super_sound_horizon(initial_scale_factor,&
                     wavenumber_k) ! DE VELOCITY PERTURBATIONS SUPER-SOUND SCALES

             Else

                write(UNIT_EXE_FILE,*) 'INITIAL CONDITIONS FOR MODE SUPER-SOUND HORIZON BUT SUB-HORIZON IN A MODEL '&
                     'WITH NON-VANISHING ANISOTROPIC STRESS ARE NOT YET IMPLEMENTED. CODE STOPPED'

                stop
                
             End if

          End if

       Else 

          initial_conditions_system(4) = -dark_energy_velocity_perturbation_super_horizon(initial_scale_factor,wavenumber_k)*&
               H0*Sqrt(Omega_m/initial_scale_factor)/(1.d0 + equation_of_state(initial_scale_factor)) ! DE VELOCITY PERTURBATIONS ON SUPER-HORIZON SCALES
          !dark_energy_velocity_perturbation_super_sound_horizon(initial_scale_factor,wavenumber_k) ! DE VELOCITY PERTURBATIONS ON SUPER-SOUND SCALES

       End if

    Else if ( (MG_parametrisation .eq. 'GR_DE') .and. (dimension_system_ode .eq. 6) ) then

       write(UNIT_EXE_FILE,*) 'Y(1) IS \delta_m, MATTER DENSITY PERTURBATION '

       initial_conditions_system(1) = matter_density_perturbation(initial_scale_factor,wavenumber_k) ! MATTER DENSITY PERTURBATIONS

       write(UNIT_EXE_FILE,*) 'Y(2) IS \delta_de, DARK ENERGY DENSITY PERTURBATION '

       If (mode_is_subhorizon) then 

          If (mode_is_subsoundhorizon) then

             initial_conditions_system(2) = dark_energy_density_perturbation_sub_sound_horizon(initial_scale_factor,wavenumber_k) ! DE VELOCITY PERTURBATIONS ON SUB-SOUND SCALES
             !de_density_perturbation_sub_sound_horizon(wavenumber_k) ! DARK ENERGY DENSITY PERTURBATIONS ON SUB-SOUND HORIZON SCALES

          Else

             If ( (e_pi .eq. 0.d0) .and. (f_pi .eq. 0.d0) ) then

                initial_conditions_system(2) = dark_energy_density_perturbation_super_sound_horizon(initial_scale_factor,&
                     wavenumber_k) ! DE DENSITY PERTURBATIONS ON SUPER-SOUND SCALES

             Else

                write(UNIT_EXE_FILE,*) 'INITIAL CONDITIONS FOR MODE SUPER-SOUND HORIZON BUT SUB-HORIZON IN A MODEL '&
                     'WITH NON-VANISHING ANISOTROPIC STRESS ARE NOT YEY IMPLEMENTED. CODE STOPPED'

                stop

             End if

          End if

       Else 

          initial_conditions_system(2) = dark_energy_density_perturbation_super_horizon(wavenumber_k) ! DE DENSITY PERTURBATIONS ON SUPER-HORIZON SCALES
          !dark_energy_density_perturbation_super_sound_horizon(initial_scale_factor,wavenumber_k) ! DE DENSITY PERTURBATIONS ON SUPER-SOUND SCALES

       End if

       write(UNIT_EXE_FILE,*) 'Y(3) IS V_m, MATTER VELOCITY PERTURBATION '

       initial_conditions_system(3) = matter_velocity_perturbation(initial_scale_factor,wavenumber_k)/wavenumber_k ! MATTER VELOCITY PERTURBATIONS

       write(UNIT_EXE_FILE,*) 'Y(4) IS V_de, DARK ENERGY VELOCITY PERTURBATION '

       If (mode_is_subhorizon) then 

          If (mode_is_subsoundhorizon) then

             initial_conditions_system(4) = -dark_energy_velocity_perturbation_sub_sound_horizon(initial_scale_factor,&
                  wavenumber_k)*H0*Sqrt(Omega_m/initial_scale_factor)/wavenumber_k ! DE VELOCITY PERTURBATIONS ON SUB-SOUND HORIZON SCALES
             !de_velocity_perturbation_sub_sound_horizon(initial_scale_factor,wavenumber_k) ! DARK ENERGY VELOCITY PERTURBATIONS ON SUB-SOUND SCALES

          Else

             If ((e_pi .eq. 0.d0) .and. (f_pi .eq. 0.d0)) then

                initial_conditions_system(4) = dark_energy_velocity_perturbation_super_sound_horizon(initial_scale_factor,&
                     wavenumber_k)*(1.d0 + equation_of_state(initial_scale_factor))/wavenumber_k ! DE VELOCITY PERTURBATIONS SUPER-SOUND SCALES

             Else

                write(UNIT_EXE_FILE,*) 'INITIAL CONDITIONS FOR MODE SUPER-SOUND HORIZON BUT SUB-HORIZON IN A MODEL '&
                     'WITH NON-VANISHING ANISOTROPIC STRESS ARE NOT YET IMPLEMENTED. CODE STOPPED'

                stop
                
             End if

          End if

       Else 

          initial_conditions_system(4) = -dark_energy_velocity_perturbation_super_horizon(initial_scale_factor,wavenumber_k)*&
               H0*Sqrt(Omega_m/initial_scale_factor)/wavenumber_k ! DE VELOCITY PERTURBATIONS ON SUPER-HORIZON SCALES
!dark_energy_velocity_perturbation_super_sound_horizon(initial_scale_factor,wavenumber_k) ! DE VELOCITY PERTURBATIONS ON SUPER-SOUND SCALES

       End if

       write(UNIT_EXE_FILE,*) 'Y(5) IS \phi, GRAVITATIONAL POTENTIAL '

       initial_conditions_system(5) = initial_condition_gravitational_potential(wavenumber_k)

       write(UNIT_EXE_FILE,*) 'Y(6) IS \psi, GRAVITATIONAL POTENTIAL '

       initial_conditions_system(6) = initial_condition_gravitational_potential(wavenumber_k)

    Else if ( ( ( (MG_parametrisation .eq. 'Savvas') .or. (MG_parametrisation .eq. 'GR_LAMBDA') ) .or. &
         ( (MG_parametrisation .eq. 'HS_Basilakos') .or. (MG_parametrisation .eq. 'Starobinsky_Basilakos')&
         ) ) .and. (approach .eq. 'CHI') ) then

       write(UNIT_EXE_FILE,*) 'Y(1) IS \delta_m, MATTER DENSITY PERTURBATION '

       initial_conditions_system(1) = matter_density_perturbation(initial_scale_factor,wavenumber_k) ! MATTER DENSITY PERTURBATIONS

       write(UNIT_EXE_FILE,*) 'Y(2) IS \theta_m, MATTER VELOCITY PERTURBATION '

       initial_conditions_system(2) = matter_velocity_perturbation(initial_scale_factor,wavenumber_k) ! MATTER VELOCITY PERTURBATIONS

       write(UNIT_EXE_FILE,*) 'Y(3) IS \phi_+ , POTENTIAL \phi_+ '

       initial_conditions_system(3) = initial_condition_gravitational_potential(wavenumber_k) ! \phi_+

       write(UNIT_EXE_FILE,*) 'Y(4) IS \chi, POTENTIAL \chi '

       initial_conditions_system(4) = 0.d0 ! \chi

    Else if ( ( (MG_parametrisation .eq. 'Savvas') .or. ( (MG_parametrisation .eq. 'HS_Basilakos') &
         .or. (MG_parametrisation .eq. 'Starobinsky_Basilakos') ) ) .and. (approach .eq. 'EF') ) then

       write(UNIT_EXE_FILE,*) 'Y(1) IS \delta_m, MATTER DENSITY PERTURBATION '

       initial_conditions_system(1) = matter_density_perturbation(initial_scale_factor,wavenumber_k) ! MATTER DENSITY PERTURBATIONS

       write(UNIT_EXE_FILE,*) 'Y(2) IS \delta_de, DARK ENERGY DENSITY PERTURBATION '

       initial_conditions_system(2) = dark_energy_density_perturbation(initial_scale_factor,initial_conditions_system(1))

       write(UNIT_EXE_FILE,*) 'Y(3) IS V_m, MATTER VELOCITY PERTURBATION '

       initial_conditions_system(3) = matter_velocity_perturbation(initial_scale_factor,wavenumber_k)!/wavenumber_k ! MATTER VELOCITY PERTURBATIONS

       write(UNIT_EXE_FILE,*) 'Y(4) IS V, DARK ENERGY VELOCITY PERTURBATION '

       initial_conditions_system(4) = dark_energy_velocity_perturbation(initial_scale_factor,initial_conditions_system(1))

       write(UNIT_EXE_FILE,*) 'Y(5) IS \phi, POTENTIAL \phi '

       initial_conditions_system(5) = initial_condition_gravitational_potential(wavenumber_k)
        
       write(UNIT_EXE_FILE,*) 'Y(6) IS \psi, POTENTIAL \psi '

       initial_conditions_system(6) = initial_condition_gravitational_potential(wavenumber_k)

    Else if ( ( (MG_parametrisation .eq. 'Savvas') .or. ( (MG_parametrisation .eq. 'HS_Basilakos') &
         .or. (MG_parametrisation .eq. 'Starobinsky_Basilakos') ) ) .and. (approach .eq. 'GI') ) then

       write(UNIT_EXE_FILE,*) 'Y(1) IS \Delta_m, MATTER DENSITY PERTURBATION '

       initial_conditions_system(1) = -matter_velocity_perturbation(initial_scale_factor,wavenumber_k)/&
            conformal_Hubble_parameter(initial_scale_factor) ! MATTER DENSITY PERTURBATIONS

       write(UNIT_EXE_FILE,*) 'Y(2) IS \Delta_de, DARK ENERGY DENSITY PERTURBATION '

       initial_conditions_system(2) = 0.d0

       write(UNIT_EXE_FILE,*) 'Y(3) IS Theta_m, MATTER VELOCITY PERTURBATION '

       initial_conditions_system(3) = 3.d0*conformal_Hubble_parameter(initial_scale_factor)*&
            matter_velocity_perturbation(initial_scale_factor,wavenumber_k)/wavenumber_k**2 ! MATTER VELOCITY PERTURBATIONS

       write(UNIT_EXE_FILE,*) 'Y(4) IS Theta_de, DARK ENERGY VELOCITY PERTURBATION '

       initial_conditions_system(4) = 0.d0

       write(UNIT_EXE_FILE,*) 'Y(5) IS Z '

       initial_conditions_system(5) = 3.d0*Omega_m*conformal_Hubble_parameter(initial_scale_factor)*&
            matter_velocity_perturbation(initial_scale_factor,wavenumber_k)/wavenumber_k**2/2.d0

    End if

  end subroutine set_initial_conditions

  !##################################
  ! SETTING INITIAL CONDITIONS ENDS
  !##################################

  !####################################
  ! WRITING ANALYTICAL SOLUTIONS STARTS
  !####################################

  subroutine write_analytical_solutions()

    Implicit none

    Integer*4 :: m

    Real*8 :: Z

    If ( MG_parametrisation .eq. 'GR_DE' ) then

       If (dimension_system_ode .eq. 4) then

          write(UNIT_EXE_FILE,*) 'WRITING ANALYTICAL SOLUTIONS FOR THE CURRENT MODE IN MATTER DOMINATED REGIME'

          write(UNIT_OUTPUT_FILE2,*) '# scale_factor    \delta_{de}^{sup-hor}    V_{de}^{sup-hor}   '//trim(' ')//&
               '\delta_{de}^{sub-sound}    V_{de}^{sub-sound}    \delta_m     V_m '

          Do m=1,101

          Z = 10**(log10(initial_scale_factor) + real(m-1)*(log10(final_scale_factor) - &
               log10(initial_scale_factor))/real(101-1))

             ! PERTURBATION VELOCITICIES BELOW FOLLOW FROM EQ. (3.3) IN MY PAPER, THAT IS, V_{de}.
             write(UNIT_OUTPUT_FILE2,89) Z,dark_energy_density_perturbation_super_horizon(wavenumber_k),&
                  -dark_energy_velocity_perturbation_super_horizon(Z,wavenumber_k)*H0*sqrt(Omega_m/Z)/wavenumber_k,&
                  dark_energy_density_perturbation_sub_sound_horizon(Z,wavenumber_k),&
                  -dark_energy_velocity_perturbation_sub_sound_horizon(Z,wavenumber_k)*H0*sqrt(Omega_m/Z)/wavenumber_k,&
                  matter_density_perturbation(Z,wavenumber_k),matter_velocity_perturbation(Z,wavenumber_k)/wavenumber_k

89           Format(E20.10,E20.10,E20.10,E20.10,E20.10,ES20.10,ES20.10)

          End do

       Else if (dimension_system_ode .eq. 6) then

          write(UNIT_EXE_FILE,*) 'WRITING ANALYTICAL SOLUTIONS FOR THE CURRENT MODE'

          write(UNIT_OUTPUT_FILE2,*) '# scale_factor    \delta_{de}^{sup-hor}    V_{de}^{sup-hor}   '//trim(' ')//&
               '\delta_{de}^{sub-sound}    V_{de}^{sub-sound}     \delta_m     V_m'

          Do m=1,101

          Z = 10**(log10(initial_scale_factor) + real(m-1)*(log10(final_scale_factor) - &
               log10(initial_scale_factor))/real(101-1))

             ! PERTURBATION VELOCITICIES BELOW FOLLOW FROM EQ. (3.3) IN MY PAPER, THAT IS, V_{de}.
             write(UNIT_OUTPUT_FILE2,88) Z,dark_energy_density_perturbation_super_horizon(wavenumber_k),&
                  -dark_energy_velocity_perturbation_super_horizon(Z,wavenumber_k)*H0*sqrt(Omega_m/Z)/wavenumber_k,&
                  dark_energy_density_perturbation_sub_sound_horizon(Z,wavenumber_k),&
                  -dark_energy_velocity_perturbation_sub_sound_horizon(Z,wavenumber_k)*H0*sqrt(Omega_m/Z)/wavenumber_k,&
                  matter_density_perturbation(Z,wavenumber_k),matter_velocity_perturbation(Z,wavenumber_k)/wavenumber_k

88           Format(E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,E20.10)

          End do

          write(UNIT_EXE_FILE,*) 'REQUIRES FURTHER IMPLEMENTATION: POTENTIALS MUST BE ADDED IN OUTPUT'

       End if

    Else if ( ( ( (MG_parametrisation .eq. 'Savvas') .or. (MG_parametrisation .eq. 'GR_LAMBDA') ) .or. &
         ( (MG_parametrisation .eq. 'HS_Basilakos') .or. (MG_parametrisation .eq. 'Starobinsky_Basilakos')&
         ) ) .and. (approach .eq. 'CHI') ) then

       write(UNIT_EXE_FILE,*) 'WRITING ANALYTICAL SOLUTIONS FOR THE CURRENT MODE IN MATTER DOMINANCE'

       write(UNIT_OUTPUT_FILE2,*) '# scale_factor    \delta_m    V_m    phi_plus    chi'

       Do m=1,101

          Z = 10**(log10(initial_scale_factor) + real(m-1)*(log10(final_scale_factor) - &
               log10(initial_scale_factor))/real(101-1))

          ! PERTURBATION VELOCITICIES BELOW FOLLOW FROM DRAFT
          write(UNIT_OUTPUT_FILE2,91) Z,matter_density_perturbation(Z,wavenumber_k),&
               matter_velocity_perturbation(Z,wavenumber_k),phi_plus(Z),chi(Z)

91        Format(E20.10,E20.10,E20.10,E20.10,E20.10)

       End do

       call system('cd figures; python plot_perturbations.py')

    Else if ( ( (MG_parametrisation .eq. 'Savvas') .or. ( (MG_parametrisation .eq. 'HS_Basilakos') .or. &
         (MG_parametrisation .eq. 'Starobinsky_Basilakos') ) ) .and. (approach .eq. 'EF') ) then

       write(UNIT_EXE_FILE,*) 'WRITING ANALYTICAL SOLUTIONS FOR THE CURRENT MODE'

       write(UNIT_OUTPUT_FILE2,*) '# scale_factor    \delta    V    \delta_m    V_m    phi    psi'!&
!            '    cs2_MD    cs2_MD_prime    angular_frequency '

       Do m=1,101

          Z = 10**(log10(initial_scale_factor) + real(m-1)*(log10(final_scale_factor) - &
               log10(initial_scale_factor))/real(101-1))

          ! PERTURBATION VELOCITICIES BELOW FOLLOW FROM DRAFT
          write(UNIT_OUTPUT_FILE2,90) Z,dark_energy_density_perturbation(Z,&
               matter_density_perturbation(Z,wavenumber_k)),&
               dark_energy_velocity_perturbation(Z,matter_density_perturbation(Z,wavenumber_k)),&
               matter_density_perturbation(Z,wavenumber_k),&
               matter_velocity_perturbation(Z,wavenumber_k),&!/wavenumber_k,&
               initial_condition_gravitational_potential(wavenumber_k),&
               initial_condition_gravitational_potential(wavenumber_k)!,&
!               sound_speed_squared(Z),&
!               sound_speed_squared_prime(Z),&
!               angular_frequency(Z)

90        Format(E20.10,E20.10,E20.10,ES20.10,ES20.10,ES20.10,ES20.10)!,ES20.10,ES20.10,ES20.10)

       End do

       call system('cd figures; python plot_perturbations_Savvas_EF.py')

    End if

  end subroutine write_analytical_solutions

  !####################################
  ! WRITING ANALYTICAL SOLUTIONS ENDS
  !####################################

  !#####################################
  ! EFFECTIVE SOUND SPEED SQUARED BEGINS 
  !#####################################

  function effective_sound_speed_squared(a,k)

    Implicit none

    Real*8 :: a,k,effective_sound_speed_squared

    If ( MG_parametrisation .eq. 'GR_DE' ) then

       effective_sound_speed_squared = sound_speed_squared(a) - 2.d0*f_pi/3.d0

    Else if (MG_parametrisation .eq. 'GR_LAMBDA') then 

       effective_sound_speed_squared = sound_speed_squared(a)

    Else if ( ( (MG_parametrisation .eq. 'Starobinsky_Basilakos') .or. &
         (MG_parametrisation .eq. 'HS_Basilakos') ) .or. (MG_parametrisation .eq. 'Savvas') ) then

       effective_sound_speed_squared = sound_speed_squared(a) - (2.d0/3.d0)*wavenumber_k**2*FR(a)/&
            (a**2*F_MG(a) - a**2*F_MG(a)**2 + wavenumber_k**2*(2.d0 - 3.d0*F_MG(a))*FR(a))

    End if

  end function effective_sound_speed_squared 

  subroutine check_effective_sound_speed_squared()

    Implicit none

    If (MG_parametrisation .eq. 'GR_DE' .and. (approach .eq. 'EF') ) then

       If ( effective_sound_speed_squared(initial_scale_factor,wavenumber_k) .ge. 0 ) then

          If ( effective_sound_speed_squared(initial_scale_factor,wavenumber_k) .gt. 1.d0) then

             write(UNIT_EXE_FILE,*) 'EFFECTIVE SOUND SPEED GREATER THAN ONE. CAUSALITY ISSUES'

             stop

          Else

             write(UNIT_EXE_FILE,*) 'EFFECTIVE SOUND SPEED SQUARED IS POSITIVE: ',&
                  effective_sound_speed_squared(initial_scale_factor,wavenumber_k)

          End if
          
       Else

          write(UNIT_EXE_FILE,*) 'EFFECTIVE SOUND SPEED SQUARED IS NEGATIVE: ',&
               effective_sound_speed_squared(initial_scale_factor,wavenumber_k)

          write(UNIT_EXE_FILE,*) 'AND WILL POSSIBLY LEAD TO INSTABILITIES IN THE PERTURBATIONS'

       End if

    Else if ( ( ( (MG_parametrisation .eq. 'Savvas')  .or. (MG_parametrisation .eq. 'HS_Basilakos') ) .or. &
         (MG_parametrisation .eq. 'Starobinsky_Basilakos')  ) .and. (approach .eq. 'EF') ) then

       If ( effective_sound_speed_squared(initial_scale_factor,wavenumber_k) .ge. 0 ) then

          write(UNIT_EXE_FILE,*) 'AT INITIAL SCALE FACTOR a = ', initial_scale_factor

          write(UNIT_EXE_FILE,*) 'AT ONSET SOUND SPEED SQUARED IN MATTER DOMINATED REGIME IS : ', &
          sound_speed_squared_in_matter_dominated_regime(initial_scale_factor)

          write(UNIT_EXE_FILE,*) 'AT ONSET SOUND SPEED SQUARED IS : ', sound_speed_squared(initial_scale_factor)

          write(UNIT_EXE_FILE,*) 'EFFECTIVE SOUND SPEED SQUARED IS : ', &
               effective_sound_speed_squared(initial_scale_factor,wavenumber_k)

       Else 

          write(UNIT_EXE_FILE,*) 'AT ONSET SOUND SPEED SQUARED IN MATTER DOMINATED REGIME IS : ', &
          sound_speed_squared_in_matter_dominated_regime(initial_scale_factor)

          write(UNIT_EXE_FILE,*) 'AT ONSET SOUND SPEED SQUARED IS : ', sound_speed_squared(initial_scale_factor)

          write(UNIT_EXE_FILE,*) 'AT ONSET EFFECTIVE SOUND SPEED SQUARED IS NEGATIVE: ',&
               effective_sound_speed_squared(initial_scale_factor,wavenumber_k)

          write(UNIT_EXE_FILE,*) 'AND WILL POSSIBLY LEAD TO INSTABILITIES IN THE PERTURBATIONS'

       End if

       If ( effective_sound_speed_squared(final_scale_factor,wavenumber_k) .ge. 0 ) then

          write(UNIT_EXE_FILE,*) 'AT FINAL SCALE FACTOR a = ', final_scale_factor

          write(UNIT_EXE_FILE,*) 'AT FINAL SCALE FACTOR  SOUND SPEED SQUARED IN MATTER DOMINATED REGIME IS : ', &
          sound_speed_squared_in_matter_dominated_regime(final_scale_factor)

          write(UNIT_EXE_FILE,*) 'AT FINAL TIME SOUND SPEED SQUARED IS : ', sound_speed_squared(final_scale_factor)

          write(UNIT_EXE_FILE,*) 'EFFECTIVE SOUND SPEED SQUARED IS : ', &
               effective_sound_speed_squared(final_scale_factor,wavenumber_k)

       Else

          write(UNIT_EXE_FILE,*) 'AT FINAL SCALE FACTOR  SOUND SPEED SQUARED IN MATTER DOMINATED REGIME IS : ', &
          sound_speed_squared_in_matter_dominated_regime(final_scale_factor)

          write(UNIT_EXE_FILE,*) 'AT FINAL TIME SOUND SPEED SQUARED IS : ', sound_speed_squared(final_scale_factor)

          write(UNIT_EXE_FILE,*) 'TODAY S EFFECTIVE SOUND SPEED SQUARED IS NEGATIVE: ', &
               effective_sound_speed_squared(final_scale_factor,wavenumber_k)

       End If

    End if

  end subroutine check_effective_sound_speed_squared

  !###################################
  ! EFFECTIVE SOUND SPEED SQUARED ENDS 
  !###################################

  !##################################
  ! POTENTIAL FOR CHI APPROACH BEGINS
  !##################################

  function chi(a)

    Implicit none 

    Real*8 :: a,chi

    chi = -2.d0*wavenumber_k**2*fMG_RR(a)*F_MG(a)*phi_plus(a)/( F_MG(a)*a**2 + &
         3.d0*wavenumber_k**2*fMG_RR(a) )

  end function chi

  function phi_plus(a)

    Implicit none

    Real*8 :: a,phi_plus

    phi_plus = -3.d0*H0**2*a**2/2.d0/wavenumber_k**2/F_MG(a)*Omega_Matter(a)*&
         matter_density_perturbation(a,wavenumber_k)

  end function phi_plus

  !################################
  ! POTENTIAL FOR CHI APPROACH ENDS
  !################################

  ! ADIABATIC SOUND SPEED FOR AN EFFECTIVE FLUID FROM F(R) PARAMETRISATION BY BASILAKOS ET AL. 

  function adiabatic_sound_speed_squared(a)

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

  ! PRESSURE PERTURBATION OVER DENSITY IN F(R) MODELS ON SUB-HORIZON SCALES 

  function pressure_perturbation_over_density(a)

    Implicit none

    Real*8 :: a,pressure_perturbation_over_density

    pressure_perturbation_over_density = (2.d0*wavenumber_k**2*FR(a)/a**2/F_MG(a))/(&
         3.d0*F_MG(a) + 9.d0*wavenumber_k**2*FR(a)/a**2 ) !+ &

!!$            ((a*conformal_Hubble_parameter(a)**2 + a**2*conformal_Hubble_parameter(a)*&
!!$            derivative_conformal_Hubble_parameter(a))*FR_prime(a) + &
!!$            a**2*conformal_Hubble_parameter(a)**2*FR_double_prime(a) )/(a**2*F_MG(a)  &
!!$            + 3.d0*wavenumber_k**2*FR(a)) !+ &
!!$
!!$            ( a*conformal_Hubble_parameter(a)**2*F_MG_prime(a) + a**2*conformal_Hubble_parameter(a)*&
!!$            derivative_conformal_Hubble_parameter(a)*F_MG_prime(a) + a**2*conformal_Hubble_parameter(a)**2*&
!!$            F_MG_double_prime(a))/(F_MG(a)*wavenumber_k**2 + 3.d0*wavenumber_k**4*FR(a)/a**2  )


  end function pressure_perturbation_over_density

  ! ANISOTROPIC STRESS FOR AN EFFECTIVE FLUID FROM F(R) PARAMETRISATION BY BASILAKOS ET AL.

  function f1(a)

    Implicit none

    Real*8 :: f1,a

    f1 = 0.d0 !FR(a)/F_MG(a)/(1.d0 - F_MG(a))

  end function f1

  function derivative_f1(a)

    Implicit none

    Real*8 :: derivative_f1,a

    derivative_f1 = 0.d0
!!$((-1.d0 + 2.d0*F_MG(a))*FR(a)*F_MG_prime(a) - &
!!$         (-1.d0 + F_MG(a))*F_MG(a)*FR_prime(a))/((-1.d0 + F_MG(a))**2*F_MG(a)**2)

  end function derivative_f1

  function f2(a)

    Implicit none

    Real*8 :: f2,a

    f2 = 0.d0 !(2.d0 - 3.d0*F_MG(a))*FR(a)/F_MG(a)/(1.d0 - F_MG(a))

  end function f2

  function derivative_f2(a)

    Implicit none

    Real*8 :: derivative_f2,a

    derivative_f2 = 0.d0 
!!$((-2 + 4*F_MG(a) - 3*F_MG(a)**2)*F_MG_prime(a)*FR(a) + &
!!$         F_MG(a)*(2 - 5*F_MG(a) + 3*F_MG(a)**2)*FR_prime(a))/&
!!$       ((-1 + F_MG(a))**2*F_MG(a)**2)

  end function derivative_f2

  function anisotropic_stress_sub_horizon(a)

    Implicit none

    Real*8 :: a,anisotropic_stress_sub_horizon

       anisotropic_stress_sub_horizon = wavenumber_k**2*FR(a)/(a**2*F_MG(a) - a**2*F_MG(a)**2 +&
            wavenumber_k**2*FR(a)*(2.d0 - 3.d0*F_MG(a)) )*dark_energy_density_perturbation(a,&
            matter_density_perturbation(a,wavenumber_k))

  end function anisotropic_stress_sub_horizon

  function anisotropic_stress_sub_horizon_MD(a)

    Implicit none

    Real*8 :: a,anisotropic_stress_sub_horizon_MD

       anisotropic_stress_sub_horizon_MD = wavenumber_k**2*FR(a)/(a**2*F_MG(a) - a**2*F_MG(a)**2 +&
            wavenumber_k**2*FR(a)*(2.d0 - 3.d0*F_MG(a)) )*dark_energy_density_perturbation_in_MD_for_f_of_R(a)

  end function anisotropic_stress_sub_horizon_MD

  function anisotropic_stress(a)

    Implicit none

    Real*8 :: a,anisotropic_stress

    If ( (MG_parametrisation .eq. 'HS_Basilakos') .or. (MG_parametrisation .eq. 'Starobinsky_Basilakos') ) then

       anisotropic_stress = wavenumber_k**2*FR(a)/(a**2*F_MG(a)**2 + &
            3.d0*wavenumber_k**2*FR(a)*F_MG(a) )

    Else if (MG_parametrisation .eq. 'GR_DE') then

       anisotropic_stress = -4.d10 !e_pi*Delta_matter(a,wavenumber_k,y1,y3) + (f_pi/(1.d0 + &
            !g_pi**2*conformal_Hubble_parameter(a)**2/wavenumber_k**2))*Delta_dark_energy(a,wavenumber_k,y2,y4)

    Else if (MG_parametrisation .eq. 'Savvas') then

       anisotropic_stress = (wavenumber_k**2*FR(a)/a**2/F_MG(a))/(1.d0 + &
            3.d0*wavenumber_k**2*FR(a)/a**2/F_MG(a))/F_MG(a)

!delta_0(wavenumber_k)*a*wavenumber_k**2*Omega_Matter(a)*fMG_RR(a)/&
!            F_MG(a)/Omega_DE(a)/(F_MG(a)*a**2 + 3.d0*wavenumber_k**2*fMG_RR(a))

    Else

       anisotropic_stress = -3.d10

    End if

  end function anisotropic_stress

  function anisotropic_stress_fr(a, params) bind(c)
    
    real(c_double), value :: a
    type(c_ptr), value :: params
    real(c_double) :: anisotropic_stress_fr

    anisotropic_stress_fr = anisotropic_stress(a)

  end function anisotropic_stress_fr

  function derivative_anisotropic_stress(a)

    real(fgsl_double) :: result, abserr
    integer(fgsl_int) :: status
    type(fgsl_function) :: pwr
    real(fgsl_double) :: derivative_anisotropic_stress,a

    If ( ( (MG_parametrisation .eq. 'HS_Basilakos') .or. (MG_parametrisation .eq. 'Savvas') ) .or. &
         (MG_parametrisation .eq. 'Starobinsky_Basilakos') ) then

       pwr = fgsl_function_init(anisotropic_stress_fr, c_null_ptr)

       status = fgsl_deriv_central (pwr, a , 1.E-8_fgsl_double, &
         result, abserr)

       derivative_anisotropic_stress = result

    Else if (MG_parametrisation .eq. 'GR') then

       derivative_anisotropic_stress = -3.d10 

    Else

       derivative_anisotropic_stress = -4.d10

    End if

  end function derivative_anisotropic_stress

  function dark_energy_pressure_perturbation(a)

    Implicit none 

    Real*8 :: dark_energy_pressure_perturbation,a

!    If ( a .lt. switch_off_pressure_perturbation_terms) then

    dark_energy_pressure_perturbation = ( 2.d0*wavenumber_k**2*FR(a)/a**2/F_MG(a) +  &
         3.d0*( 1.d0 + 5.d0*wavenumber_k**2*FR(a)/a**2/F_MG(a) )*(&
         a*conformal_Hubble_parameter(a)**2*F_MG_prime(a) + a**2*conformal_Hubble_parameter(a)*&
         derivative_conformal_Hubble_parameter(a)*F_MG_prime(a) + a**2*conformal_Hubble_parameter(a)**2*&
         F_MG_double_prime(a) )/wavenumber_k**2)/(1.d0 + 3.d0*wavenumber_k**2*FR(a)/a**2/F_MG(a) )/3.d0/F_MG(a)

!2.d0*wavenumber_k**2*FR(a)/(3.d0*a**2*F_MG(a)**2 + &
!            9.d0*wavenumber_k**2*FR(a)*F_MG(a) ) + &

!            ((a*conformal_Hubble_parameter(a)**2 + a**2*conformal_Hubble_parameter(a)*&
!            derivative_conformal_Hubble_parameter(a))*FR_prime(a) + &
!            a**2*conformal_Hubble_parameter(a)**2*FR_double_prime(a) )/(a**2*F_MG(a)  &
!            + 3.d0*wavenumber_k**2*FR(a)) + &

!            ( a*conformal_Hubble_parameter(a)**2*F_MG_prime(a) + a**2*conformal_Hubble_parameter(a)*&
!            derivative_conformal_Hubble_parameter(a)*F_MG_prime(a) + a**2*conformal_Hubble_parameter(a)**2*&
!            F_MG_double_prime(a))/(F_MG(a)*wavenumber_k**2 + 3.d0*wavenumber_k**4*FR(a)/a**2  )

!    Else

!       dark_energy_pressure_perturbation = 2.d0*wavenumber_k**2*FR(a)/(3.d0*a**2*F_MG(a)**2 + &
!            9.d0*wavenumber_k**2*FR(a)*F_MG(a) ) 

!    End if

  end function dark_energy_pressure_perturbation

  function second_dark_energy_pressure_perturbation(a,params) bind(c)

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


  ! SCALE FACTOR AT WHICH THE MODE CROSSES THE HUBBLE HORIZON IN A FLAT MATTER DOMINATED UNIVERSE

  function scale_factor_horizon_crossing(k)

    Implicit none

    Real*8 :: k,scale_factor_horizon_crossing

    scale_factor_horizon_crossing = H0**2*Omega_m/(speedL**2*k**2)

  end function scale_factor_horizon_crossing

  ! SCALE FACTOR AT WHICH THE MODE CROSSES THE EFFECTIVE SOUND HORIZON IN A FLAT MATTER DOMINATED UNIVERSE

  function scale_factor_effective_sound_horizon(a,k) 

    Implicit none

    Real*8 :: a,k,scale_factor_effective_sound_horizon

    scale_factor_effective_sound_horizon = H0**2*Omega_m/(effective_sound_speed_squared(a,k)*k**2)

  end function scale_factor_effective_sound_horizon

  ! GAUGE-INVARIANT COMOVING DENSITY CONTRAST (MATTER) 

  function Delta_matter(a,k,y1,y3) 

    Implicit none

    Real*8 :: a,k,y1,y3,Delta_matter

    If (MG_parametrisation .eq. 'GR_DE') then

       Delta_matter = y1 + 3.d0*conformal_Hubble_parameter(a)*y3/k**2

    Else

       Delta_matter = y1 + 3.d0*conformal_Hubble_parameter(a)*y3/k**2

    End if

  end function Delta_matter

  ! GAUGE-INVARIANT COMOVING DENSITY CONTRAST (DARK ENERGY) 

  function Delta_dark_energy(a,k,y2,y4) 

    Implicit none

    Real*8 :: a,k,y2,y4,Delta_dark_energy

    If (MG_parametrisation .eq. 'GR_DE') then

       Delta_dark_energy = y2 + 3.d0*conformal_Hubble_parameter(a)*(1.d0 + equation_of_state(a))*y4/k**2

    Else

       Delta_dark_energy = y2 + 3.d0*conformal_Hubble_parameter(a)*y4/k**2 ! CAREFUL WITH DIFFERENCES IN y4

    End if

  end function Delta_dark_energy

  ! DARK ENERGY ANISOTROPIC STRESS 

  function sigma(a,y1,y2,y3,y4)

    Implicit none

    Real*8 :: a,y1,y2,y3,y4,sigma,anisotropic_stress

    anisotropic_stress = e_pi*Delta_matter(a,wavenumber_k,y1,y3) + (f_pi/(1.d0 + &
         g_pi**2*conformal_Hubble_parameter(a)**2/wavenumber_k**2))*Delta_dark_energy(a,wavenumber_k,y2,y4)

    If (MG_parametrisation .eq. 'GR_DE') then

       sigma = 2.d0*anisotropic_stress/(3.d0*(1.d0+equation_of_state(a)))

    Else

       sigma = 2.d0*anisotropic_stress/(3.d0*(1.d0 + equation_of_state(a) ))

    End if

  end function sigma

  !####################################################################################################################################################
  ! EXPRESSIONS FOR THE POTENTIALS AND THEIR INITIAL CONDITIONS (SET IN MATTER DOMINANCE) FOLLOW SOLUTIONS IN THE PAPER BY GUILLERMO, MARTIN, AND LUKAS
  !####################################################################################################################################################


  ! POTENTIAL \Phi

  function phi(a,k,y1,y2,y3,y4)

    Implicit none

    Real*8 :: a,k,y1,y2,y3,y4,phi

    If (MG_parametrisation .eq. 'GR_DE') then

       phi = -((3.d0*H0**2)/(2.d0*k**2))*(Omega_m*Delta_matter(a,k,y1,y3)/a + &
            (1 - Omega_m)*Delta_dark_energy(a,k,y2,y4)/a**(1.d0 + 3.d0*equation_of_state(a)))

    Else

       phi = -(3.d0*H0**2/2.d0/k**2)*(Omega_m*Delta_matter(a,k,y1,y3)/a + &
            Omega_DE(a)*a**2*Delta_dark_energy(a,k,y2,y4) )

    End if

  end function phi

  ! POTENTIAL \Psi 

  function psi(a,k,y1,y2,y3,y4) 

    Implicit none

    Real*8 :: a,k,y1,y2,y3,y4,psi,anisotropic_stress

    anisotropic_stress = e_pi*Delta_matter(a,wavenumber_k,y1,y3) + (f_pi/(1.d0 + &
         g_pi**2*conformal_Hubble_parameter(a)**2/wavenumber_k**2))*Delta_dark_energy(a,wavenumber_k,y2,y4)

    If (MG_parametrisation .eq. 'GR_DE') then

       psi = -9.d0*H0**2*(1.d0+w0_fld)*sigma(a,y1,y2,y3,y4)*(1.d0-&
            Omega_m)/(2.d0*k**2*a**(1.d0+3.d0*equation_of_state(a))) + phi(a,k,y1,y2,y3,y4)

    Else

       psi = phi(a,k,y1,y2,y3,y4) - 3.d0*H0**2/k**2*anisotropic_stress*Omega_DE(a)*a**2

    End if

  end function psi

  !################################################################################################################
  ! SOLUTIONS FOR DARK ENERGY PERTURBATIONS IN MATTER DOMINATED REGIME WHEN CONSIDERING THE DARK ENERGY ANISOTROPIC
  ! STRESS MODEL DEFINED ABOVE
  !################################################################################################################

  ! DARK ENERGY DENSITY PERTURBATIONS ON SUPER-HORIZON SCALES 

  function dark_energy_density_perturbation_super_horizon(k)

    Implicit none

    Real*8 :: k,dark_energy_density_perturbation_super_horizon

    dark_energy_density_perturbation_super_horizon = delta_0(k)*3.d0*H0**2*Omega_m/k**2*(&
         4.d0*e_pi - 3.d0*(1.d0 + equation_of_state(1.d0)))/(4.d0*f_pi - 3.d0)

  end function dark_energy_density_perturbation_super_horizon

  ! DARK ENERGY VELOCITY PERTURBATIONS ON SUPER-HORIZON SCALES. NOTE THE VARIABLE USED HERE IS V_{de} = -(1+w)*\theta/conformal_H 

  function dark_energy_velocity_perturbation_super_horizon(a,k)

    Implicit none 

    Real*8 :: a,k,dark_energy_velocity_perturbation_super_horizon

    dark_energy_velocity_perturbation_super_horizon = a*delta_0(k)*(4.d0*e_pi - &
         3.d0*(1.d0 + equation_of_state(1.d0) ))/(4.d0*f_pi - 3.d0)

  end function dark_energy_velocity_perturbation_super_horizon

  ! DARK ENERGY DENSITY PERTURBATIONS ON SUB-SOUND HORIZON SCALES

  function dark_energy_density_perturbation_sub_sound_horizon(a,k)

    Implicit none

    Real*8 :: a,k,dark_energy_density_perturbation_sub_sound_horizon

    dark_energy_density_perturbation_sub_sound_horizon = a*delta_0(k)*(2.d0*e_pi/(3.d0*&
         effective_sound_speed_squared(a,k)) + 3.d0*H0**2*Omega_m*(1.d0 + &
         equation_of_state(a))/(2.d0*a*effective_sound_speed_squared(a,k)*k**2))

  end function dark_energy_density_perturbation_sub_sound_horizon

  ! DARK ENERGY VELOCITY PERTURBATIONS ON SUB-SOUND HORIZON SCALES

  function dark_energy_velocity_perturbation_sub_sound_horizon(a,k)

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

    Implicit none

    Real*8 :: a,B,k,d01

    B = initial_condition_gravitational_potential(k)

    d01 = -(2.d0*a**(1.d0-2.d0*f_pi)*(3.d0+7.d0*f_pi+4.d0*f_pi**2)*k**2*delta_0(k))/(9.d0*f_pi*H0**2*(1.d0-Omega_m))

  end function d01

  ! MATTER DENSITY PERTURBATION

  function dm3(a,k) 

    Implicit none 

    Real*8 ::  a,k,dm3

    dm3 = (d01(0.6d0,k)*a**(2.d0*f_pi)/2.d0)*( -(18.d0*a**2*f_pi*H0**4*(1.d0&
         -Omega_m)**2)/((2.d0+8.d0*f_pi**2/3.d0 + 14.d0*f_pi/3.d0)*k**4) - 3.d0*(6.d0*(2.d0+8.d0*f_pi**2 + 10.d0*f_pi)*f_pi + &
         (13.d0+16.d0*f_pi**2+32.d0*f_pi)*g_pi**2)*H0**2*(1.d0 - Omega_m)/((2.d0+192.d0*f_pi**2/9.d0 + 32.d0*f_pi**3/3.d0 + &
         38.d0*f_pi/3.d0)*g_pi**2*k**2))         

  end function dm3

  ! MATTER VELOCITY PERTURBATION

  function Vm3(a,k)

    Implicit none

    Real*8 :: a,k,Vm3

    Vm3 = -d01(0.6d0,k)*a**(2.d0*f_pi)*3.d0*f_pi*H0**2*(1.d0-Omega_m)/(k**2*(2.d0+8.d0*f_pi**2/3.d0 + 14.d0*f_pi/3.d0))

  end function Vm3

  ! DARK ENERGY DENSITY PERTURBATION 

  function dd3(a,k)

    Implicit none 

    Real*8 :: a,k,dd3

    dd3 = d01(0.6d0,k)*a**(2.d0*f_pi)*(3.d0+2.d0*f_pi)*H0**2*(1.d0-Omega_m)/(k**2*(1.d0+4.d0*f_pi/3.d0))

  end function dd3

  ! DARK ENERGY VELOCITY PERTURBATION 

  function Vd3(a,k)

    Implicit none 

    Real*8 :: a,k,Vd3

    Vd3 = d01(0.6d0,k)*a**(-2.d0+2.d0*f_pi)

  end function Vd3

  !###############################################################################################################
  ! SOLUTIONS FOR BOTH DARK ENERGY DENSITY AND BOTH DARK ENERGY VELOCITY PERTURBATIONS IN THE f(R) PARAMETRISATION
  ! BY BASILAKOS ET AL. IN MATTER DOMINATED REGIME
  !###############################################################################################################

  function dark_energy_density_perturbation(a,y1)

    Implicit none

    Real*8 :: dark_energy_density_perturbation,a,y1

    dark_energy_density_perturbation = ((1.d0 - F_MG(a) + &
         wavenumber_k**2*(2.d0 - 3.d0*F_MG(a))*FR(a)/a**2/F_MG(a))/(1.d0 + &
         3.d0*wavenumber_k**2*FR(a)/a**2/F_MG(a)))*Omega_Matter(a)*y1/F_MG(a)/Omega_DE(a)

!!$-delta_0(wavenumber_k)*b_fR*wavenumber_k**2*&
!!$         (1.d0-Omega_m)*a**5/3.d0/Omega_m**2/H0**2 - 52.d0*a**4*b_fR*delta_0(wavenumber_k)*&
!!$         (1.d0-Omega_m)/35.d0/Omega_m - 81.d0*a**2*b_fR*delta_0(wavenumber_k)*(1.d0-Omega_m)*&
!!$         (2.d0*Omega_m - 9.d0*a**3*b_fR*(1.d0-Omega_m))*H0**4/5.d0/wavenumber_k**4

         

  end function dark_energy_density_perturbation

  function dark_energy_velocity_perturbation(a,y1)

    Implicit none

    Real*8 :: dark_energy_velocity_perturbation,a,y1

    dark_energy_velocity_perturbation = ( 1.d0 + 6.d0*wavenumber_k**2*FR(a)/a**2/F_MG(a)  )/( 1.d0 + &
         3.d0*wavenumber_k**2*FR(a)/a**2/F_MG(a) )*( a*conformal_Hubble_parameter(a)*F_MG_prime(a)/F_MG(a)/2.d0 )*&
         Omega_Matter(a)*y1/Omega_DE(a)

!!$12.d0*delta_0(wavenumber_k)*(1.d0-Omega_m)*b_fR*a**(7.d0/2.d0)*H0/&
!!$         5.d0/Sqrt(Omega_m) + H0**3*108.d0*a**(5.d0/2.d0)*b_fR*delta_0(wavenumber_k)*(1.d0-Omega_m)*Sqrt(Omega_m)/&
!!$         13.d0/wavenumber_k**2 

         !(conformal_Hubble_parameter(a)*wavenumber_k**2*&
         !FR_prime(a)/F_MG(a)/a  &
         !+ a*conformal_Hubble_parameter(a)*F_MG_prime(a)/F_MG(a)/2.d0&
         !)/(1.d0 + 3.d0*wavenumber_k**2*FR(a)/a**2/F_MG(a))

         

  end function dark_energy_velocity_perturbation
  

  function derivative_matter_perturbation(X,y1,y2,y3,y4,y5,y6)

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
!!$    
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

    Implicit none

    Integer*4,parameter :: number_points = 10000
    Integer*4 :: index
    Real*8,dimension(number_points) :: scale_factor

    open(UNIT_TEST,file='./output/test.txt')

    write(UNIT_TEST,*) '# a  \rho_DE  \rho_LAMBDA  \rho_DE_prime  W_DE  W_LAMBDA  W_DE_prime  ca2  cs2  c  ceff2'&
         ' \delta_m  V_m  \delta_de  V_de  \pi_de/\delta_de  '&
         !\delta P_de  \delta \rho_de'&
         !  \delta P_de/\rho_de'&
         ' F FR H'

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
            anisotropic_stress(scale_factor(index)),&
            !dark_energy_pressure_perturbation_over_dark_energy_density(scale_factor(index),&
            !matter_density_perturbation(scale_factor(index),wavenumber_k))*dark_energy_density(scale_factor(index)),&
            !dark_energy_density_perturbation(scale_factor(index),matter_density_perturbation(scale_factor(index),&
            !wavenumber_k))*dark_energy_density(scale_factor(index)),&
            !dark_energy_pressure_perturbation_over_dark_energy_density(scale_factor(index),&
            !matter_density_perturbation(scale_factor(index),wavenumber_k)),&
            F_MG(scale_factor(index)),&
            FR(scale_factor(index)),conformal_Hubble_parameter(scale_factor(index))/scale_factor(index)

99     FORMAT(E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,&
            !E20.10,E20.10,&
            !E20.10,
            E20.10,E20.10,E20.10,E20.10)

    End do

    close(UNIT_TEST)

  end subroutine test_approximations

  subroutine test_function()

    Implicit none

    Integer*4,parameter :: number_points = 1000
    Integer*4 :: index
    Real*8,dimension(number_points) :: scale_factor

    open(UNIT_TEST,file='./output/functions.txt')

    write(UNIT_TEST,*) '# a  anisotropic_stress_MD  anisotropic_stress_sub_horizon_MD   anisotropic_stress_sub'&
         'cs2MD   cs2MD_sub_horizon'

!w(a)  w_prime(a)  H(a)/H0  cs2(a)  ceff2(a)  Omega_m(a)  Omega_DE(a)'&
!         '  pressure_perturbation_over_density  dm_th  Vm_th  dde_th  Vde_th  \pi(a)  Geff/GN  Qeff'&
!         '  ca2(a) x'

    Do index=1,number_points

       scale_factor(index) = 10**(log10(initial_scale_factor) + real(index-1)*(log10(final_scale_factor) - &
            log10(initial_scale_factor))/real(number_points-1))

       write(UNIT_TEST,99) scale_factor(index),&
            anisotropic_stress(scale_factor(index)),&
            anisotropic_stress_sub_horizon_MD(scale_factor(index)),&
            anisotropic_stress_sub_horizon(scale_factor(index)),&
            sound_speed_squared_in_matter_dominated_regime(scale_factor(index)),&            
            sound_speed_squared(scale_factor(index))
!!$            equation_of_state(scale_factor(index)),&
!!$            derivative_equation_of_state(scale_factor(index)),&
!!$            conformal_Hubble_parameter(scale_factor(index))/scale_factor(index)/H0,&
!!$            sound_speed_squared(scale_factor(index)),&
!!$            effective_sound_speed_squared(scale_factor(index),wavenumber_k),&
!!$            Omega_m/scale_factor(index)**3,&
!!$            Omega_DE(scale_factor(index)),&
!!$            dark_energy_pressure_perturbation(scale_factor(index)),&
!!$            matter_density_perturbation(scale_factor(index),wavenumber_k)/delta_0(wavenumber_k),&
!!$            matter_velocity_perturbation(scale_factor(index),wavenumber_k)/delta_0(wavenumber_k),&
!!$            dark_energy_density_perturbation(scale_factor(index),&
!!$            matter_density_perturbation(scale_factor(index),wavenumber_k))/delta_0(wavenumber_k),&
!!$            dark_energy_velocity_perturbation(scale_factor(index),&
!!$            matter_density_perturbation(scale_factor(index),wavenumber_k))/delta_0(wavenumber_k),&
!!$            anisotropic_stress(scale_factor(index)),&
!!$            Geff_over_GN(scale_factor(index)),&
!!$            Qeff(scale_factor(index)),&
!!$            adiabatic_sound_speed_squared(scale_factor(index)),&
            !F_MG_savvas_prime(scale_factor(index))
!            Lambda/(ricci_scalar(scale_factor(index))-3.d0*Lambda)
            !ricci_scalar_prime(scale_factor(index))
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

99     FORMAT(ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10)!,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,&
!            ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10)!,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,ES20.10,&
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

End module perturbations
