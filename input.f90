Module input

  use fiducial

  Implicit none

Contains

  subroutine test_input_parameters()

    Implicit none

    Logical :: dir_exist ! CHECK EXISTENCE OF FILES
    
    inquire(file='./output',exist=dir_exist)

    If (dir_exist) then

       continue

    Else

       call system('mkdir output')

    End if

    inquire(file='./figures',exist=dir_exist)

    If (dir_exist) then

       inquire(file='./figures/background',exist=dir_exist)

       If (dir_exist) then

          continue

       Else

          call system('mkdir figures/background')

       End if

       inquire(file='./figures/perturbations',exist=dir_exist)

       If (dir_exist) then

          continue

       Else

          call system('mkdir figures/perturbations')

       End if

    Else

       call system('mkdir figures')

       call system('mkdir figures/background')

       call system('mkdir figures/perturbations')

    End if

    open(UNIT_EXE_FILE,file=Execution_information)

    If (MG_parametrisation .eq. 'GR_LAMBDA') then

       open(UNIT_OUTPUT_FILE,file=NUMERICAL_SOLUTION_GR_LAMBDA)

    Else if (MG_parametrisation .eq. 'GR_DE') then

       open(UNIT_OUTPUT_FILE,file=NUMERICAL_SOLUTION_GR_DE)

    Else

       open(UNIT_OUTPUT_FILE,file=NUMERICAL_SOLUTION)

    End if

    open(UNIT_OUTPUT_FILE2,file=ANALYTICAL_SOLUTION)

    write(UNIT_EXE_FILE,*) 'STARTING ANALYSIS. COMMON PARAMETERS FOR CURRENT RUN ARE AS FOLLOWS: '

    write(UNIT_EXE_FILE,*) ' '

    If (Omega_m .gt. 1.d0) then

       write(UNIT_EXE_FILE,*) 'MATTER DENSITY PARAMETER IS GREATER THAN 1. ANALYSIS ASSUMES FLAT UNIVERSE'

       stop

    Else if (Omega_m .lt. 0.d0) then

       write(UNIT_EXE_FILE,*) 'MATTER DENSITY PARAMETER IS NEGATIVE'

       stop
       
    Else

       write(UNIT_EXE_FILE,*) 'MATTER DENSITY PARAMETER: ', Omega_m

    End if

    If (wavenumber_k .lt. 0.d0) then

       write(UNIT_EXE_FILE,*) 'WAVENUMBER MUST BE A POSITIVE NUMBER'

       stop

    Else

       If (approach .eq. 'GI') then 

          write(UNIT_EXE_FILE,*) 'DIMENSIONLESS WAVENUMBER IS: ', dimensionless_wavenumber_K

       Else

          write(UNIT_EXE_FILE,*) 'WAVENUMBER: ', wavenumber_k, ' Mpc^{-1}'

       End if

    End if

    write(UNIT_EXE_FILE,*) 'HUBBLE CONSTANT IS: ', H0, ' Mpc^{-1}' 

    If (MG_parametrisation .eq. 'GR_DE') then

       write(UNIT_EXE_FILE,*) ''

       If ( (dimension_system_ode .eq. 4) .and. (approach .eq. 'EF') ) then

          write(UNIT_EXE_FILE,*) 'THIS ANALYSIS ASSUMES GENERAL RELATIVITY AND DARK ENERGY IS MODELLED AS A FLUID: ' 

          write(UNIT_OUTPUT_FILE,*) '# scale_factor    \delta_m    \delta_de    V_m    V_de '

       Else if ( (dimension_system_ode .eq. 6) .and. (approach .eq. 'EF') ) then

          write(UNIT_OUTPUT_FILE,*) '# scale_factor    \delta_m    \delta_de    V_m    V_de    \phi    \psi'

          write(UNIT_EXE_FILE,*) 'THIS ANALYSIS ASSUMES GENERAL RELATIVITY AND DARK ENERGY IS MODELLED AS A FLUID: ' 

          write(UNIT_EXE_FILE,*) 'SYSTEM OF EQUATIONS INCLUDES POTENTIALS'

       Else

          write(UNIT_EXE_FILE,*) 'THIS PARAMETRISATION MUST USE THE EFFECTIVE FLUID APPROACH'

          write(UNIT_EXE_FILE,*) 'DIMENSION OF THE SYSTEM OF DIFFERENTIAL EQUATIONS MUST BE 4 OR 6:'

          write(UNIT_EXE_FILE,*) '\delta_m, V_m, \delta_de, V_de, \phi, \psi.'

          stop

       End if
          
       write(UNIT_EXE_FILE,*) 'CONSTANT EQUATION OF STATE: ', w0_fld

       write(UNIT_EXE_FILE,*) 'CONSTANT "SOUND SPEED SQUARED": ', cs2_fld 

       write(UNIT_EXE_FILE,*) 'ANISOSTROPIC STRESS GIVEN BY THE MODEL IN 1402.5993:'

       write(UNIT_EXE_FILE,*) 'DEA PARAMETER e_\pi : ', e_pi

       write(UNIT_EXE_FILE,*) 'DEA PARAMETER f_\pi : ', f_pi

       write(UNIT_EXE_FILE,*) 'DEA PARAMETER g_\pi : ', g_pi

    Else if (MG_parametrisation .eq. 'GR_LAMBDA') then

       write(UNIT_EXE_FILE,*) ''

       If (dimension_system_ode .eq. 4) then

          write(UNIT_OUTPUT_FILE,*) '# scale_factor    \delta_m    V_m    \phi_+    \chi'

          continue

       Else

          write(UNIT_EXE_FILE,*) 'DIMENSION OF THE SYSTEM OF DIFFERENTIAL EQUATIONS MUST BE 4:'

          write(UNIT_EXE_FILE,*) '\delta_m, V_m, \Phi_+, \chi.'

          stop

       End if

       If (approach .eq. 'EF') then

          write(UNIT_EXE_FILE,*) 'EFFECTIVE FLUID APPROACH IS NOT IMPLEMENTED FOR A COSMOLOGICAL CONSTANT'

          write(UNIT_EXE_FILE,*) 'BY DEFINITION, IT DOES NOT HAVE PERTURBATIONS'

          stop

       Else

          write(UNIT_EXE_FILE,*) 'USING MODIFICATIONS TO GENERAL RELATIVITY WRITTEN IN TERMS OF \PHI_+ AND \CHI'

       End if

       write(UNIT_EXE_FILE,*) 'THIS ANALYSIS ASSUMES GENERAL RELATIVITY AND A COSMOLOGICAL CONSTANT'

    Else if (MG_parametrisation .eq. 'Savvas') then

       write(UNIT_EXE_FILE,*) ''

       If ( (dimension_system_ode .eq. 6) .and. (approach .eq. 'EF') ) then
          
          write(UNIT_OUTPUT_FILE,*) '# scale_factor    \delta_m    \delta_de    V_m    V_de    \phi    \psi'

          write(UNIT_EXE_FILE,*) 'USING EFFECTIVE FLUID APPROACH'

       Else if ( (dimension_system_ode .eq. 5) .and. (approach .eq. 'GI') ) then

          write(UNIT_OUTPUT_FILE,*) '# scale_factor    \Delta_m    \Delta_de    Theta_m    Theta_de    Z'

          write(UNIT_EXE_FILE,*) 'USING EFFECTIVE FLUID APPROACH AND GAUGE INVARIANT FORMALISM'

       Else if ( (dimension_system_ode .eq. 4) .and. (approach .eq. 'CHI') ) then

          write(UNIT_OUTPUT_FILE,*) '# scale_factor    \delta_m    V_m    \phi_+    \chi'

          write(UNIT_EXE_FILE,*) 'USING MODIFICATIONS TO GENERAL RELATIVITY WRITTEN IN TERMS OF \PHI_+ AND \CHI'

       Else

          write(UNIT_EXE_FILE,*) 'DIMENSION OF THE SYSTEM OF DIFFERENTIAL EQUATIONS MUST BE 4, 5 OR 6:'

          stop

       End if

       write(UNIT_EXE_FILE,*) 'THIS ANALYSIS ASSUMES A MODIFICATION TO GENERAL RELATIVITY GIVEN BY THE f(R)'

       write(UNIT_EXE_FILE,*) 'PARAMETRISATION IN 1309.1055. THE EQUATION OF STATE DOES NOT EVOLVE WITH TIME; IT IS -1'

       write(UNIT_EXE_FILE,*) 'THE PARAMETER alpha IN THIS PARAMETRISATION IS CHOSEN SO THAT THE CONDITION'

       write(UNIT_EXE_FILE,*) 'f_R(a=1) = ', fR0, ' BE SATISFIED'

    Else if (MG_parametrisation .eq. 'HS_Basilakos') then

       write(UNIT_EXE_FILE,*) ''

       If ( (dimension_system_ode .eq. 6) .and. (approach .eq. 'EF') ) then

          write(UNIT_OUTPUT_FILE,*) '# scale_factor    \delta_m    \delta_de    V_m    V_de    \phi    \psi  '

          write(UNIT_EXE_FILE,*) 'USING EFFECTIVE FLUID APPROACH'

       Else if ( (dimension_system_ode .eq. 4) .and. (approach .eq. 'CHI') ) then

          write(UNIT_OUTPUT_FILE,*) '# scale_factor    \delta_m    V_m    \phi_+    \chi'

          write(UNIT_EXE_FILE,*) 'USING MODIFICATIONS TO GENERAL RELATIVITY WRITTEN IN TERMS OF \PHI_+ AND \CHI'

       Else

          write(UNIT_EXE_FILE,*) 'DIMENSION OF THE SYSTEM OF DIFFERENTIAL EQUATIONS MUST BE 4 OR 6:'

          stop

       End if

       write(UNIT_EXE_FILE,*) 'THIS ANALYSIS ASSUMES A MODIFICATION TO GENERAL RELATIVITY GIVEN BY THE f(R)'

       write(UNIT_EXE_FILE,*) 'PARAMETRISATION IN 1302.6051 (HU & SAWICKI MODEL). THE EQUATION OF STATE EVOLVES WITH TIME'

       write(UNIT_EXE_FILE,*) 'THE PARAMETER b IN THIS PARAMETRISATION IS CHOSEN TO BE:'

       write(UNIT_EXE_FILE,*) 'b : ', b_fR

    Else if (MG_parametrisation .eq. 'Starobinsky_Basilakos') then

       write(UNIT_EXE_FILE,*) ''

       If ( (dimension_system_ode .eq. 6) .and. (approach .eq. 'EF') ) then

          write(UNIT_OUTPUT_FILE,*) '# scale_factor        \delta_m             \delta_de             V_m        '//trim(' ')//&
               'V_de            \phi                 \psi  \delta_m_prime'

          write(UNIT_EXE_FILE,*) 'USING EFFECTIVE FLUID APPROACH'

       Else if (dimension_system_ode .eq. 4) then

          write(UNIT_OUTPUT_FILE,*) '# scale_factor    \delta_m    V_m    \phi_+    \chi'

          write(UNIT_EXE_FILE,*) 'USING MODIFICATIONS TO GENERAL RELATIVITY WRITTEN IN TERMS OF \PHI_+ AND \CHI'

          write(UNIT_EXE_FILE,*) 'DIMENSION OF THE SYSTEM OF DIFFERENTIAL EQUATIONS MUST BE 6:'

          write(UNIT_EXE_FILE,*) '\delta_m, V_m, \delta_de, V_de, \Phi, \Psi'

          write(UNIT_EXE_FILE,*) 'IN THE CURRENT IMPLEMENTATION'

          stop

       Else

          write(UNIT_EXE_FILE,*) 'DIMENSION OF THE SYSTEM OF DIFFERENTIAL EQUATIONS MUST BE 4 OR 6:'

          stop

       End if

       write(UNIT_EXE_FILE,*) 'THIS ANALYSIS ASSUMES A MODIFICATION TO GENERAL RELATIVITY GIVEN BY THE f(R)'

       write(UNIT_EXE_FILE,*) 'PARAMETRISATION IN 1302.6051 (STAROBINSKY MODEL). THE EQUATION OF STATE EVOLVES WITH TIME'

       write(UNIT_EXE_FILE,*) 'THE PARAMETER b IN THIS PARAMETRISATION IS CHOSEN TO BE:'

       write(UNIT_EXE_FILE,*) 'b : ', b_fR

       write(UNIT_EXE_FILE,*) 'THIS MODEL IS NOT IMPLEMENTED YET'

       stop

    Else

       write(UNIT_EXE_FILE,*) 'MG_PARAMETRISATION IS NOT IMPLEMENTED. CHECK fiducial.f90. CODE WILL STOP'

       stop

    End if

    If (initial_scale_factor .gt. final_scale_factor) then

       write(UNIT_EXE_FILE,*) 'INITIAL SCALE FACTOR MUST BE LESS THAN FINAL SCALE FACTOR'

       stop

    End if
    
  end subroutine test_input_parameters

End Module input
