Program time_evolution

  !################
  ! REQUIRED MODULES
  !################

  use fgsl
  use fiducial 
  use functions
  
 
  !#################
  ! DEFINE VARIABLES 
  !#################

  IMPLICIT REAL*8 (A-H,O-Z)  

  Integer :: m   ! VARIABLE HOLDING INDEX IN LOOPS

  Logical :: dir_exist ! CHECK EXISTENCE OF FILES

  ! PARAMETERS FOR RADAU5 (FULL JACOBIAN)                             
  
  PARAMETER (ND=dimension_system_ode,LWORK=4*ND*ND+12*ND+20,LIWORK=3*ND+20) ! ND SEEMS TO BE THE NUMBER OF VARIABLES WHICH IN OUR CASE COINCIDES WITH 
  ! MATTER AND DARK ENERGY PERTURBATIONS: \delta_cdm, \delta_de, \theta_m, \theta_de
  DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK),RPAR(number_of_parameters)      ! Y IS AN ARRAY HOLDING THE PERTURBATION VARIABLES, RPAR HOLDS PARAMETERS 
  ! IN THE SYSTEM OF DIFFERENTIAL EQUATIONS 
  EXTERNAL RHSPER,JRHSPER,SOLOUT,MAS                         ! SUBROUTINES FROM THE SOLVER
       
  !########################################### 
  ! ASSIGMENTS AND INITIALIZATION OF VARIABLES
  !###########################################

  inquire(file='./output',exist=dir_exist)

  If (dir_exist) then

     continue

  Else

     call system('mkdir output')

  End if

  open(UNIT_EXE_FILE,file=Execution_information)

  open(UNIT_OUTPUT_FILE,file=NUMERICAL_SOLUTION)

  open(UNIT_OUTPUT_FILE2,file=ANALYTICAL_SOLUTION)

!  call test_function()
!  call test_approximations()

!  stop

  write(UNIT_EXE_FILE,*) 'STARTING ANALYSIS. PARAMETERS FOR CURRENT RUN ARE AS FOLLOWS: '

  If (MG_parametrisation .eq. 'GR') then

     write(UNIT_EXE_FILE,*) 'EQUATION OF STATE: ', w0_fld

     write(UNIT_EXE_FILE,*) 'MATTER DENSITY PARAMETER: ', Omega_m

     write(UNIT_EXE_FILE,*) 'SOUND SPEED SQUARED: ', cs2_fld 

     write(UNIT_EXE_FILE,*) 'WAVENUMBER: ', wavenumber_k, ' Mpc^{-1}'

     write(UNIT_EXE_FILE,*) 'DEA PARAMETER e_\pi : ', e_pi

     write(UNIT_EXE_FILE,*) 'DEA PARAMETER f_\pi : ', f_pi

     write(UNIT_EXE_FILE,*) 'DEA PARAMETER g_\pi : ', g_pi

  Else

     write(UNIT_EXE_FILE,*) 'EQUATION OF STATE EVOLVES WITH TIME'

     write(UNIT_EXE_FILE,*) 'MATTER DENSITY PARAMETER: ', Omega_m

     write(UNIT_EXE_FILE,*) 'SOUND SPEED SQUARED EVOLVES WITH TIME'

     write(UNIT_EXE_FILE,*) 'WAVENUMBER: ', wavenumber_k, ' Mpc^{-1}'

     write(UNIT_EXE_FILE,*) 'MG PARAMETER b : ', b_fR

     write(UNIT_EXE_FILE,*) 'DEA EVOLVES WITH TIME : '

  End if

  If (MG_parametrisation .eq. 'HS_Basilakos') then

     write(UNIT_EXE_FILE,*) 'CURRENT ANALYSIS USES HU & SAWICKI PARAMETRISATION FROM BASILAKOS ET AL.'

  Else if (MG_parametrisation .eq. 'Starobinsky_Basilakos') then

     write(UNIT_EXE_FILE,*) 'CURRENT ANALYSIS USES STAROBINSKY PARAMETRISATION FROM BASILAKOS ET AL.'

     write(UNIT_EXE_FILE,*) 'NOT IMPLEMENTED YET'

     stop

  Else if (MG_parametrisation .eq. 'GR') then

     write(UNIT_EXE_FILE,*) 'CURRENT ANALYSIS ASSUMES GENERAL RELATIVITY'

  Else

     write(UNIT_EXE_FILE,*) 'SELECTED MG PARAMETRISATION IS NOT IMPLEMENTED. CODE WILL STOP'

     stop

  End if

  If ( effective_sound_speed_squared(initial_scale_factor,wavenumber_k) .ge. 0 ) Then

     write(UNIT_EXE_FILE,*) 'AT INITIAL SCALE FACTOR a = ', initial_scale_factor

     write(UNIT_EXE_FILE,*) 'EFFECTIVE SOUND SPEED SQUARED FOR CURRENT MODEL IS : ', &
          effective_sound_speed_squared(initial_scale_factor,wavenumber_k)

     write(UNIT_EXE_FILE,*) 'AT FINAL SCALE FACTOR a = ', final_scale_factor

     write(UNIT_EXE_FILE,*) 'EFFECTIVE SOUND SPEED SQUARED FOR CURRENT MODEL IS : ', &
          effective_sound_speed_squared(final_scale_factor,wavenumber_k)

  Else

     write(UNIT_EXE_FILE,*) 'AT ONSET EFFECTIVE SOUND SPEED SQUARED FOR CURRENT MODEL IS NEGATIVE: ',&
          effective_sound_speed_squared(initial_scale_factor,wavenumber_k)

     write(UNIT_EXE_FILE,*) 'AND WILL POSSIBLY LEAD TO INSTABILITIES IN THE PERTURBATIONS'

     write(UNIT_EXE_FILE,*) 'TODAY EFFECTIVE SOUND SPEED SQUARED FOR CURRENT MODEL IS: ', &
          effective_sound_speed_squared(final_scale_factor,wavenumber_k)

  End If
  
  !##############################################################################################################################
  ! INPUT PARAMETERS FOR THE SUBROUTINE SOLVING THE SYSTEM. LOOK INSIDE 'radau5.f90 TO GET THE CORRECT MEANING FOR EACH PARAMETER
  !##############################################################################################################################

  N = dimension_system_ode ! DIMENSION OF THE SYSTEM 
                                 
  IJAC = 0 ! COMPUTE THE JACOBIAN -> 0: JACOBIAN IS INTERNALLY COMPUTED BY FINITE DIFFERENCES 
        
  MLJAC = N ! JACOBIAN IS A FULL MATRIX                                         

  IMAS = IMAS_RADAU ! DIFFERENTIAL EQUATION IS IN EXPLICIT FORM -> 0: M IS THE IDENTITY                         
        
  IOUT = 1 ! OUTPUT ROUTINE IS USED DURING INTEGRATION                         
         
  X = initial_scale_factor  ! INITIAL VALUE OF THE SCALE FACTOR
             
  write(UNIT_EXE_FILE,*) 'INITIAL CONDITIONS ARE SET AT SCALE FACTOR : ', X

  write(UNIT_EXE_FILE,*) 'SYSTEM OF DIFFERENTIAL EQUATIONS IS WRITTEN IN THE CODE AS FOLLOWS: '

  write(UNIT_EXE_FILE,*) 'Y(1) IS \delta_m, MATTER DENSITY PERTURBATION '

  write(UNIT_EXE_FILE,*) 'Y(2) IS \delta_de, DARK ENERGY DENSITY PERTURBATION '

  If (MG_parametrisation .eq. 'GR') then

     write(UNIT_EXE_FILE,*) 'Y(3) IS \theta_m, MATTER VELOCITY PERTURBATION '

     write(UNIT_EXE_FILE,*) 'Y(4) IS \theta_de, DARK ENERGY VELOCITY PERTURBATION '

  Else

     write(UNIT_EXE_FILE,*) 'Y(3) IS \theta_m, MATTER VELOCITY PERTURBATION '

     write(UNIT_EXE_FILE,*) 'Y(4) IS (1+wDE)*\theta_de, DARK ENERGY VELOCITY PERTURBATION '

     write(UNIT_EXE_FILE,*) 'Y(5) IS \phi, POTENTIAL \phi '

     write(UNIT_EXE_FILE,*) 'Y(6) IS \psi, POTENTIAL \psi '

  End if

  write(UNIT_EXE_FILE,*) 'THE CONFORMAL HUBBLE PARAMETER H(a) AT ', X, ' IS : ',conformal_Hubble_parameter(X), ' Mpc^{-1}'

  write(UNIT_EXE_FILE,*) 'THE WAVENUMBER K CORRESPONDING TO THE HORIZON AT INITIAL SCALE FACTOR ', X, ' IS : ',&
       conformal_Hubble_parameter(X)/speedL, ' Mpc^{-1}'

  If (wavenumber_k .lt. conformal_Hubble_parameter(X)/speedL) then

     write(UNIT_EXE_FILE,*) 'CURRENT MODE STARTS BEYOND THE HORIZON (SUPER-HORIZON)'

  Else

     write(UNIT_EXE_FILE,*) 'CURRENT MODE STARTS INSIDE THE HORIZON (SUB-HORIZON)'

  End if

  If ( wavenumber_k/ks .lt. lower_limit_ks ) then

     write(UNIT_EXE_FILE,*) 'MODE IS GREATER THAN SOUND HORIZON AT MATTER-RADIATION EQUALITY'

  Else if ( wavenumber_k/ks .gt. upper_limit_ks ) then

     write(UNIT_EXE_FILE,*) 'MODE IS SMALLER THAN SOUND HORIZON AT MATTER-RADIATION EQUALITY'

  Else

     write(UNIT_EXE_FILE,*) 'MODE DOES NOT SATISFY CONDITION FOR INITIAL POTENTIAL'

     stop

  End if

  write(UNIT_EXE_FILE,*) 'INITIAL VALUE FOR POTENTIAL IS: ',initial_condition_gravitational_potential(wavenumber_k)
       
  write(UNIT_EXE_FILE,*) 'THE CONFORMAL HUBBLE PARAMETER AT THE PRESENT TIME IS : ',conformal_Hubble_parameter(final_scale_factor),&
       ' Mpc^{-1}'

  write(UNIT_EXE_FILE,*) 'THE WAVENUMBER K CORRESPONDING TO THE HORIZON AT THE PRESENT TIME IS : ', &
       conformal_Hubble_parameter(final_scale_factor)/speedL, ' Mpc^{-1}'
  
  write(UNIT_EXE_FILE,*) 'THE WAVENUMBER K FOR THE CURRENT MODE IS : ', wavenumber_k, ' Mpc^{-1}'

  If (wavenumber_k .lt. conformal_Hubble_parameter(final_scale_factor)/speedL) then

     write(UNIT_EXE_FILE,*) 'CURRENT MODE ENDS BEYOND THE HORIZON (SUPER-HORIZON)'

  Else

     write(UNIT_EXE_FILE,*) 'CURRENT MODE ENDS INSIDE THE HORIZON (SUB-HORIZON)'

     write(UNIT_EXE_FILE,*) 'THE MODE CROSSES THE HORIZON (IN MATTER DOMINANCE) AT SCALE FACTOR : ', &
          scale_factor_horizon_crossing(wavenumber_k) 

     If (MG_parametrisation .eq. 'GR') then

        write(UNIT_EXE_FILE,*) 'THE MODE CROSSES THE EFFECTIVE SOUND HORIZON (IN MATTER DOMINANCE) AT SCALE FACTOR : ', &
             scale_factor_effective_sound_horizon(X,wavenumber_k)

     Else

        continue

     End if

  End if
  
  If (MG_parametrisation .eq. 'GR') then

     write(UNIT_OUTPUT_FILE,*) '# scale_factor        \delta_m             \delta_de             v_m        '//trim(' ')//&
          'v_de            \phi                 \psi'

  Else

     write(UNIT_OUTPUT_FILE,*) '# scale_factor        \delta_m             \delta_de             V_m        '//trim(' ')//&
          'V_de            \phi                 \psi  \delta_m_prime'

  End if
  
  !######################################################################################
  ! SETTING INITIAL CONDITIONS FOR THE PERTURBATIONS USING DOMENICO AND MARTIN'S SOLUTION
  !######################################################################################

  Y(1) = matter_density_perturbation(X,wavenumber_k) ! MATTER DENSITY PERTURBATIONS
  
  Y(3) = matter_velocity_perturbation(X,wavenumber_k) ! MATTER VELOCITY PERTURBATIONS

  If ( (effective_sound_speed_squared(X,wavenumber_k) .eq. 0.d0) .or. &
       ( (wavenumber_k/conformal_Hubble_parameter(X))**2 .lt. &
       abs( 1.d0/effective_sound_speed_squared(X,wavenumber_k) )  ) ) then

     If (MG_parametrisation .eq. 'GR') then

        Y(2) = dark_energy_density_perturbation_super_sound_horizon(X,wavenumber_k) ! DARK ENERGY DENSITY PERTURBATIONS ON SUPER-SOUND HORIZON SCALES

        Y(4) = dark_energy_velocity_perturbation_super_sound_horizon(X,wavenumber_k) ! DARK ENERGY VELOCITY PERTURBATIONS ON SUPER-SOUND HORIZON SCALES

     Else

        Y(2) = dark_energy_density_perturbation(X,Y(1)) ! DARK ENERGY DENSITY PERTURBATIONS ON SUPER-SOUND HORIZON SCALES
        !Y(2) = dark_energy_density_perturbation_super_sound_horizon(X,wavenumber_k) ! DARK ENERGY DENSITY PERTURBATIONS ON SUPER-SOUND HORIZON SCALES`
        !Y(4) = (1.d0 + equation_of_state(X))*dark_energy_velocity_perturbation_super_sound_horizon(X,wavenumber_k) ! DARK ENERGY VELOCITY PERTURBATIONS ON SUPER-SOUND HORIZON SCALES
        Y(4) = dark_energy_velocity_perturbation(X,Y(1)) ! DARK ENERGY VELOCITY PERTURBATIONS ON SUPER-SOUND HORIZON SCALES

     End if

     write(UNIT_EXE_FILE,*) 'CURRENT MODEL IS SUPER-SOUND-HORIZON AT THE INITIAL SCALE FACTOR'

  Else if ( (wavenumber_k/conformal_Hubble_parameter(X))**2 .gt. abs( 1.d0/effective_sound_speed_squared(X,wavenumber_k) ) ) Then

     If (MG_parametrisation .eq. 'GR') then

        Y(2) = de_density_perturbation_sub_sound_horizon(wavenumber_k) ! DARK ENERGY DENSITY PERTURBATIONS ON SUB-SOUND HORIZON SCALES

        Y(4) = de_velocity_perturbation_sub_sound_horizon(X,wavenumber_k) ! DARK ENERGY VELOCITY PERTURBATIONS ON SUB-SOUND HORIZON SCALES

     Else

        Y(2) = dark_energy_density_perturbation(X,Y(1)) ! DARK ENERGY DENSITY PERTURBATIONS ON SUB HORIZON SCALES

        !Y(4) = (1.d0 + equation_of_state(X))*de_velocity_perturbation_sub_sound_horizon(X,wavenumber_k) ! DARK ENERGY VELOCITY PERTURBATIONS ON SUB-SOUND HORIZON SCALES
        Y(4) = dark_energy_velocity_perturbation(X,Y(1)) ! DARK ENERGY VELOCITY PERTURBATIONS ON SUB HORIZON SCALES

     End if

     write(UNIT_EXE_FILE,*) 'CURRENT MODE IS SUB-SOUND-HORIZON AT THE INITIAL SCALE FACTOR'

  End If

  Y(5) = initial_condition_gravitational_potential(wavenumber_k)

  Y(6) = initial_condition_gravitational_potential(wavenumber_k)

  !################################################################################################################
  ! SOLUTIONS ON BOTH SUPER-HORIZON AND SUB-HORIZON SCALES FOR DARK ENERGY PERTURBATIONS IN MATTER DOMINATED REGIME
  !################################################################################################################

  If (MG_parametrisation .eq. 'GR') then

     If ( effective_sound_speed_squared(X,wavenumber_k) .ge. 0) then

        write(UNIT_EXE_FILE,*) 'WRITING ANALYTICAL SOLUTIONS FOR THE CURRENT MODE'

        write(UNIT_OUTPUT_FILE2,*) '# scale_factor    \delta_{de}^{sup-hor}    v_{de}^{sup-hor}   '//trim(' ')//&
             '\delta_{de}^{sub-sound}    v_{de}^{sub-sound}'

        Do m=1,101

           Z=10.d0**(-5.d0+Real(m-1)/100.d0*5.d0)

           ! PERTURBATION VELOCITICIES BELOW FOLLOW FROM EQ. (3.3) IN MY PAPER, THAT IS, v_{de}.
           write(UNIT_OUTPUT_FILE2,89) Z,dark_energy_density_perturbation_super_horizon(wavenumber_k),&
                -dark_energy_velocity_perturbation_super_horizon(Z,wavenumber_k)*H0*sqrt(Omega_m)/((1.d0+&
                w0_fld)*sqrt(Z)*wavenumber_k),&
                dark_energy_density_perturbation_sub_sound_horizon(Z,wavenumber_k),&
                -dark_energy_velocity_perturbation_sub_sound_horizon(Z,wavenumber_k)*H0*sqrt(Omega_m)/(wavenumber_k*(1.d0+&
                w0_fld)*sqrt(Z))

89         Format(E20.10,E20.10,E20.10,E20.10,E20.10)

        End do

     Else

        write(UNIT_EXE_FILE,*) 'EFFECTIVE SQUARE SOUND SPEED IS NEGATIVE. NOT ANALYTICAL SOLUTIONS WRITTEN'

     End If

  Else

        write(UNIT_EXE_FILE,*) 'WRITING ANALYTICAL SOLUTIONS FOR THE CURRENT MODE'

        write(UNIT_OUTPUT_FILE2,*) '# scale_factor  \delta_sub_horizon  V_sub_horizon  \delta_m  V_m  '&
             'phi  psi'

        Do m=1,100

           Z = 10**(log10(initial_scale_factor) + real(m-1)*(log10(final_scale_factor) - &
                log10(initial_scale_factor))/real(100-1))

           ! PERTURBATION VELOCITICIES BELOW FOLLOW FROM DRAFT
           write(UNIT_OUTPUT_FILE2,90) Z,dark_energy_density_perturbation(Z,&
                matter_density_perturbation(Z,wavenumber_k)),&
                dark_energy_velocity_perturbation(Z,&
                matter_density_perturbation(Z,wavenumber_k)),&
                matter_density_perturbation(Z,wavenumber_k),&
                matter_velocity_perturbation(Z,wavenumber_k),&
                phi(Z,wavenumber_k,matter_density_perturbation(Z,wavenumber_k),&
                dark_energy_density_perturbation(Z,&
                matter_density_perturbation(Z,wavenumber_k)),matter_velocity_perturbation(Z,wavenumber_k),&
                dark_energy_velocity_perturbation(Z,matter_density_perturbation(Z,wavenumber_k))),&
                psi(Z,wavenumber_k,matter_density_perturbation(Z,wavenumber_k),&
                dark_energy_density_perturbation(Z,&
                matter_density_perturbation(Z,wavenumber_k)),matter_velocity_perturbation(Z,wavenumber_k),&
                dark_energy_velocity_perturbation(Z,matter_density_perturbation(Z,wavenumber_k)))

90         Format(E20.10,E20.10,E20.10,ES20.10,ES20.10,ES20.10,ES20.10)

        End do

  End if

  !################################################################
  ! ENDPOINT INTEGRATION, TOLERANCE FOR SOLUTIONS, INITIAL STEPSIZE
  !################################################################

  XEND = final_scale_factor ! ENDPOINT OF INTEGRATION                                           

  RTOL = 1.d-14             ! REQUIRED TOLERANCE
  ATOL = 1.0d-10*RTOL       ! REQUIRED TOLERANCE 
  ITOL = 0                  ! REQUIRED TOLERANCE                                     
                                               
  H = 1.d-25                ! INITIAL STEP SIZE

  DO I=1,20                 ! SET DEFAULT VALUES                                                

     IWORK(I)=0 

     WORK(I)=0.D0 

  END DO

  !##############################################################
  ! CALL SUBROUTINE TO SOLVE THE SYSTEM OF DIFFERENTIAL EQUATIONS
  !##############################################################

!  CALL RADAU5(N,RHSPER,X,Y,XEND,H,RTOL,ATOL,ITOL,JRHSPER,IJAC,MLJAC,MUJAC,RHSPER,IMAS,MLMAS,MUMAS,&
!       SOLOUT,IOUT,WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)         

  CALL RADAU5(N,RHSPER,X,Y,XEND,H,RTOL,ATOL,ITOL,JRHSPER,IJAC,MLJAC,MUJAC,MAS,IMAS,MLMAS,MUMAS,&
       SOLOUT,IOUT,WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)         

  !##############################################################

  close(UNIT_EXE_FILE)

  close(UNIT_OUTPUT_FILE)

  close(UNIT_OUTPUT_FILE2)

End Program time_evolution                                          

SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N,RPAR,IPAR,IRTRN) 
  ! --- PRINTS SOLUTION AT EQUIDISTANT OUTPUT-POINTS BY USING "CONTR5"

  use functions
  use fiducial    

  IMPLICIT REAL*8 (A-H,O-Z)
  DIMENSION Y(N),CONT(LRC),RPAR(number_of_parameters) 
  COMMON /INTERN/XOUT 



  IF (NR.EQ.1) THEN 

     If (MG_parametrisation .eq. 'GR') then

        WRITE (UNIT_OUTPUT_FILE,99) X,Y(1),Y(2),Y(3)/wavenumber_k,Y(4)/wavenumber_k,&                        
             phi(X,wavenumber_k,Y(1),Y(2),Y(3),Y(4)),psi(X,wavenumber_k,Y(1),Y(2),Y(3),Y(4)),&
             derivative_matter_perturbation(X,Y(1),Y(2),Y(3),Y(4),Y(5),Y(6))

     Else

        WRITE (UNIT_OUTPUT_FILE,99) X,Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),&
             derivative_matter_perturbation(X,Y(1),Y(2),Y(3),Y(4),Y(5),Y(6))
        
     End if

     XOUT = initial_scale_factor 

  ELSE 

10   CONTINUE 

     IF (X.GE.XOUT) THEN 

        ! --- CONTINUOUS OUTPUT FOR RADAU5                                      

        If (MG_parametrisation .eq. 'GR') then

           WRITE (UNIT_OUTPUT_FILE,99) XOUT,CONTR5(1,XOUT,CONT,LRC),                &
                CONTR5(2,XOUT,CONT,LRC),                     &    
                CONTR5(3,XOUT,CONT,LRC)/wavenumber_k,                     &
                CONTR5(4,XOUT,CONT,LRC)/wavenumber_k, &                     
                phi(XOUT,wavenumber_k,CONTR5(1,XOUT,CONT,LRC),&
                CONTR5(2,XOUT,CONT,LRC),CONTR5(3,XOUT,CONT,LRC),CONTR5(4,XOUT,CONT,LRC)),&
                psi(XOUT,wavenumber_k,CONTR5(1,XOUT,CONT,LRC),CONTR5(2,XOUT,CONT,LRC),&
                CONTR5(3,XOUT,CONT,LRC),CONTR5(4,XOUT,CONT,LRC)),&
                derivative_matter_perturbation(XOUT,CONTR5(1,XOUT,CONT,LRC),&
                CONTR5(2,XOUT,CONT,LRC),CONTR5(3,XOUT,CONT,LRC),CONTR5(4,XOUT,CONT,LRC),&
                CONTR5(5,XOUT,CONT,LRC),CONTR5(6,XOUT,CONT,LRC))

        Else

           WRITE (UNIT_OUTPUT_FILE,99) XOUT,CONTR5(1,XOUT,CONT,LRC),CONTR5(2,XOUT,CONT,LRC),&    
                CONTR5(3,XOUT,CONT,LRC),CONTR5(4,XOUT,CONT,LRC),CONTR5(5,XOUT,CONT,LRC),&
                CONTR5(6,XOUT,CONT,LRC),derivative_matter_perturbation(XOUT,CONTR5(1,XOUT,CONT,LRC),&
                CONTR5(2,XOUT,CONT,LRC),CONTR5(3,XOUT,CONT,LRC),CONTR5(4,XOUT,CONT,LRC),&
                CONTR5(5,XOUT,CONT,LRC),CONTR5(6,XOUT,CONT,LRC))  

        End if

        If (XOUT .lt. 1.d-4) then

           XOUT = XOUT + 1.0d-6              

        Else if (XOUT .lt. 1.d-3) then
           
           XOUT = XOUT + 1.d-5 

        Else if (XOUT .lt. 1.d-2) then
           
           XOUT = XOUT + 1.d-4 
           
        Else if (XOUT .lt. 1.d-1) then

           XOUT = XOUT + 1.d-3 

        Else if (XOUT .le. 1.d0) then

           XOUT = XOUT + 1.d-3 
           
        End if

        GOTO 10 

     END IF

  END IF

99 FORMAT(E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,ES20.10)
 
  RETURN 

END SUBROUTINE SOLOUT                                   

subroutine RHSPER(N,X,Y,F,RPAR,IPAR)
           
  use fiducial
  use functions

  IMPLICIT REAL*8 (A-H,O-Z)
  Dimension Y(N),F(N),RPAR(number_of_parameters)
 
           
  !####################################################################################################
  ! THIS SUBROUTINE COMPUTES THE RIGHT HAND SIDE OF THE SYSTEM Y'=F(X,Y). IN OUR CASE IT CORRESPONDS TO:
  ! \delta_m' = ...,  \delta_de' = ..., \theta_m' = ...., \theta_de' = ....
  !####################################################################################################

  If (MG_parametrisation .eq. 'GR') then

     F(1) = (X**(-3.d0 - 6.d0*w0_fld)*(9.d0*H0**3*wavenumber_k**4*X**(1.d0 + 3.d0*w0_fld)*(-2.d0*(-1.d0 + &
          Omega_m)*e_pi + Omega_m*X**(3.d0*w0_fld))*(1.d0 + Omega_m*(-1.d0&
          + X**(3.d0*w0_fld)))*Y(1) - 9.d0*H0**3*(-1.d0 + Omega_m)*wavenumber_k**2*(1.d0 + Omega_m*(-1.d0 + &
          X**(3.d0*w0_fld)))*(wavenumber_k**2*(1.d0 + 2.d0*f_pi)*X**(1.d0 + 3.d0*w0_fld)&
          + 2.d0*H0**2*g_pi*(1.d0 + Omega_m*(-1.d0 + X**(3.d0*w0_fld))))*Y(2) + Sqrt(X**(-1.d0 - &
          3.d0*w0_fld)*(1.d0 + Omega_m*(-1.d0 + X**(3.d0*w0_fld))))*(wavenumber_k**2*X**(1.d0&
          + 3.d0*w0_fld)*(9.d0*H0**2*Omega_m*wavenumber_k**2*X**(1.d0 + 6.d0*w0_fld) - 2.d0*wavenumber_k**4*X**(2.d0 + &
          6.d0*w0_fld) + 27.d0*H0**4*(2.d0*e_pi - 2.d0*Omega_m*e_pi +&
          Omega_m*X**(3.d0*w0_fld))*(1.d0 + Omega_m*(-1.d0 + X**(3.d0*w0_fld))))*Y(3) - 9.d0*H0**2*(1.d0 + &
          w0_fld)*(-1.d0 + Omega_m)*(wavenumber_k**4*X**(2.d0 + 6.d0*w0_fld) +&
          3.d0*H0**2*wavenumber_k**2*(1.d0 + 2.d0*f_pi)*X**(1.d0 + 3.d0*w0_fld)*(1.d0 + Omega_m*(-1.d0 + X**(3.d0*w0_fld))) + &
          6.d0*H0**4*g_pi*(1.d0 + Omega_m*(-1.d0 +&
          X**(3.d0*w0_fld)))**2)*Y(4))))/(2.d0*H0*wavenumber_k**6*(1.d0 + Omega_m*(-1.d0 + X**(3.d0*w0_fld))))

     F(2) = (X**(-4.d0 - 9.d0*w0_fld)*(3.d0*H0*wavenumber_k**2*X**(1.d0 + 3.d0*w0_fld)*(6.d0*H0**4*(1.d0 + w0_fld)*(-1.d0 + &
          Omega_m)**2*g_pi + 2.d0*(w0_fld - &
          cs2_fld)*wavenumber_k**4*X**(2.d0 + 6.d0*w0_fld) - 3.d0*H0**2*(1.d0 + w0_fld)*(-1.d0 + &
          Omega_m)*X**(3.d0*w0_fld)*(2.d0*H0**2*Omega_m*g_pi + wavenumber_k**2*(1.d0 &
          + 2.d0*f_pi)*X))*Sqrt(X**(-1.d0 - 3.d0*w0_fld)*(1.d0 + Omega_m*(-1.d0 + X**(3.d0*w0_fld))))*Y(2) + &
          (1.d0 + w0_fld)*(9.d0*H0**3*wavenumber_k**4*X**(2.d0 + & 
          6.d0*w0_fld)*(-2.d0*(-1.d0 + Omega_m)*e_pi + Omega_m*X**(3.d0*w0_fld))*Sqrt(X**(-1.d0 - &
          3.d0*w0_fld)*(1.d0 + Omega_m*(-1.d0 + X**(3.d0*w0_fld))))*Y(1)&
          + 9.d0*H0**2*wavenumber_k**2*X**(1.d0 + 3.d0*w0_fld)*(Omega_m*wavenumber_k**2*X**(1.d0 + 6.d0*w0_fld) + &
          3.d0*H0**2*(-2.d0*(-1.d0 + Omega_m)*e_pi + Omega_m*X**(3.d0*w0_fld))*(1.d0&
          + Omega_m*(-1.d0 + X**(3.d0*w0_fld))))*Y(3) + (-2.d0*wavenumber_k**6*X**(3.d0 + 9.d0*w0_fld) + &
          9.d0*H0**2*wavenumber_k**4*X**(2.d0 + 6.d0*w0_fld)*((1.d0 - Omega_m)*(1.d0 +&
          3.d0*w0_fld - 2.d0*cs2_fld) + 2.d0*Omega_m*(w0_fld - cs2_fld)*X**(3.d0*w0_fld)) - &
          27.d0*H0**4*(1.d0 + w0_fld)*(-1.d0 + Omega_m)*wavenumber_k**2*(1.d0 + & 
          2.d0*f_pi)*X**(1.d0 + 3.d0*w0_fld)*(1.d0 + Omega_m*(-1.d0 + X**(3.d0*w0_fld))) - 54.d0*H0**6*(1.d0 + &
          w0_fld)*(-1.d0 + Omega_m)*g_pi*(1.d0 + &
          Omega_m*(-1.d0 + X**(3.d0*w0_fld)))**2)*Y(4))))/(2.d0*H0*wavenumber_k**6*Sqrt(X**(-1.d0 - &
          3.d0*w0_fld)*(1.d0 + Omega_m*(-1.d0 + X**(3.d0*w0_fld)))))


     F(3) = (X**(-3.d0 - 6.d0*w0_fld)*(-3.d0*H0*wavenumber_k**4*X**(1.d0 + 3.d0*w0_fld)*(2.d0*e_pi - &
          2.d0*Omega_m*e_pi + Omega_m*X**(3.d0*w0_fld))*Y(1) + &
          3.d0*H0*(-1.d0 + Omega_m)*wavenumber_k**2*(wavenumber_k**2*(1.d0 + 2.d0*f_pi)*X**(1.d0 + 3.d0*w0_fld) + &
          2.d0*H0**2*g_pi*(1.d0 + Omega_m*(-1.d0 + X**(3.d0*w0_fld))))*Y(2) - &
          Sqrt(X**(-1.d0 - 3.d0*w0_fld)*(1.d0 + Omega_m*(-1.d0 + X**(3.d0*w0_fld))))*(wavenumber_k**2*X**(1.d0 + &
          6.d0*w0_fld)*(9.d0*H0**2*Omega_m + 2.d0*wavenumber_k**2*X)*Y(3) + &
          18.d0*H0**4*(1.d0 + w0_fld)*(-1.d0 + Omega_m)**2*g_pi*Y(4) - 9.d0*H0**2*(-1.d0 + &
          Omega_m)*X**(3.d0*w0_fld)*(2.d0*wavenumber_k**2*e_pi*X*Y(3) + &
          (1.d0 + w0_fld)*(2.d0*H0**2*Omega_m*g_pi + wavenumber_k**2*(1.d0 + &
          2.d0*f_pi)*X)*Y(4)))))/(2.d0*wavenumber_k**4*Sqrt((Omega_m - (-1.d0 + Omega_m)/X**(3.d0*w0_fld))/X))

     F(4) = -(X**(-3.d0 - 6.d0*w0_fld)*(wavenumber_k**4*X**(1.d0 + 3.d0*w0_fld)*(4.d0*wavenumber_k**2*e_pi*X**(1.d0 + &
          3.d0*w0_fld) + 9.d0*H0**2*(1.d0 + w0_fld)*(2.d0*e_pi&
          - 2.d0*Omega_m*e_pi + Omega_m*X**(3.d0*w0_fld)))*Y(1) + wavenumber_k**2*(-2.d0*wavenumber_k**4*(3.d0*cs2_fld - &
          2.d0*f_pi)*X**(2.d0 + 6.d0*w0_fld) + H0**2*wavenumber_k**2*X**(1.d0&
          + 3.d0*w0_fld)*((1.d0 - Omega_m)*(9.d0*(1.d0 + w0_fld)*(1.d0 + 2.d0*f_pi) + 4.d0*g_pi) + &
          4.d0*Omega_m*g_pi*X**(3.d0*w0_fld)) - 18.d0*H0**4*(1.d0 + w0_fld)*(-1.d0&
          + Omega_m)*g_pi*(1.d0 + Omega_m*(-1.d0 + X**(3.d0*w0_fld))))*Y(2) + 3.d0*H0*Sqrt(X**(-1.d0 - &
          3.d0*w0_fld)*(1.d0 + Omega_m*(-1.d0 + X**(3.d0*w0_fld))))*(18.d0*H0**4*(1.d0&
          + w0_fld)**2*(-1.d0 + Omega_m)**2*g_pi*Y(4) + wavenumber_k**2*X**(1.d0 + 6.d0*w0_fld)*((9.d0*H0**2*(1.d0 + &
          w0_fld)*Omega_m + 4.d0*wavenumber_k**2*e_pi*X)*Y(3) + &
          2.d0*(1.d0 + w0_fld)*(2.d0*H0**2*Omega_m*g_pi + wavenumber_k**2*(1.d0 - 3.d0*cs2_fld + 2.d0*f_pi)*X)*Y(4)) - &
          H0**2*(1.d0 + w0_fld)*(-1.d0 + &
          Omega_m)*X**(3.d0*w0_fld)*(18.d0*wavenumber_k**2*e_pi*X*Y(3) + (18.d0*H0**2*(1.d0 + w0_fld)*Omega_m*g_pi + &
          wavenumber_k**2*(9.d0*(1.d0 + w0_fld)*(1.d0 + 2.d0*f_pi) +&
          4.d0*g_pi)*X)*Y(4)))))/(6.d0*H0*(1.d0 + w0_fld)*wavenumber_k**4*Sqrt((Omega_m - (-1.d0 + Omega_m)/X**(3.d0*w0_fld))/X))

     Else if (MG_parametrisation .eq. 'HS_Basilakos') then

        F(1) = -(3*H0**2*Omega_m*Y(1) + (-3*H0**2*Omega_m + 3*X*conformal_Hubble_parameter(X)**2)*Y(2) + &
          2*X*(wavenumber_k**2*Y(5) + conformal_Hubble_parameter(X)*(Y(3) + 3*conformal_Hubble_parameter(X)*&
          Y(6))))/(2.*X**2*conformal_Hubble_parameter(X)**2)

        F(2) = (2*H0**2*Omega_m*X*conformal_Hubble_parameter(X)*Y(4) - 2*X**2*conformal_Hubble_parameter(X)**3*Y(4) + &
         H0**2*Omega_m*(1 + equation_of_state(X))*(3*H0**2*Omega_m*(Y(1) - Y(2)) + 2*wavenumber_k**2*X*Y(5)) + &
         3*X**2*conformal_Hubble_parameter(X)**4*((-1 + equation_of_state(X))*Y(2) - 2*(1 + equation_of_state(X))*Y(6)) +& 
         X*conformal_Hubble_parameter(X)**2*(-3*H0**2*Omega_m*(1 + 2*dark_energy_pressure_perturbation(X) + &
         equation_of_state(X))*Y(1) + 6*H0**2*Omega_m*Y(2) - 2*(1 + equation_of_state(X))*(wavenumber_k**2*X*Y(5) - &
         3*H0**2*Omega_m*Y(6))))/(2.*X**2*conformal_Hubble_parameter(X)**2*(-(H0**2*Omega_m) + X*conformal_Hubble_parameter(X)**2)) 

        F(3) = -Y(3)/X + wavenumber_k**2*Y(6)/X/conformal_Hubble_parameter(X)

        F(4) = -((conformal_Hubble_parameter(X)*(1 - 3*equation_of_state(X))*Y(4) + (wavenumber_k**2*&
             ((H0**2*Omega_m*(-2*anisotropic_stress(X,0.d0,0.d0,0.d0,0.d0) + 3*dark_energy_pressure_perturbation(X))*Y(1))/&
             (H0**2*Omega_m - X*conformal_Hubble_parameter(X)**2) - 3*(1 + equation_of_state(X))*Y(6)))/3.)/&
             (X*conformal_Hubble_parameter(X)))

        F(5) = -(3*H0**2*Omega_m*Y(1) + (-3*H0**2*Omega_m + 3*X*conformal_Hubble_parameter(X)**2)*Y(2) + &
             2*wavenumber_k**2*X*Y(5) + 6*X*conformal_Hubble_parameter(X)**2*Y(6))/(6.*X**2*conformal_Hubble_parameter(X)**2)

        F(6) = (-3*H0**2*Omega_m*wavenumber_k**2*Y(1) + 3*wavenumber_k**2*(H0**2*Omega_m - &
             X*conformal_Hubble_parameter(X)**2)*Y(2) - 2*(9*X*conformal_Hubble_parameter(X)**4*(-2*anisotropic_stress(X,&
             0.d0,0.d0,0.d0,0.d0) + X*derivative_anisotropic_stress(X)) + 18*X**2*anisotropic_stress(X,0.d0,0.d0,0.d0,0.d0)*&
             conformal_Hubble_parameter(X)**3*derivative_conformal_Hubble_parameter(X) + wavenumber_k**4*X*Y(5) + &
             3*conformal_Hubble_parameter(X)**2*(3*H0**2*Omega_m*(3*anisotropic_stress(X,0.d0,0.d0,0.d0,0.d0) - &
             X*derivative_anisotropic_stress(X)) + wavenumber_k**2*X*(2*Y(5) - Y(6)))))/(6.*wavenumber_k**2*X**2*&
             conformal_Hubble_parameter(X)**2)

     Else if (MG_parametrisation .eq. 'Starobinsky_Basilakos') then

        write(UNIT_EXE_FILE,*) 'STAROBINSKY PARAMETRISATION IS NOT IMPLEMENTED YET'

        stop

     Else

        write(UNIT_EXE_FILE,*) 'UNKNOWN MG PARAMETRISATION'

        stop

     End if

end subroutine RHSPER                                          

subroutine JRHSPER(N,X,Y,DFY,LDFY,RPAR,IPAR)

           use functions
           use fiducial

           IMPLICIT REAL*8 (A-H,O-Z)
           Dimension Y(N),DFY(N,N),RPAR(number_of_parameters)

end subroutine JRHSPER

Subroutine MAS(N,AM,LMAS,RPAR,IPAR,X)
  
  use functions
  use fiducial 

  IMPLICIT REAL*8 (A-H,O-Z)
  Real*8,dimension(N,N) :: AM
  Real*8 :: X
  Dimension RPAR(number_of_parameters)
                        
  AM(1,1) = (9.d0*H0**2*Omega_m*conformal_Hubble_parameter(X))/(2.d0*wavenumber_k**2) + &
       X*conformal_Hubble_parameter(X)

  AM(1,2) = (9.d0*H0**2*(1.d0 - Omega_m)*conformal_Hubble_parameter(X))/&
       (2.d0*wavenumber_k**2*X**(3.d0*w0_fld))

  AM(1,3) = (27.d0*H0**2*Omega_m*conformal_Hubble_parameter(X)**2)/(2.d0*wavenumber_k**4)

  AM(1,4) = (27.d0*H0**2*(1.d0 - Omega_m)*(1.d0 + w0_fld)*conformal_Hubble_parameter(X)**2)/&
       (2.d0*wavenumber_k**4*X**(3.d0*w0_fld))

  AM(2,1) = (9.d0*H0**2*Omega_m*(1.d0 + w0_fld)*conformal_Hubble_parameter(X))/(2.d0*wavenumber_k**2)

  AM(2,2) = X*conformal_Hubble_parameter(X) + (9.d0*H0**2*(1.d0 - Omega_m)*(1.d0 + &
       w0_fld)*conformal_Hubble_parameter(X))/(2.d0*wavenumber_k**2*X**(3.d0*w0_fld))

  AM(2,3) = (27.d0*H0**2*Omega_m*(1.d0 + w0_fld)*conformal_Hubble_parameter(X)**2)/(2.d0*wavenumber_k**4)

  AM(2,4) = (27.d0*H0**2*(1.d0 - Omega_m)*(1.d0 + w0_fld)**2*conformal_Hubble_parameter(X)**2)/&
       (2.d0*wavenumber_k**4*X**(3.d0*w0_fld))

  AM(3,1) = 0.d0

  AM(3,2) = 0.d0

  AM(3,3) = X*conformal_Hubble_parameter(X)

  AM(3,4) = 0.d0

  AM(4,1) = 0.d0

  AM(4,2) = 0.d0

  AM(4,3) = 0.d0

  AM(4,4) = X*conformal_Hubble_parameter(X)

end subroutine MAS
