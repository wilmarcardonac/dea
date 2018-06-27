Program time_evolution

  !################
  ! REQUIRED MODULES
  !################

  use fgsl
  use fiducial
  use input
  use background
  use initial_conditions
  use perturbations
   
  !#################
  ! DEFINE VARIABLES 
  !#################

  IMPLICIT REAL*8 (A-H,O-Z)  

  Integer :: m   ! VARIABLE HOLDING INDEX IN LOOPS

  ! PARAMETERS FOR RADAU5 (FULL JACOBIAN)                             
  
  PARAMETER (ND=dimension_system_ode,LWORK=4*ND*ND+12*ND+20,LIWORK=3*ND+20) ! ND SEEMS TO BE THE NUMBER OF VARIABLES WHICH IN OUR CASE COINCIDES WITH 
  ! MATTER AND DARK ENERGY PERTURBATIONS: \delta_cdm, \delta_de, \theta_m, \theta_de. ALSO POTENTIALS: \Phi, \Psi, \Phi_+, \chi
  DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK),RPAR(number_of_parameters)      ! Y IS AN ARRAY HOLDING THE PERTURBATION VARIABLES, RPAR HOLDS PARAMETERS 
  ! IN THE SYSTEM OF DIFFERENTIAL EQUATIONS 
  EXTERNAL RHSPER,JRHSPER,SOLOUT,MAS                         ! SUBROUTINES FROM THE SOLVER
       
  !########################################### 
  ! ASSIGMENTS AND INITIALIZATION OF VARIABLES
  !###########################################

  call test_input_parameters()

  call check_effective_sound_speed_squared()

!  call test_function()
!  call test_approximations()
!  call compute_background()
!  stop
  !##############################################################################################################################
  ! INPUT PARAMETERS FOR THE SUBROUTINE SOLVING THE SYSTEM. LOOK INSIDE 'radau5.f90 TO GET THE CORRECT MEANING FOR EACH PARAMETER
  !##############################################################################################################################

  N = dimension_system_ode ! DIMENSION OF THE SYSTEM 
                                 
  IJAC = 0 ! COMPUTE THE JACOBIAN -> 0: JACOBIAN IS INTERNALLY COMPUTED BY FINITE DIFFERENCES 
        
  MLJAC = N ! JACOBIAN IS A FULL MATRIX                                         

  IMAS = IMAS_RADAU ! DIFFERENTIAL EQUATION IS IN EXPLICIT FORM -> 0: M IS THE IDENTITY                         
        
  IOUT = 1 ! OUTPUT ROUTINE IS USED DURING INTEGRATION                         

  X = initial_scale_factor  ! STARTING POINT INTEGRATION

  XEND = final_scale_factor ! ENDPOINT OF INTEGRATION                                           

  RTOL = 1.d-6              ! REQUIRED TOLERANCE
  ATOL = 1.d-6*RTOL!  0*RTOL ! REQUIRED TOLERANCE 
  ITOL = 0                   ! REQUIRED TOLERANCE                                     
                                               
  H = 1.d-8                  ! INITIAL STEP SIZE

  DO I=1,20                 ! SET DEFAULT VALUES                                                

     IWORK(I)=0 

     WORK(I)=0.D0 

  END DO

  !#################################################
  ! SETTING INITIAL CONDITIONS FOR THE PERTURBATIONS
  !#################################################

  call set_initial_conditions()

  Do m=1,dimension_system_ode

     Y(m) = initial_conditions_system(m)

  End do

  !################################################
  ! ANALYTICAL SOLUTIONS IN MATTER DOMINATED REGIME
  !################################################

  call write_analytical_solutions()

  !##############################################################
  ! CALL SUBROUTINE TO SOLVE THE SYSTEM OF DIFFERENTIAL EQUATIONS
  !##############################################################

!  CALL RADAU5(N,RHSPER,X,Y,XEND,H,RTOL,ATOL,ITOL,JRHSPER,IJAC,MLJAC,MUJAC,RHSPER,IMAS,MLMAS,MUMAS,&
!       SOLOUT,IOUT,WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)         

  CALL RADAU5(N,RHSPER,X,Y,XEND,H,RTOL,ATOL,ITOL,JRHSPER,IJAC,MLJAC,MUJAC,MAS,IMAS,MLMAS,MUMAS,&
       SOLOUT,IOUT,WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)         

  !##############
  ! CLOSING FILES
  !##############

  close(UNIT_EXE_FILE)

  close(UNIT_OUTPUT_FILE)

  close(UNIT_OUTPUT_FILE2)

End Program time_evolution                                          

!##################
! MAIN PROGRAM ENDS
!##################

!############################
! OUTPUTING SUBROUTINE STARTS
!############################

SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N,RPAR,IPAR,IRTRN) 
  ! --- PRINTS SOLUTION AT EQUIDISTANT OUTPUT-POINTS BY USING "CONTR5"

  use perturbations
  use fiducial    

  IMPLICIT REAL*8 (A-H,O-Z)
  DIMENSION Y(N),CONT(LRC),RPAR(number_of_parameters)
  COMMON /INTERN/XOUT 

  If (NR.EQ.1) THEN 

     If (MG_parametrisation .eq. 'GR_DE') then
        
        If (dimension_system_ode .eq. 6) then

           WRITE (UNIT_OUTPUT_FILE,98) X,Y(1),Y(2),Y(3)/wavenumber_k,Y(4)/wavenumber_k,Y(5),Y(6)

        Else if (dimension_system_ode .eq. 4) then
           
           WRITE (UNIT_OUTPUT_FILE,100) X,Y(1),Y(2),Y(3)/wavenumber_k,Y(4)/wavenumber_k                        

        End if

     Else if ( ( ( (MG_parametrisation .eq. 'Savvas') .or. (MG_parametrisation .eq. 'GR_LAMBDA') ) .or. &
         ( (MG_parametrisation .eq. 'HS_Basilakos') .or. (MG_parametrisation .eq. 'Starobinsky_Basilakos')&
         ) ) .and. (approach .eq. 'CHI') ) then
        
        WRITE (UNIT_OUTPUT_FILE,100) X,Y(1),Y(2)/wavenumber_k,Y(3),Y(4)

     Else if ( ( ( MG_parametrisation .eq. 'Savvas' ) .or. ( (MG_parametrisation .eq. 'HS_Basilakos') &
          .or. (MG_parametrisation .eq. 'Starobinsky_Basilakos') ) ) .and. (approach .eq. 'EF')  ) then

        WRITE (UNIT_OUTPUT_FILE,99) X,Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),&
             derivative_matter_perturbation(X,Y(1),Y(2),Y(3),Y(4),Y(5),Y(6))
        
     End if

     XOUT = initial_scale_factor 

  Else 

10   CONTINUE 

     IF (X.GE.XOUT) THEN 

        ! --- CONTINUOUS OUTPUT FOR RADAU5                                      

        If (MG_parametrisation .eq. 'GR_DE') then

           If (dimension_system_ode .eq. 6) then

              WRITE (UNIT_OUTPUT_FILE,98) XOUT,CONTR5(1,XOUT,CONT,LRC),CONTR5(2,XOUT,CONT,LRC),&    
                   CONTR5(3,XOUT,CONT,LRC)/wavenumber_k,CONTR5(4,XOUT,CONT,LRC)/wavenumber_k,&                     
                   CONTR5(5,XOUT,CONT,LRC),CONTR5(6,XOUT,CONT,LRC)

           Else if (dimension_system_ode .eq. 4) then

              WRITE (UNIT_OUTPUT_FILE,100) XOUT,CONTR5(1,XOUT,CONT,LRC),CONTR5(2,XOUT,CONT,LRC),&    
                   CONTR5(3,XOUT,CONT,LRC)/wavenumber_k,CONTR5(4,XOUT,CONT,LRC)/wavenumber_k
              
           End if
           
        Else if ( ( ( (MG_parametrisation .eq. 'Savvas') .or. (MG_parametrisation .eq. 'GR_LAMBDA') ) .or. &
         ( (MG_parametrisation .eq. 'HS_Basilakos') .or. (MG_parametrisation .eq. 'Starobinsky_Basilakos')&
         ) ) .and. (approach .eq. 'CHI') ) then

           WRITE (UNIT_OUTPUT_FILE,100) XOUT,CONTR5(1,XOUT,CONT,LRC),CONTR5(2,XOUT,CONT,LRC)/wavenumber_k,&    
                CONTR5(3,XOUT,CONT,LRC),CONTR5(4,XOUT,CONT,LRC)

        Else if ( ( ( MG_parametrisation .eq. 'Savvas' ) .or. ( (MG_parametrisation .eq. 'HS_Basilakos') &
          .or. (MG_parametrisation .eq. 'Starobinsky_Basilakos') ) ) .and. (approach .eq. 'EF')  ) then

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

98 FORMAT(E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,E20.10)

99 FORMAT(E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,ES20.10)

100 FORMAT(E20.10,E20.10,E20.10,E20.10,E20.10)
 
  RETURN 

END SUBROUTINE SOLOUT    

!##########################
! OUTPUTING SUBROUTINE ENDS
!##########################
                               
!############################################################################
! SUBROUTINE THAT COMPUTES RIGHT HAND SIDE THE THE SYSTEM OF EQUATIONS STARTS
!############################################################################

subroutine RHSPER(N,X,Y,F,RPAR,IPAR)
           
  use fiducial
  use background
  use perturbations

  IMPLICIT REAL*8 (A-H,O-Z)
  Dimension Y(N),F(N),RPAR(number_of_parameters)
 
           
  !####################################################################################################
  ! THIS SUBROUTINE COMPUTES THE RIGHT HAND SIDE OF THE SYSTEM Y'=F(X,Y). IN OUR CASE IT CORRESPONDS TO:
  ! \delta_m' = ...,  \delta_de' = ..., \theta_m' = ...., \theta_de' = ....
  !####################################################################################################

  If (MG_parametrisation .eq. 'GR_DE') then

     If (dimension_system_ode .eq. 4) then

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

     Else if (dimension_system_ode .eq. 6) then

        F(1) = -(3.d0*H0**2*X**2*(Omega_Matter(X)*Y(1) + Omega_DE(X)*Y(2)) + 2.d0*(wavenumber_k**2*Y(5) + &
             conformal_Hubble_parameter(X)*(Y(3) + 3.d0*conformal_Hubble_parameter(X)*Y(6))))/&
             (2.d0*X*conformal_Hubble_parameter(X)**2)

        F(2) = (-2.d0*wavenumber_k**2*conformal_Hubble_parameter(X)*(1.d0 + equation_of_state(X))*Y(4) - &
             18.d0*conformal_Hubble_parameter(X)**3*(1.d0 + equation_of_state(X))*(-equation_of_state(X) + &
             sound_speed_squared(X))*Y(4) + wavenumber_k**2*(1.d0 + equation_of_state(X))*&
             (-3.d0*H0**2*X**2*(Omega_Matter(X)*Y(1) + Omega_DE(X)*Y(2)) - 2.d0*wavenumber_k**2*Y(5)) -& 
             6.d0*wavenumber_k**2*conformal_Hubble_parameter(X)**2*((-equation_of_state(X) + &
             sound_speed_squared(X))*Y(2) + (1.d0 + equation_of_state(X))*Y(6)))/&
             (2.d0*wavenumber_k**2*X*conformal_Hubble_parameter(X)**2)

        F(3) = -Y(3)/X + wavenumber_k**2*Y(6)/X/conformal_Hubble_parameter(X)

        F(4) = -((-((wavenumber_k**2*sound_speed_squared(X)*Y(2))/(1.d0 + equation_of_state(X))) +& 
             conformal_Hubble_parameter(X)*(1.d0 - 3.d0*equation_of_state(X))*Y(4) + & 
             3.d0*conformal_Hubble_parameter(X)*(equation_of_state(X) - sound_speed_squared(X))*Y(4) - &
             wavenumber_k**2*((-2.d0*(e_pi*(Y(1) + (3.d0*conformal_Hubble_parameter(X)*Y(3))/wavenumber_k**2) + &
             (f_pi*(wavenumber_k**2*Y(2) + 3.d0*conformal_Hubble_parameter(X)*(1.d0 + equation_of_state(X))*Y(4)&
             ))/(wavenumber_k**2 + g_pi**2*conformal_Hubble_parameter(X)**2)))&
             /(3.d0*(1.d0 + equation_of_state(X))) + Y(6)))/(X*conformal_Hubble_parameter(X)))

        F(5) = -(3.d0*H0**2*X**2*(Omega_Matter(X)*Y(1) + Omega_DE(X)*Y(2)) + 2.d0*wavenumber_k**2*Y(5) + &
             6.d0*conformal_Hubble_parameter(X)**2*Y(6))/(6.d0*X*conformal_Hubble_parameter(X)**2)

        F(6) = (27.d0*H0**4*wavenumber_k**2*X**4*(wavenumber_k**2 + g_pi**2*conformal_Hubble_parameter(X)**2)*&
             (e_pi*g_pi**2*conformal_Hubble_parameter(X)**2 + wavenumber_k**2*(e_pi + f_pi + &
             f_pi*equation_of_state(X)))*Omega_DE(X)**2*Y(2) + 3.d0*H0**2*X**2*Omega_DE(X)*&
             (-(wavenumber_k**4*(wavenumber_k**4 + conformal_Hubble_parameter(X)**2*&
             (2.d0*wavenumber_k**2*(-6.d0*f_pi**2 + g_pi**2 + 9.d0*f_pi*equation_of_state(X)) + &
             g_pi**2*conformal_Hubble_parameter(X)*(-12.d0*f_pi*X*derivative_conformal_Hubble_parameter(X) + &
             conformal_Hubble_parameter(X)*(g_pi**2 + 18.d0*f_pi*equation_of_state(X)))))*Y(2)) + &
             3.d0*(6.d0*e_pi*g_pi**4*conformal_Hubble_parameter(X)**7*Y(3) - &
             6.d0*e_pi*g_pi**4*X*conformal_Hubble_parameter(X)**6*derivative_conformal_Hubble_parameter(X)*Y(3) + &
             2.d0*wavenumber_k**6*conformal_Hubble_parameter(X)*(e_pi*Y(3) + f_pi*(1.d0 + equation_of_state(X))*Y(4)) + &
             2.d0*wavenumber_k**4*conformal_Hubble_parameter(X)**3*(e_pi*(3.d0 + 6.d0*f_pi + 2.d0*g_pi**2)*Y(3) +& 
             f_pi*(3.d0 + 6.d0*f_pi + g_pi**2 - 9.d0*equation_of_state(X))*(1.d0 + equation_of_state(X))*Y(4)) + &
             2.d0*g_pi**2*wavenumber_k**2*conformal_Hubble_parameter(X)**5*(e_pi*(6.d0 + 6.d0*f_pi + g_pi**2)*Y(3) -& 
             3.d0*f_pi*(1.d0 + equation_of_state(X))*(-1.d0 + 3.d0*equation_of_state(X))*Y(4)) + &
             wavenumber_k**6*(e_pi + f_pi + f_pi*equation_of_state(X))*(3.d0*H0**2*X**2*Omega_Matter(X)*Y(1) + &
             2.d0*wavenumber_k**2*Y(5)) + g_pi**2*wavenumber_k**2*conformal_Hubble_parameter(X)**4*&
             (e_pi*(4.d0*f_pi*wavenumber_k**2 + 3.d0*g_pi**2*H0**2*X**2*Omega_Matter(X))*Y(1) +& 
             6.d0*X*derivative_conformal_Hubble_parameter(X)*(-2.d0*e_pi*Y(3) + f_pi*(1.d0 + &
             equation_of_state(X))*Y(4)) + 2.d0*e_pi*g_pi**2*wavenumber_k**2*Y(5)) + &
             wavenumber_k**4*conformal_Hubble_parameter(X)**2*((4.d0*e_pi*f_pi*wavenumber_k**2 + &
             3.d0*g_pi**2*H0**2*X**2*(2.d0*e_pi + f_pi + f_pi*equation_of_state(X))*Omega_Matter(X))*Y(1) - &
             6.d0*X*derivative_conformal_Hubble_parameter(X)*(e_pi*Y(3) + f_pi*(1.d0 + equation_of_state(X))*Y(4)) + &
             2.d0*g_pi**2*wavenumber_k**2*(2.d0*e_pi + f_pi + f_pi*equation_of_state(X))*Y(5)))) - &
             (wavenumber_k**2 + g_pi**2*conformal_Hubble_parameter(X)**2)*&
             (18.d0*H0**2*X**3*conformal_Hubble_parameter(X)**2*derivative_Omega_DE(X)*(f_pi*wavenumber_k**4*Y(2) + &
             3.d0*e_pi*conformal_Hubble_parameter(X)*(wavenumber_k**2 + g_pi**2*conformal_Hubble_parameter(X)**2)*Y(3) + &
             wavenumber_k**2*(e_pi*(wavenumber_k**2 + g_pi**2*conformal_Hubble_parameter(X)**2)*Y(1) + &
             3.d0*f_pi*conformal_Hubble_parameter(X)*(1.d0 + equation_of_state(X))*Y(4))) + &
             wavenumber_k**4*(wavenumber_k**2 + g_pi**2*conformal_Hubble_parameter(X)**2)*&
             (3.d0*H0**2*X**2*Omega_Matter(X)*Y(1) + 2.d0*(wavenumber_k**2 + 6.d0*conformal_Hubble_parameter(X)**2)*Y(5) - &
             6.d0*conformal_Hubble_parameter(X)**2*Y(6))))/(6.d0*wavenumber_k**4*X*conformal_Hubble_parameter(X)**2*&
             (wavenumber_k**2 + g_pi**2*conformal_Hubble_parameter(X)**2)**2)

     End if

!!$  Else if (MG_parametrisation .eq. 'HS_Basilakos') then
!!$
!!$     F(1) = -(3*H0**2*Omega_m*Y(1) + 3*H0**2*X**3*Omega_DE(X)*Y(2) + 2*X*(wavenumber_k**2*Y(5) + &
!!$          conformal_Hubble_parameter(X)*(Y(3) + 3*conformal_Hubble_parameter(X)*Y(6))))/&
!!$          (2.*X**2*conformal_Hubble_parameter(X)**2)

!!$-(3*H0**2*Omega_m*Y(1) + (-3*H0**2*Omega_m + 3*X*conformal_Hubble_parameter(X)**2)*Y(2) + &
!!$          2*X*(wavenumber_k**2*Y(5) + conformal_Hubble_parameter(X)*(Y(3) + 3*conformal_Hubble_parameter(X)*&
!!$          Y(6))))/(2.*X**2*conformal_Hubble_parameter(X)**2)

!!$     F(2) = -(2*X**3*conformal_Hubble_parameter(X)*Omega_DE(X)*Y(4) + X**2*(1 + &
!!$          equation_of_state(X))*Omega_DE(X)*(3*H0**2*(Omega_m*Y(1) + X**3*Omega_DE(X)*Y(2)) + &
!!$          2*wavenumber_k**2*X*Y(5)) + 6*conformal_Hubble_parameter(X)**2*(Omega_m*dark_energy_pressure_perturbation(X)*Y(1) + &
!!$          X**3*Omega_DE(X)*(Y(6) + equation_of_state(X)*(-Y(2) + Y(6)))))/(2.*X**4*conformal_Hubble_parameter(X)**2*Omega_DE(X))

!!$(2*H0**2*Omega_m*X*conformal_Hubble_parameter(X)*Y(4) - 2*X**2*conformal_Hubble_parameter(X)**3*Y(4) + &
!!$         H0**2*Omega_m*(1 + equation_of_state(X))*(3*H0**2*Omega_m*(Y(1) - Y(2)) + 2*wavenumber_k**2*X*Y(5)) + &
!!$         3*X**2*conformal_Hubble_parameter(X)**4*((-1 + equation_of_state(X))*Y(2) - 2*(1 + equation_of_state(X))*Y(6)) +& 
!!$         X*conformal_Hubble_parameter(X)**2*(-3*H0**2*Omega_m*(1 + 2*dark_energy_pressure_perturbation(X) + &
!!$         equation_of_state(X))*Y(1) + 6*H0**2*Omega_m*Y(2) - 2*(1 + equation_of_state(X))*(wavenumber_k**2*X*Y(5) - &
!!$         3*H0**2*Omega_m*Y(6))))/(2.*X**2*conformal_Hubble_parameter(X)**2*(-(H0**2*Omega_m) + X*conformal_Hubble_parameter(X)**2)) 

!!$     F(3) = -Y(3)/X + wavenumber_k**2*Y(6)/X/conformal_Hubble_parameter(X)
!!$
!!$     F(4) = (Omega_m*wavenumber_k**2*(-2*anisotropic_stress(X,0.d0,0.d0,0.d0,0.d0) + &
!!$          3*dark_energy_pressure_perturbation(X))*Y(1) + &
!!$          3*X**3*Omega_DE(X)*(conformal_Hubble_parameter(X)*(-1 + 3*equation_of_state(X))*Y(4) + wavenumber_k**2*(1 + &
!!$          equation_of_state(X))*Y(6)))/(3.*X**4*conformal_Hubble_parameter(X)*Omega_DE(X))

!!$-((conformal_Hubble_parameter(X)*(1 - 3*equation_of_state(X))*Y(4) + (wavenumber_k**2*&
!!$             ((H0**2*Omega_m*(-2*anisotropic_stress(X,0.d0,0.d0,0.d0,0.d0) + 3*dark_energy_pressure_perturbation(X))*Y(1))/&
!!$             (H0**2*Omega_m - X*conformal_Hubble_parameter(X)**2) - 3*(1 + equation_of_state(X))*Y(6)))/3.)/&
!!$             (X*conformal_Hubble_parameter(X)))

!!$     F(5) = -(3*H0**2*Omega_m*Y(1) + 3*H0**2*X**3*Omega_DE(X)*Y(2) + 2*wavenumber_k**2*X*Y(5) + &
!!$          6*X*conformal_Hubble_parameter(X)**2*Y(6))/(6.*X**2*conformal_Hubble_parameter(X)**2)

!!$-(3*H0**2*Omega_m*Y(1) + (-3*H0**2*Omega_m + 3*X*conformal_Hubble_parameter(X)**2)*Y(2) + &
!!$             2*wavenumber_k**2*X*Y(5) + 6*X*conformal_Hubble_parameter(X)**2*Y(6))/(6.*X**2*conformal_Hubble_parameter(X)**2)

!!$     F(6) = -((3*H0**2*(Omega_m*Y(1) + X**3*Omega_DE(X)*Y(2)) + 2*wavenumber_k**2*X*Y(5))/conformal_Hubble_parameter(X)**2 + &
!!$          (6*X*(3*H0**2*X**3*(anisotropic_stress(X,0.d0,0.d0,0.d0,0.d0)*derivative_Omega_DE(X) + &
!!$          derivative_anisotropic_stress(X)*Omega_DE(X)) + wavenumber_k**2*(2*Y(5) - Y(6))))/wavenumber_k**2)/(6.*X**2)

!!$(-3*H0**2*Omega_m*wavenumber_k**2*Y(1) + 3*wavenumber_k**2*(H0**2*Omega_m - &
!!$             X*conformal_Hubble_parameter(X)**2)*Y(2) - 2*(9*X*conformal_Hubble_parameter(X)**4*(-2*anisotropic_stress(X,&
!!$             0.d0,0.d0,0.d0,0.d0) + X*derivative_anisotropic_stress(X)) + 18*X**2*anisotropic_stress(X,0.d0,0.d0,0.d0,0.d0)*&
!!$             conformal_Hubble_parameter(X)**3*derivative_conformal_Hubble_parameter(X) + wavenumber_k**4*X*Y(5) + &
!!$             3*conformal_Hubble_parameter(X)**2*(3*H0**2*Omega_m*(3*anisotropic_stress(X,0.d0,0.d0,0.d0,0.d0) - &
!!$             X*derivative_anisotropic_stress(X)) + wavenumber_k**2*X*(2*Y(5) - Y(6)))))/(6.*wavenumber_k**2*X**2*&
!!$             conformal_Hubble_parameter(X)**2)

!!$  Else if (MG_parametrisation .eq. 'Starobinsky_Basilakos') then
!!$
!!$     write(UNIT_EXE_FILE,*) 'STAROBINSKY PARAMETRISATION IS NOT IMPLEMENTED YET'
!!$
!!$     stop

  Else if  ( (MG_parametrisation .eq. 'GR_LAMBDA') .and. (approach .eq. 'CHI') ) then 

     F(1) = -(3.d0*H0**2*X**2*Omega_Matter(X)*Y(1) + 2.d0*conformal_Hubble_parameter(X)*Y(2) + &
          2.d0*(wavenumber_k**2 + 3.d0*conformal_Hubble_parameter(X)**2)*Y(3))/&
       (2.d0*X*conformal_Hubble_parameter(X)**2)

     F(2) = -Y(2)/X  + wavenumber_k**2*Y(3)/X/conformal_Hubble_parameter(X)

     F(3) = -(3.d0*H0**2*X**2*Omega_Matter(X)*Y(1) + 2.d0*(wavenumber_k**2 + &
          3*conformal_Hubble_parameter(X)**2)*Y(3))/(6.d0*X*conformal_Hubble_parameter(X)**2)

     F(4) = 0.d0

  Else if ( ( ( MG_parametrisation .eq. 'Savvas' ) .or. (  (MG_parametrisation .eq. 'HS_Basilakos') .or. &
       (MG_parametrisation .eq. 'Starobinsky_Basilakos') ) ) .and. (approach .eq. 'CHI')  ) then

     F(1) = -Y(2)/conformal_Hubble_parameter(X)/X - 3.d0*H0**2*Omega_Matter(X)*&
          Delta_matter(X,wavenumber_k,Y(1),Y(3))/conformal_Hubble_parameter(X)**2/fMG_R_prime(X) - ( 3.d0/X + &
          2.d0*wavenumber_k**2*F_MG(X)/conformal_Hubble_parameter(X)**2/fMG_R_prime(X)/X**2 )*Y(3) + &
          ( 3.d0/2.d0/F_MG(X)/X - 3.d0*(X*derivative_conformal_Hubble_parameter(X) - &
          conformal_Hubble_parameter(X))/conformal_Hubble_parameter(X)/X**2/fMG_R_prime(X))*Y(4)

     F(2) = wavenumber_k**2*Y(3)/conformal_Hubble_parameter(X)/X - Y(2)/X - &
          wavenumber_k**2*Y(4)/conformal_Hubble_parameter(X)/2.d0/F_MG(X)/X

     F(3) = 3.d0*H0**2*Omega_Matter(X)*X*Y(2)/conformal_Hubble_parameter(X)/F_MG(X)/wavenumber_k**2 - &
          ( 1.d0/X + fMG_R_prime(X)/2.d0/F_MG(X) )*Y(3) + 3.d0*fMG_R_prime(X)*Y(4)/4.d0/F_MG(X)**2

     F(4) = -2.d0*H0**2*Omega_Matter(X)*F_MG(X)/conformal_Hubble_parameter(X)**2/fMG_R_prime(X)*&
          Delta_matter(X,wavenumber_k,Y(1),Y(3)) - &
          3.d0*H0**2*Omega_Matter(X)*X*Y(2)/wavenumber_k**2/conformal_Hubble_parameter(X) + &
          ( 1.d0/X - fMG_R_prime(X)/2.d0/F_MG(X) - 2.d0*F_MG(X)/fMG_R_prime(X)/X*(X*&
          derivative_conformal_Hubble_parameter(X) - conformal_Hubble_parameter(X))/&
          conformal_Hubble_parameter(X) )*Y(4) + ( fMG_R_prime(X) - &
          4.d0*F_MG(X)**2*wavenumber_k**2/3.d0/fMG_R_prime(X)/X**2/conformal_Hubble_parameter(X)**2 )*Y(3)

  Else

     write(UNIT_EXE_FILE,*) 'UNKNOWN MG PARAMETRISATION'

     stop

  End if

end subroutine RHSPER                                          

!##########################################################################
! SUBROUTINE THAT COMPUTES RIGHT HAND SIDE THE THE SYSTEM OF EQUATIONS ENDS
!##########################################################################

subroutine JRHSPER(N,X,Y,DFY,LDFY,RPAR,IPAR)

!           use functions
           use fiducial

           IMPLICIT REAL*8 (A-H,O-Z)
           Dimension Y(N),DFY(N,N),RPAR(number_of_parameters)

end subroutine JRHSPER

Subroutine MAS(N,AM,LMAS,RPAR,IPAR,X)

  use background
  use perturbations
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
