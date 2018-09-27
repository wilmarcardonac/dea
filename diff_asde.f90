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
  
  Real*8 :: R0,x0,y0
  Real*8,parameter :: b2 = (-7.d0+sqrt(73.d0))/12.d0
  Real*8,parameter :: a2 = 1.d0 + b2
  Real*8,parameter :: b = 1.5d0 + b2
  Real*8,parameter :: c = 13.d0/6.d0 + 2.d0*b2
  
  R0 = ricci_scalar(final_scale_factor)

  x0 = Lambda/(R0-3.d0*Lambda)

  y0 = R0/(R0-3.d0*Lambda)

  alpha = -fR0/(b2*H0**2/R0*y0**a2*(Lambda/R0)**b2*fgsl_sf_hyperg_2f1(a2,b,c,x0)) 

!  print *, 'Hello', -final_scale_factor*fMG_R_prime(final_scale_factor)/E_H(final_scale_factor)/(1.d0 + fMG_R(final_scale_factor)) 

  print *, fgsl_sf_hyperg_2f1(1.d0,1.d0,1.d0,5.d-1)
 
  stop
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

  RTOL = 1.d-5              ! REQUIRED TOLERANCE
  ATOL = 1.d-1*RTOL!  0*RTOL ! REQUIRED TOLERANCE 
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

  !#############################################
  ! COMPARING ANALYTICAL AND NUMERICAL SOLUTIONS
  !#############################################

  If (MG_parametrisation .eq. 'GR_LAMBDA') then

     call system('cd figures; python plot_numerical_analytical_GR_LAMBDA.py')

  Else if (MG_parametrisation .eq. 'GR_DE') then

     call system('cd figures; python plot_numerical_analytical_GR_DE.py')

  Else if ( ( (MG_parametrisation .eq. 'Savvas') .or. (MG_parametrisation .eq. 'HS_Basilakos') ) .or. &
       (MG_parametrisation .eq. 'Starobinsky_Basilakos') ) then

     If (approach .eq. 'GI') then 

        call system('cd figures; python plot_numerical_analytical_Savvas_GI.py')

     Else

        call system('cd figures; python plot_numerical_analytical_Savvas.py')
        
     End if

  Else

     continue

  End if

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

           WRITE (UNIT_OUTPUT_FILE,98) X,Y(1),Y(2),Y(3),Y(4),Y(5),Y(6)

        Else if (dimension_system_ode .eq. 4) then
           
           WRITE (UNIT_OUTPUT_FILE,100) X,Y(1),Y(2),Y(3)/wavenumber_k,Y(4)/wavenumber_k*(1.d0 + equation_of_state(X))                        

        End if

     Else if ( ( ( (MG_parametrisation .eq. 'Savvas') .or. (MG_parametrisation .eq. 'GR_LAMBDA') ) .or. &
         ( (MG_parametrisation .eq. 'HS_Basilakos') .or. (MG_parametrisation .eq. 'Starobinsky_Basilakos')&
         ) ) .and. (approach .eq. 'CHI') ) then
        
        WRITE (UNIT_OUTPUT_FILE,100) X,Y(1),Y(2),Y(3),Y(4)

     Else if ( ( ( MG_parametrisation .eq. 'Savvas' ) .or. ( (MG_parametrisation .eq. 'HS_Basilakos') &
          .or. (MG_parametrisation .eq. 'Starobinsky_Basilakos') ) ) .and. (approach .eq. 'EF')  ) then

        WRITE (UNIT_OUTPUT_FILE,98) X,Y(1),Y(2),Y(3),Y(4),Y(5),Y(6)
        
     Else if ( ( ( MG_parametrisation .eq. 'Savvas' ) .or. ( (MG_parametrisation .eq. 'HS_Basilakos') &
          .or. (MG_parametrisation .eq. 'Starobinsky_Basilakos') ) ) .and. (approach .eq. 'GI')  ) then

        WRITE (UNIT_OUTPUT_FILE,101) X,Y(1),Y(2),Y(3),Y(4),Y(5)

     End if

     XOUT = initial_scale_factor 

  Else 

10   CONTINUE 

     IF (X.GE.XOUT) THEN 

        ! --- CONTINUOUS OUTPUT FOR RADAU5                                      

        If (MG_parametrisation .eq. 'GR_DE') then

           If (dimension_system_ode .eq. 6) then

              WRITE (UNIT_OUTPUT_FILE,98) XOUT,CONTR5(1,XOUT,CONT,LRC),CONTR5(2,XOUT,CONT,LRC),&    
                   CONTR5(3,XOUT,CONT,LRC),CONTR5(4,XOUT,CONT,LRC),CONTR5(5,XOUT,CONT,LRC),&
                   CONTR5(6,XOUT,CONT,LRC)

           Else if (dimension_system_ode .eq. 4) then

              WRITE (UNIT_OUTPUT_FILE,100) XOUT,CONTR5(1,XOUT,CONT,LRC),CONTR5(2,XOUT,CONT,LRC),&    
                   CONTR5(3,XOUT,CONT,LRC)/wavenumber_k,&
                   CONTR5(4,XOUT,CONT,LRC)/wavenumber_k*(1.d0 + equation_of_state(X))
              
           End if
           
        Else if ( ( ( (MG_parametrisation .eq. 'Savvas') .or. (MG_parametrisation .eq. 'GR_LAMBDA') ) .or. &
         ( (MG_parametrisation .eq. 'HS_Basilakos') .or. (MG_parametrisation .eq. 'Starobinsky_Basilakos')&
         ) ) .and. (approach .eq. 'CHI') ) then

           WRITE (UNIT_OUTPUT_FILE,100) XOUT,CONTR5(1,XOUT,CONT,LRC),CONTR5(2,XOUT,CONT,LRC),&!/wavenumber_k,&    
                CONTR5(3,XOUT,CONT,LRC),CONTR5(4,XOUT,CONT,LRC)

        Else if ( ( ( MG_parametrisation .eq. 'Savvas' ) .or. ( (MG_parametrisation .eq. 'HS_Basilakos') &
          .or. (MG_parametrisation .eq. 'Starobinsky_Basilakos') ) ) .and. (approach .eq. 'EF')  ) then

           WRITE (UNIT_OUTPUT_FILE,98) XOUT,CONTR5(1,XOUT,CONT,LRC),CONTR5(2,XOUT,CONT,LRC),&    
                CONTR5(3,XOUT,CONT,LRC),CONTR5(4,XOUT,CONT,LRC),CONTR5(5,XOUT,CONT,LRC),&
                CONTR5(6,XOUT,CONT,LRC)

        Else if ( ( ( MG_parametrisation .eq. 'Savvas' ) .or. ( (MG_parametrisation .eq. 'HS_Basilakos') &
          .or. (MG_parametrisation .eq. 'Starobinsky_Basilakos') ) ) .and. (approach .eq. 'GI')  ) then

           WRITE (UNIT_OUTPUT_FILE,98) XOUT,CONTR5(1,XOUT,CONT,LRC),CONTR5(2,XOUT,CONT,LRC),&    
                CONTR5(3,XOUT,CONT,LRC),CONTR5(4,XOUT,CONT,LRC),CONTR5(5,XOUT,CONT,LRC)

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

101 FORMAT(E20.10,E20.10,E20.10,E20.10,E20.10,E20.10)
 
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

        !\delta_m
        F(1) = -(3*H0**2*X**2*(Omega_Matter(X)*Y(1) + Omega_DE(X)*Y(2)) + 2*wavenumber_k*(conformal_Hubble_parameter(X)*&
             Y(3) + wavenumber_k*Y(5)) + 6*conformal_Hubble_parameter(X)**2*Y(6))/(2.*X*conformal_Hubble_parameter(X)**2)

        !\delta
        F(2) = (-2*wavenumber_k**2*conformal_Hubble_parameter(X)*Y(4) + 18*conformal_Hubble_parameter(X)**3*&
             (equation_of_state(X) - sound_speed_squared(X))*Y(4) + wavenumber_k*(1 + equation_of_state(X))*&
             (-3*H0**2*X**2*(Omega_Matter(X)*Y(1) + Omega_DE(X)*Y(2)) - 2*wavenumber_k**2*Y(5)) - & 
             6*wavenumber_k*conformal_Hubble_parameter(X)**2*((-equation_of_state(X) + &
             sound_speed_squared(X))*Y(2) + (1 + equation_of_state(X))*Y(6)))/&
             (2.*wavenumber_k*X*conformal_Hubble_parameter(X)**2)

        !V_m
        F(3) = (-(conformal_Hubble_parameter(X)*Y(3)) + wavenumber_k*Y(6))/(X*conformal_Hubble_parameter(X))

        !V
        F(4) = -((conformal_Hubble_parameter(X)*(2*e_pi*Y(3) + (1 + (2*f_pi*wavenumber_k**2)/&
             (wavenumber_k**2 + g_pi**2*conformal_Hubble_parameter(X)**2) - 3*sound_speed_squared(X))*Y(4)) +& 
             (wavenumber_k*(2*e_pi*Y(1) + ((2*f_pi*wavenumber_k**2)/(wavenumber_k**2 + &
             g_pi**2*conformal_Hubble_parameter(X)**2) - 3*sound_speed_squared(X))*Y(2) - &
             3*(1 + equation_of_state(X))*Y(6)))/3.)/(X*conformal_Hubble_parameter(X)))

        !phi
        F(5) = -(3*H0**2*X**2*(Omega_Matter(X)*Y(1) + Omega_DE(X)*Y(2)) + 2*wavenumber_k**2*Y(5) + &
             6*conformal_Hubble_parameter(X)**2*Y(6))/(6.*X*conformal_Hubble_parameter(X)**2)

        !psi
        F(6) = (54*e_pi*g_pi**4*H0**2*X**2*conformal_Hubble_parameter(X)**7*(-(X*derivative_Omega_DE(X)) + &
             Omega_DE(X))*Y(3) + 18*H0**2*wavenumber_k**6*X**2*conformal_Hubble_parameter(X)*Omega_DE(X)*&
             (e_pi*Y(3) + f_pi*Y(4)) - 18*H0**2*wavenumber_k**3*X**2*conformal_Hubble_parameter(X)**3*&
             (3*wavenumber_k*X*derivative_Omega_DE(X)*(e_pi*Y(3) + f_pi*Y(4)) - &
             Omega_DE(X)*(2*f_pi*g_pi**2*X*derivative_conformal_Hubble_parameter(X)*Y(2) + e_pi*(3 + 6*f_pi + &
             2*g_pi**2)*wavenumber_k*Y(3) + f_pi*wavenumber_k*(3 + 6*f_pi + g_pi**2 - &
             9*equation_of_state(X))*Y(4))) - 18*g_pi**2*H0**2*wavenumber_k**2*X**2*&
             conformal_Hubble_parameter(X)**5*(3*X*derivative_Omega_DE(X)*(2*e_pi*Y(3) + f_pi*Y(4)) + &
             Omega_DE(X)*(-(e_pi*(6 + 6*f_pi + g_pi**2)*Y(3)) + 3*f_pi*(-1 + 3*equation_of_state(X))*Y(4))) + &
             wavenumber_k**5*(-wavenumber_k**2 + 9*H0**2*X**2*(e_pi + f_pi + f_pi*equation_of_state(X))*&
             Omega_DE(X))*(3*H0**2*X**2*(Omega_Matter(X)*Y(1) + Omega_DE(X)*Y(2)) + 2*wavenumber_k**2*Y(5)) - &
             g_pi**2*wavenumber_k*conformal_Hubble_parameter(X)**4*(-27*e_pi*g_pi**2*H0**4*X**4*Omega_DE(X)**2*Y(2) - &
             3*H0**2*X**2*Omega_DE(X)*(3*e_pi*(4*f_pi*wavenumber_k**2 + 3*g_pi**2*H0**2*X**2*Omega_Matter(X))*&
             Y(1) - wavenumber_k**2*(g_pi**2 + 18*f_pi*equation_of_state(X))*Y(2) + 6*wavenumber_k*(3*X*&
             derivative_conformal_Hubble_parameter(X)*(-2*e_pi*Y(3) + f_pi*Y(4)) + e_pi*g_pi**2*wavenumber_k*Y(5))) + &
             wavenumber_k**2*(3*g_pi**2*H0**2*X**2*Omega_Matter(X)*Y(1) + 18*H0**2*X**3*derivative_Omega_DE(X)*&
             (2*e_pi*Y(1) + f_pi*Y(2)) + 2*wavenumber_k**2*((12 + g_pi**2)*Y(5) - 6*Y(6)))) + &
             wavenumber_k**3*conformal_Hubble_parameter(X)**2*(27*g_pi**2*H0**4*X**4*(2*e_pi + f_pi + &
             f_pi*equation_of_state(X))*Omega_DE(X)**2*Y(2) + 3*H0**2*X**2*Omega_DE(X)*&
             (2*wavenumber_k**2*(6*f_pi**2 - g_pi**2 - 9*f_pi*equation_of_state(X))*Y(2) + 3*((4*e_pi*f_pi*&
             wavenumber_k**2 + 3*g_pi**2*H0**2*X**2*(2*e_pi + f_pi + f_pi*equation_of_state(X))*&
             Omega_Matter(X))*Y(1) + 2*wavenumber_k*(-3*X*derivative_conformal_Hubble_parameter(X)*&
             (e_pi*Y(3) + f_pi*Y(4)) + g_pi**2*wavenumber_k*(2*e_pi + f_pi + f_pi*equation_of_state(X))*Y(5)&
             ))) - 2*wavenumber_k**2*(3*g_pi**2*H0**2*X**2*Omega_Matter(X)*Y(1) + &
             9*H0**2*X**3*derivative_Omega_DE(X)*(e_pi*Y(1) + f_pi*Y(2)) + wavenumber_k**2*(2*(3 + g_pi**2)*Y(5) &
             - 3*Y(6)))) - 6*g_pi**4*conformal_Hubble_parameter(X)**6*(3*e_pi*H0**2*X**3*(wavenumber_k*&
             derivative_Omega_DE(X)*Y(1) + 3*derivative_conformal_Hubble_parameter(X)*Omega_DE(X)*Y(3)) + &
             2*wavenumber_k**3*Y(5) - wavenumber_k**3*Y(6)))/(6.*wavenumber_k**3*X*conformal_Hubble_parameter(X)**2*&
             (wavenumber_k**2 + g_pi**2*conformal_Hubble_parameter(X)**2)**2)

     End if
     
  Else if  ( (MG_parametrisation .eq. 'GR_LAMBDA') .and. (approach .eq. 'CHI') ) then 

     !delta_m
     F(1) = -(3.d0*H0**2*X**2*Omega_Matter(X)*Y(1) + 2.d0*conformal_Hubble_parameter(X)*Y(2) + &
          2.d0*(wavenumber_k**2 + 3.d0*conformal_Hubble_parameter(X)**2)*Y(3))/&
       (2.d0*X*conformal_Hubble_parameter(X)**2)

     F(2) = -Y(2)/X  + wavenumber_k**2*Y(3)/X/conformal_Hubble_parameter(X)

     !\phi_+
     F(3) = -(3.d0*H0**2*X**2*Omega_Matter(X)*Y(1) + 2.d0*(wavenumber_k**2 + &
          3*conformal_Hubble_parameter(X)**2)*Y(3))/(6.d0*X*conformal_Hubble_parameter(X)**2)

     !\chi
     F(4) = 0.d0

  Else if ( ( ( MG_parametrisation .eq. 'Savvas' ) .or. (  (MG_parametrisation .eq. 'HS_Basilakos') .or. &
       (MG_parametrisation .eq. 'Starobinsky_Basilakos') ) ) .and. (approach .eq. 'CHI')  ) then

     If ( abs(fMG(X)/ricci_scalar(X)) .gt. switch_GR_equations ) then 

        F(1) = -Y(2)/conformal_Hubble_parameter(X)/X - 3.d0*H0**2*Omega_Matter(X)*&
             Delta_matter(X,wavenumber_k,Y(1),Y(3))/conformal_Hubble_parameter(X)**2/fMG_R_prime(X) - ( 3.d0/X + &
             2.d0*wavenumber_k**2*F_MG(X)/conformal_Hubble_parameter(X)**2/fMG_R_prime(X)/X**2 )*Y(3) + &
             ( 3.d0/2.d0/F_MG(X)/X - 3.d0*(X*derivative_conformal_Hubble_parameter(X) - &
             conformal_Hubble_parameter(X))/conformal_Hubble_parameter(X)/X**2/fMG_R_prime(X))*Y(4)

        F(2) = wavenumber_k**2*Y(3)/conformal_Hubble_parameter(X)/X - Y(2)/X - &
             wavenumber_k**2*Y(4)/conformal_Hubble_parameter(X)/2.d0/F_MG(X)/X

        F(3) = 3.d0*H0**2*Omega_Matter(X)*X*Y(2)/2.d0/conformal_Hubble_parameter(X)/F_MG(X)/wavenumber_k**2 - &
             ( 1.d0/X + fMG_R_prime(X)/2.d0/F_MG(X) )*Y(3) + 3.d0*fMG_R_prime(X)*Y(4)/4.d0/F_MG(X)**2

        F(4) = -2.d0*H0**2*Omega_Matter(X)*Delta_matter(X,wavenumber_k,Y(1),Y(3))*F_MG(X)/&
             conformal_Hubble_parameter(X)**2/fMG_R_prime(X) - & 
             3.d0*H0**2*X*Omega_Matter(X)*Y(2)/wavenumber_k**2/conformal_Hubble_parameter(X) + &
             ( 1.d0/X - fMG_R_prime(X)/2.d0/F_MG(X) - 2.d0*F_MG(X)/fMG_R_prime(X)/X**2*(X*&
             derivative_conformal_Hubble_parameter(X) - conformal_Hubble_parameter(X))/&
             conformal_Hubble_parameter(X) )*Y(4) + ( fMG_R_prime(X) - &
             4.d0*F_MG(X)**2*wavenumber_k**2/3.d0/fMG_R_prime(X)/X**2/conformal_Hubble_parameter(X)**2 )*Y(3)

     Else

        F(1) = -(3.d0*H0**2*X**2*Omega_Matter(X)*Y(1) + 2.d0*conformal_Hubble_parameter(X)*Y(2) + &
             2.d0*(wavenumber_k**2 + 3.d0*conformal_Hubble_parameter(X)**2)*Y(3))/&
             (2.d0*X*conformal_Hubble_parameter(X)**2)

        F(2) = -Y(2)/X  + wavenumber_k**2*Y(3)/X/conformal_Hubble_parameter(X)

        F(3) = -(3.d0*H0**2*X**2*Omega_Matter(X)*Y(1) + 2.d0*(wavenumber_k**2 + &
             3*conformal_Hubble_parameter(X)**2)*Y(3))/(6.d0*X*conformal_Hubble_parameter(X)**2)

        F(4) = 0.d0

     End if

  Else if ( ( MG_parametrisation .eq. 'Savvas' ) .and. (approach .eq. 'EF') ) then

     F(1) = -(3*H0**2*X**2*(Omega_Matter(X)*Y(1) + Omega_DE(X)*Y(2)) + 2*wavenumber_k**2*Y(5) + &
          2*conformal_Hubble_parameter(X)*(Y(3) + 3*conformal_Hubble_parameter(X)*Y(6)))/(2.*X*&
          conformal_Hubble_parameter(X)**2)

     F(2) = -(2*conformal_Hubble_parameter(X)*Omega_DE(X)*Y(4) + (1 + equation_of_state(X))*Omega_DE(X)*&
          (3*H0**2*X**2*(Omega_Matter(X)*Y(1) + Omega_DE(X)*Y(2)) + 2*wavenumber_k**2*Y(5)) + &
          6*conformal_Hubble_parameter(X)**2*(dark_energy_pressure_perturbation(X)*Omega_Matter(X)*Y(1) + &
          Omega_DE(X)*(Y(6) + equation_of_state(X)*(-Y(2) + Y(6)))))/(2.*X*conformal_Hubble_parameter(X)**2*Omega_DE(X)) 

     F(3) = -((conformal_Hubble_parameter(X)*Y(3) - wavenumber_k**2*Y(6))/(X*conformal_Hubble_parameter(X)))

     F(4) = (3*conformal_Hubble_parameter(X)*(-1 + 3*equation_of_state(X))*Omega_DE(X)*Y(4) + &
          wavenumber_k**2*((-2*anisotropic_stress(X) + 3*dark_energy_pressure_perturbation(X))*Omega_Matter(X)*Y(1) + &
          3*(1 + equation_of_state(X))*Omega_DE(X)*Y(6)))/(3.*X*conformal_Hubble_parameter(X)*Omega_DE(X))

     F(5) = -(3*H0**2*X**2*Omega_Matter(X)*Y(1) + 3*H0**2*X**2*Omega_DE(X)*Y(2) + 2*wavenumber_k**2*Y(5) + &
          6*conformal_Hubble_parameter(X)**2*Y(6))/(6.*X*conformal_Hubble_parameter(X)**2)

     F(6) = (27*H0**2*X**4*anisotropic_stress(X)*Omega_Matter(X)**2*Y(1) + 3*H0**2*X**2*Omega_DE(X)*&
          (-wavenumber_k**2 + 9*X**2*anisotropic_stress(X)*Omega_Matter(X))*Y(2) - 2*((wavenumber_k**4 + &
          6*wavenumber_k**2*conformal_Hubble_parameter(X)**2)*Y(5) - 3*conformal_Hubble_parameter(X)**2*&
          (-3*X**3*anisotropic_stress(X)*derivative_Omega_Matter(X)*Y(1) + wavenumber_k**2*Y(6))) + &
          3*X**2*Omega_Matter(X)*(-((H0**2*wavenumber_k**2 + 6*X*conformal_Hubble_parameter(X)**2*&
          derivative_anisotropic_stress(X))*Y(1)) + 6*anisotropic_stress(X)*(wavenumber_k**2*Y(5) + &
          conformal_Hubble_parameter(X)*(Y(3) + 3*conformal_Hubble_parameter(X)*Y(6)))))/&
          (6.*wavenumber_k**2*X*conformal_Hubble_parameter(X)**2)

     If (F(1)**2 .ge. 0.d0 ) then

        continue

     Else
        
        print *, 'NaN', F(:) 

        stop

     End if

!!$     ! \delta_m_prime STARTS
!!$     F(1) = -(3*H0**2*X**2*(Omega_Matter(X)*Y(1) + Omega_DE(X)*Y(2)) + &
!!$          2*wavenumber_k*(conformal_Hubble_parameter(X)*Y(3) + wavenumber_k*Y(5)) +& 
!!$          6*conformal_Hubble_parameter(X)**2*Y(6))/&
!!$          (2.*X*conformal_Hubble_parameter(X)**2)
     ! \delta_m prime ENDS

     ! \delta_prime STARTS
!!$     F(2) = (-2*wavenumber_k**2*conformal_Hubble_parameter(X)*Y(4) + &
!!$         18*conformal_Hubble_parameter(X)**3*&
!!$          (equation_of_state(X) - sound_speed_squared(X))*Y(4) + &
!!$         wavenumber_k*(1 + equation_of_state(X))*&
!!$          (-3*H0**2*X**2*(Omega_Matter(X)*Y(1) + Omega_DE(X)*Y(2)) -& 
!!$            2*wavenumber_k**2*Y(5)) - &
!!$         6*wavenumber_k*conformal_Hubble_parameter(X)**2*&
!!$          ((-equation_of_state(X) + sound_speed_squared(X))*Y(2) +& 
!!$            (1 + equation_of_state(X))*Y(6)))/&
!!$       (2.*wavenumber_k*X*conformal_Hubble_parameter(X)**2)
     ! \delta_prime ENDS

     ! V_m_prime STARTS
!!$     F(3) = (-(conformal_Hubble_parameter(X)*Y(3)) + wavenumber_k*Y(6))/(X*conformal_Hubble_parameter(X))
     ! V_m_prime ENDS

     ! V_prime STARTS
!!$     F(4) = ((3*wavenumber_k*X**2*sound_speed_squared(X) + &
!!$            wavenumber_k**3*(-2*f1(X) + 3*f2(X)*sound_speed_squared(X)))*Y(2) +& 
!!$         3*conformal_Hubble_parameter(X)*(X**2 + wavenumber_k**2*f2(X))*&
!!$          (-1 + 3*sound_speed_squared(X))*Y(4) + &
!!$         3*wavenumber_k*(1 + equation_of_state(X))*(X**2 + wavenumber_k**2*f2(X))*Y(6))&
!!$        /(3.*X*conformal_Hubble_parameter(X)*(X**2 + wavenumber_k**2*f2(X)))
     ! V_prime ENDS

     ! \phi_prime STARTS
!!$     F(5) = -(3*H0**2*X**2*(Omega_Matter(X)*Y(1) + Omega_DE(X)*Y(2)) + &
!!$          2*wavenumber_k**2*Y(5) + 6*conformal_Hubble_parameter(X)**2*Y(6))/&
!!$       (6.*X*conformal_Hubble_parameter(X)**2)
     ! \phi_prime ENDS

     ! \psi_prime STARTS
!!$     F(6) = (-(wavenumber_k*(X**2 + wavenumber_k**2*f2(X))*&
!!$            (18*H0**2*X**3*conformal_Hubble_parameter(X)**2*derivative_f1(X)*&
!!$            Omega_DE(X)*Y(2) + (X**2 + wavenumber_k**2*f2(X))*&
!!$            (3*H0**2*X**2*(Omega_Matter(X)*Y(1) + Omega_DE(X)*Y(2)) +& 
!!$            2*(wavenumber_k**2 + 6*conformal_Hubble_parameter(X)**2)*Y(5) -& 
!!$            6*conformal_Hubble_parameter(X)**2*Y(6)))) +& 
!!$            9*H0**2*X**2*f1(X)*(2*wavenumber_k**2*conformal_Hubble_parameter(X)*&
!!$            (X**2 + wavenumber_k**2*f2(X))*Omega_DE(X)*Y(4) + &
!!$            18*conformal_Hubble_parameter(X)**3*(sound_speed_squared(X) - equation_of_state(X))*&
!!$            (X**2 + wavenumber_k**2*f2(X))*Omega_DE(X)*Y(4) + &
!!$            wavenumber_k*(1 + equation_of_state(X))*(X**2 + wavenumber_k**2*f2(X))*&
!!$            Omega_DE(X)*(3*H0**2*X**2*(Omega_Matter(X)*Y(1) + Omega_DE(X)*Y(2)) + &
!!$            2*wavenumber_k**2*Y(5)) + 2*wavenumber_k*conformal_Hubble_parameter(X)**2*&
!!$            (-(X*derivative_Omega_DE(X)*(X**2 + wavenumber_k**2*f2(X))*Y(2)) + &
!!$            Omega_DE(X)*((X**2*(2 + 3*sound_speed_squared(X) - 3*equation_of_state(X)) + &
!!$            wavenumber_k**2*(X*derivative_f2(X) + 3*(sound_speed_squared(X) - &
!!$            equation_of_state(X))*f2(X)))*Y(2) + 3*(1 + equation_of_state(X))*(X**2 + &
!!$            wavenumber_k**2*f2(X))*Y(6)))))/&
!!$            (6.*wavenumber_k*X*conformal_Hubble_parameter(X)**2*&
!!$            (X**2 + wavenumber_k**2*f2(X))**2)
     ! \psi_prime ENDS
!!$
!!$     F(1) = -(3*H0**2*X**2*(Omega_Matter(X)*Y(1) + Omega_DE(X)*Y(2)) + 2*wavenumber_k**2*Y(5) +& 
!!$          2*conformal_Hubble_parameter(X)*(Y(3) + 3*conformal_Hubble_parameter(X)*Y(6)))/(2.*X*&
!!$          conformal_Hubble_parameter(X)**2)
!!$
!!$     F(2) = -(3*H0**2*X**4*(1 + equation_of_state(X))*Omega_DE(X)**3*(H0**2*X**2*Omega_DE(X) + &
!!$             (-X**2 + X**2*F_MG(X) + 2*wavenumber_k**2*fMG_RR(X))*Omega_Matter(X)*Y(1))*(3*H0**2*X**2*&
!!$             Omega_Matter(X)*Y(1) + 3*H0**2*X**2*Omega_DE(X)*Y(2) + 2*wavenumber_k**2*Y(5)) + &
!!$             conformal_Hubble_parameter(X)**2*Omega_DE(X)*(-18*H0**4*X**6*Omega_DE(X)**3*&
!!$             (equation_of_state(X)*(Y(2) - Y(6)) - Y(6)) + 3*H0**2*X**4*Omega_DE(X)**2*Omega_Matter(X)*Y(1)*&
!!$             ((2*(5*H0**2*X**4 - 3*X**2*equation_of_state(X))*F_MG(X) + X**2*(6*equation_of_state(X) + &
!!$             H0**2*X*(3*X**2*fMG_R_prime(X) + 2*(-5*X + (2*wavenumber_k**2 + &
!!$             9*X**4*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR_prime(X)))) + &
!!$             2*fMG_RR(X)*(-6*wavenumber_k**2*equation_of_state(X) + H0**2*X**2*(11*wavenumber_k**2 + &
!!$             114*X**4*derivative_conformal_Hubble_parameter(X)**2 + &
!!$             12*X**5*derivative_conformal_Hubble_parameter(X)*second_derivative_conformal_Hubble_parameter(X))))*Y(2) + &
!!$             6*(1 + equation_of_state(X))*(-X**2 + X**2*F_MG(X) + 2*wavenumber_k**2*fMG_RR(X))*Y(6)) + &
!!$             H0**2*X**2*Omega_DE(X)*Omega_Matter(X)*Y(1)*(2*wavenumber_k**2*((2*X**2*(-1 + 5*X**2)*F_MG(X) + &
!!$             X**2*(3*X**3*fMG_R_prime(X) + 2*(1 - 5*X**2 + (2*wavenumber_k**2*X + &
!!$             9*X**5*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR_prime(X))) + &
!!$             fMG_RR(X)*(2*wavenumber_k**2*(-4 + 11*X**2) + 228*X**6*derivative_conformal_Hubble_parameter(X)**2 + &
!!$             24*X**7*derivative_conformal_Hubble_parameter(X)*second_derivative_conformal_Hubble_parameter(X)))*Y(5) + &
!!$             2*(-X**2 + X**2*F_MG(X) + 2*wavenumber_k**2*fMG_RR(X))*Y(6)) + 3*X**2*Omega_Matter(X)*&
!!$             (Y(1)*(X**3*(3*X**2*fMG_R_prime(X) + 2*(-5*X + (2*wavenumber_k**2 + &
!!$             9*X**4*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR_prime(X)))*(H0**2 - Y(2)) + &
!!$             10*X**4*F_MG(X)**2*Y(2) + 8*wavenumber_k**2*fMG_RR(X)**2*(7*wavenumber_k**2 + &
!!$             72*X**4*derivative_conformal_Hubble_parameter(X)**2 + 6*X**5*derivative_conformal_Hubble_parameter(X)*&
!!$             second_derivative_conformal_Hubble_parameter(X))*Y(2) + 2*X**2*fMG_RR(X)*&
!!$             (H0**2*(11*wavenumber_k**2 + 114*X**4*derivative_conformal_Hubble_parameter(X)**2 + &
!!$             12*X**5*derivative_conformal_Hubble_parameter(X)*second_derivative_conformal_Hubble_parameter(X)) + &
!!$             (-23*wavenumber_k**2 + 4*wavenumber_k**2*X*fMG_R_prime(X) - 6*X**3*derivative_conformal_Hubble_parameter(X)**2*&
!!$             (19*X - 3*wavenumber_k**2*fMG_RR_prime(X)) - 12*X**5*derivative_conformal_Hubble_parameter(X)*&
!!$             second_derivative_conformal_Hubble_parameter(X))*Y(2)) + X**2*F_MG(X)*(10*H0**2*X**2 + &
!!$             (X*(-20*X + 3*X**2*fMG_R_prime(X) + 2*(2*wavenumber_k**2 + &
!!$             9*X**4*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR_prime(X)) + fMG_RR(X)*(46*wavenumber_k**2 + &
!!$             228*X**4*derivative_conformal_Hubble_parameter(X)**2 + 24*X**5*derivative_conformal_Hubble_parameter(X)*&
!!$             second_derivative_conformal_Hubble_parameter(X)))*Y(2))) + 3*(2*X**2*F_MG(X) + &
!!$             6*(wavenumber_k**2 + 10*X**4*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR(X)+ X*(-2*X + &
!!$             X**2*fMG_R_prime(X) - 4*wavenumber_k**2*fMG_RR_prime(X)))*Y(2)*((-X**2 + X**2*F_MG(X) + &
!!$             4*wavenumber_k**2*fMG_RR(X))*Y(5) + (X**2 - X**2*F_MG(X) - 2*wavenumber_k**2*fMG_RR(X))*Y(6)))) + &
!!$             Omega_Matter(X)**2*Y(1)*(3*H0**2*X**4*Omega_Matter(X)*Y(1)*((10*X**4*F_MG(X)**2 - &
!!$             X**3*(3*X**2*fMG_R_prime(X) + 2*(-5*X + (2*wavenumber_k**2 + 9*X**4*derivative_conformal_Hubble_parameter(X)**2)*&
!!$             fMG_RR_prime(X))) + 8*wavenumber_k**2*fMG_RR(X)**2*(7*wavenumber_k**2 + &
!!$             72*X**4*derivative_conformal_Hubble_parameter(X)**2 + 6*X**5*derivative_conformal_Hubble_parameter(X)*&
!!$             second_derivative_conformal_Hubble_parameter(X)) - 2*X**2*fMG_RR(X)*(23*wavenumber_k**2 - &
!!$             4*wavenumber_k**2*X*fMG_R_prime(X) + 6*X**3*derivative_conformal_Hubble_parameter(X)**2*&
!!$             (19*X - 3*wavenumber_k**2*fMG_RR_prime(X)) + 12*X**5*derivative_conformal_Hubble_parameter(X)*&
!!$             second_derivative_conformal_Hubble_parameter(X)) + X**2*F_MG(X)*(X*(-20*X + 3*X**2*fMG_R_prime(X) + &
!!$             2*(2*wavenumber_k**2 + 9*X**4*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR_prime(X)) + &
!!$             fMG_RR(X)*(46*wavenumber_k**2 + 228*X**4*derivative_conformal_Hubble_parameter(X)**2 + &
!!$             24*X**5*derivative_conformal_Hubble_parameter(X)*second_derivative_conformal_Hubble_parameter(X))))*Y(1) + &
!!$             3*(2*X**2*F_MG(X) + 6*(wavenumber_k**2 + 10*X**4*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR(X) + &
!!$             X*(-2*X + X**2*fMG_R_prime(X) - 4*wavenumber_k**2*fMG_RR_prime(X)))*((-X**2 + X**2*F_MG(X) + & 
!!$             4*wavenumber_k**2*fMG_RR(X))*Y(5) + (X**2 - X**2*F_MG(X) - 2*wavenumber_k**2*fMG_RR(X))*Y(6))) + &
!!$             2*wavenumber_k**2*(3*X**2*(-2*X*derivative_conformal_Hubble_parameter(X)*fMG_RR(X)*Y(3) + (2*X**2*F_MG(X) + &
!!$             6*(wavenumber_k**2 + 10*X**4*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR(X) + &
!!$             X*(-2*X + X**2*fMG_R_prime(X) - 4*wavenumber_k**2*fMG_RR_prime(X)))*Y(5))*((-X**2 + X**2*F_MG(X) + &
!!$             4*wavenumber_k**2*fMG_RR(X))*Y(5) + (X**2 - X**2*F_MG(X) - 2*wavenumber_k**2*fMG_RR(X))*Y(6)) + &
!!$             Y(1)*((2*X**4*(-1 + 5*X**2)*F_MG(X)**2 - X**4*(3*X**3*fMG_R_prime(X) + 2*(1 - 5*X**2 + &
!!$             (2*wavenumber_k**2*X + 9*X**5*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR_prime(X))) + &
!!$             8*wavenumber_k**2*fMG_RR(X)**2*(wavenumber_k**2*(-2 + 7*X**2) + &
!!$             72*X**6*derivative_conformal_Hubble_parameter(X)**2 + 6*X**7*derivative_conformal_Hubble_parameter(X)*&
!!$             second_derivative_conformal_Hubble_parameter(X)) + fMG_RR(X)*(2*wavenumber_k**2*X**2*(6 - 23*X**2) + &
!!$             8*wavenumber_k**2*X**5*fMG_R_prime(X) + 12*X**7*derivative_conformal_Hubble_parameter(X)**2*(-19*X + &
!!$             3*wavenumber_k**2*fMG_RR_prime(X)) - 24*X**9*derivative_conformal_Hubble_parameter(X)*&
!!$             second_derivative_conformal_Hubble_parameter(X)) + F_MG(X)*(X**4*(4 - 20*X**2 + 3*X**3*fMG_R_prime(X) + &
!!$             2*(2*wavenumber_k**2*X + 9*X**5*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR_prime(X)) + &
!!$             2*fMG_RR(X)*(wavenumber_k**2*X**2*(-6 + 23*X**2) + 114*X**8*derivative_conformal_Hubble_parameter(X)**2 + &
!!$             12*X**9*derivative_conformal_Hubble_parameter(X)*second_derivative_conformal_Hubble_parameter(X))))*Y(5) + &
!!$             2*(X**2 - X**2*F_MG(X) - 2*wavenumber_k**2*fMG_RR(X))**2*Y(6))))) + &
!!$             2*X**3*conformal_Hubble_parameter(X)*Omega_DE(X)*(3*H0**4*X**3*Omega_DE(X)**3*Y(4) + &
!!$             3*H0**2*X*Omega_DE(X)**2*Omega_Matter(X)*Y(1)*(H0**2*X**3*derivative_conformal_Hubble_parameter(X)*&
!!$             (-X**2 + X**2*F_MG(X) + 3*(wavenumber_k**2 + X**4*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR(X))*Y(2) + &
!!$             (-X**2 + X**2*F_MG(X) + 2*wavenumber_k**2*fMG_RR(X))*Y(4)) + &
!!$             derivative_conformal_Hubble_parameter(X)*Omega_Matter(X)**2*Y(1)*(3*H0**2*X**2*Omega_Matter(X)*Y(1) + &
!!$             2*wavenumber_k**2*Y(5))*((X**4 + X**4*F_MG(X)**2 - (5*wavenumber_k**2*X**2 + &
!!$             3*X**6*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR(X) + (4*wavenumber_k**4 + &
!!$             6*wavenumber_k**2*X**4*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR(X)**2 + &
!!$             F_MG(X)*(-2*X**4 + (5*wavenumber_k**2*X**2 + 3*X**6*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR(X)))*Y(1) + &
!!$             3*wavenumber_k**2*fMG_RR(X)*((X**2 - X**2*F_MG(X) - 4*wavenumber_k**2*fMG_RR(X))*Y(5) + &
!!$             (-X**2 + X**2*F_MG(X) + 2*wavenumber_k**2*fMG_RR(X))*Y(6))) + &
!!$             H0**2*X**2*derivative_conformal_Hubble_parameter(X)*Omega_DE(X)*Omega_Matter(X)*Y(1)*(2*wavenumber_k**2*(-X**2 + &
!!$             X**2*F_MG(X) + 3*(wavenumber_k**2 + X**4*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR(X))*Y(5) + &
!!$             3*Omega_Matter(X)*(Y(1)*(X**4*F_MG(X)**2*Y(2) + 2*wavenumber_k**2*(2*wavenumber_k**2 + &
!!$             3*X**4*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR(X)**2*Y(2) + X**4*(-H0**2 + Y(2)) + &
!!$             X**2*fMG_RR(X)*(3*H0**2*(wavenumber_k**2 + X**4*derivative_conformal_Hubble_parameter(X)**2) - &
!!$             (5*wavenumber_k**2 + 3*X**4*derivative_conformal_Hubble_parameter(X)**2)*Y(2)) + F_MG(X)*(H0**2*X**4 + &
!!$             (-2*X**4 + 5*wavenumber_k**2*X**2*fMG_RR(X) + 3*X**6*derivative_conformal_Hubble_parameter(X)**2*&
!!$             fMG_RR(X))*Y(2))) + 3*wavenumber_k**2*fMG_RR(X)*Y(2)*((X**2 - X**2*F_MG(X) - 4*wavenumber_k**2*fMG_RR(X))*Y(5) + &
!!$             (-X**2 + X**2*F_MG(X) + 2*wavenumber_k**2*fMG_RR(X))*Y(6))))) - 36*X**4*conformal_Hubble_parameter(X)**6*&
!!$             (12*fMG_RR(X) + X*(X*fMG_RR_double_prime(X) + 13*fMG_RR_prime(X)))*Omega_Matter(X)*Y(1)*(-2*H0**2*X**2*&
!!$             Omega_DE(X)**2*Y(5) + X*Omega_Matter(X)*Y(1)*((X**2 - X**2*F_MG(X) - 4*wavenumber_k**2*fMG_RR(X))*Y(5) + &
!!$             (-X**2 + X**2*F_MG(X) + 2*wavenumber_k**2*fMG_RR(X))*Y(6))*derivative_Omega_DE(X) + &
!!$             Omega_DE(X)*(Omega_Matter(X)*(3*Y(6)*((X**2 - X**2*F_MG(X) - 4*wavenumber_k**2*fMG_RR(X))*Y(5) + &
!!$             (-X**2 + X**2*F_MG(X) + 2*wavenumber_k**2*fMG_RR(X))*Y(6)) + Y(1)*((2*X**2 - 2*X**2*F_MG(X) + &
!!$             X**3*fMG_R_prime(X) - 16*wavenumber_k**2*fMG_RR(X) + 4*wavenumber_k**2*X*fMG_RR_prime(X))*Y(5) + &
!!$             (-(X**3*fMG_R_prime(X)) + 2*wavenumber_k**2*fMG_RR(X) - 2*wavenumber_k**2*X*fMG_RR_prime(X))*Y(6))) + &
!!$             X*Y(1)*((-X**2 + X**2*F_MG(X) + 4*wavenumber_k**2*fMG_RR(X))*Y(5) + (X**2 - X**2*F_MG(X) - &
!!$             2*wavenumber_k**2*fMG_RR(X))*Y(6))*derivative_Omega_Matter(X))) + &
!!$             6*X**2*conformal_Hubble_parameter(X)**3*Omega_Matter(X)*Y(1)*(3*H0**4*X**7*Omega_DE(X)**3*&
!!$             (X*(X*derivative_conformal_Hubble_parameter(X)*fMG_RR_double_prime(X) + &
!!$             fMG_RR_prime(X)*(41*derivative_conformal_Hubble_parameter(X) + &
!!$             2*X*second_derivative_conformal_Hubble_parameter(X))) + fMG_RR(X)*(95*derivative_conformal_Hubble_parameter(X) + &
!!$             X*(23*second_derivative_conformal_Hubble_parameter(X) + X*third_derivative_conformal_Hubble_parameter(X))))*Y(2) + &
!!$             H0**2*X**3*Omega_DE(X)**2*(6*wavenumber_k**2*X**2*fMG_RR(X)**2*Omega_Matter(X)*Y(2)*&
!!$             ((135*derivative_conformal_Hubble_parameter(X) + X*(29*second_derivative_conformal_Hubble_parameter(X) + &
!!$             X*third_derivative_conformal_Hubble_parameter(X)))*Y(1) + 6*(20*derivative_conformal_Hubble_parameter(X) + &
!!$             3*X*second_derivative_conformal_Hubble_parameter(X))*(2*Y(5) - Y(6))) + X*(2*wavenumber_k**2*&
!!$             (X**3*derivative_conformal_Hubble_parameter(X)*fMG_RR_double_prime(X) + fMG_RR_prime(X)*((-2 + &
!!$             41*X**2)*derivative_conformal_Hubble_parameter(X) + 2*X**3*second_derivative_conformal_Hubble_parameter(X)))*Y(5) + &
!!$             3*X**4*Omega_Matter(X)*((X*derivative_conformal_Hubble_parameter(X)*fMG_RR_double_prime(X) + &
!!$             fMG_RR_prime(X)*(41*derivative_conformal_Hubble_parameter(X) + &
!!$             2*X*second_derivative_conformal_Hubble_parameter(X)))*Y(1)*(H0**2 + (-1 + F_MG(X))*Y(2)) + &
!!$             27*derivative_conformal_Hubble_parameter(X)*(-1 + F_MG(X))*fMG_RR_prime(X)*Y(2)*Y(5)) + &
!!$             derivative_conformal_Hubble_parameter(X)*(2*X - 2*X**2*fMG_R_prime(X) + 2*wavenumber_k**2*fMG_RR_prime(X) + &
!!$             81*X**4*fMG_RR_prime(X)*Omega_Matter(X)*Y(2) - X*F_MG(X)*(2 + 81*X**3*fMG_RR_prime(X)*Omega_Matter(X)*Y(2)))*&
!!$             Y(6)) + fMG_RR(X)*(2*(wavenumber_k**2*((4 + 95*X**2)*derivative_conformal_Hubble_parameter(X) + &
!!$             X**3*(23*second_derivative_conformal_Hubble_parameter(X) + &
!!$             X*third_derivative_conformal_Hubble_parameter(X)))*Y(5) - derivative_conformal_Hubble_parameter(X)*&
!!$             (2*wavenumber_k**2 + 3*X**4*derivative_conformal_Hubble_parameter(X)**2)*Y(6)) + &
!!$             3*X**3*Omega_Matter(X)*(Y(1)*(derivative_conformal_Hubble_parameter(X)*(95*H0**2*X + (-95*X + 95*X*F_MG(X) + & 
!!$             2*wavenumber_k**2*X*fMG_RR_double_prime(X) + 100*wavenumber_k**2*fMG_RR_prime(X))*Y(2)) + &
!!$             X*(X**2*third_derivative_conformal_Hubble_parameter(X)*(H0**2 + (-1 + F_MG(X))*Y(2)) + &
!!$             second_derivative_conformal_Hubble_parameter(X)*(23*H0**2*X + (-23*X + 23*X*F_MG(X) + &
!!$             4*wavenumber_k**2*fMG_RR_prime(X))*Y(2)))) + 6*Y(2)*((2*derivative_conformal_Hubble_parameter(X)*&
!!$             (-10*X + 10*X*F_MG(X) + 9*wavenumber_k**2*fMG_RR_prime(X)) + 3*X**2*(-1 + F_MG(X))*&
!!$             second_derivative_conformal_Hubble_parameter(X))*Y(5) + (derivative_conformal_Hubble_parameter(X)*&
!!$             (20*X - 20*X*F_MG(X) - 9*wavenumber_k**2*fMG_RR_prime(X)) - 3*X**2*(-1 + &
!!$             F_MG(X))*second_derivative_conformal_Hubble_parameter(X))*Y(6))))) + 2*wavenumber_k**2*X**2*&
!!$             derivative_conformal_Hubble_parameter(X)*fMG_RR(X)*Omega_Matter(X)*Y(1)*((X**2 - X**2*F_MG(X) - &
!!$             4*wavenumber_k**2*fMG_RR(X))*Y(5) + (-X**2 + X**2*F_MG(X) + 2*wavenumber_k**2*fMG_RR(X))*Y(6))*&
!!$             derivative_Omega_DE(X) + Omega_DE(X)*(3*H0**2*X**5*Omega_Matter(X)**2*Y(1)*(X**3*(-1 + &
!!$             F_MG(X))*((X*derivative_conformal_Hubble_parameter(X)*fMG_RR_double_prime(X) + &
!!$             fMG_RR_prime(X)*(41*derivative_conformal_Hubble_parameter(X) + &
!!$             2*X*second_derivative_conformal_Hubble_parameter(X)))*Y(1) + 27*derivative_conformal_Hubble_parameter(X)*&
!!$             fMG_RR_prime(X)*(Y(5) - Y(6))) + 2*wavenumber_k**2*fMG_RR(X)**2*&
!!$             ((135*derivative_conformal_Hubble_parameter(X) + 29*X*second_derivative_conformal_Hubble_parameter(X) + &
!!$             X**2*third_derivative_conformal_Hubble_parameter(X))*Y(1) + 6*(20*derivative_conformal_Hubble_parameter(X) + &
!!$             3*X*second_derivative_conformal_Hubble_parameter(X))*(2*Y(5) - Y(6))) + &
!!$             X*fMG_RR(X)*((derivative_conformal_Hubble_parameter(X)*(-95*X + 95*X*F_MG(X) + &
!!$             2*wavenumber_k**2*X*fMG_RR_double_prime(X) + 100*wavenumber_k**2*fMG_RR_prime(X)) + &
!!$             X*((-23*X + 23*X*F_MG(X) + 4*wavenumber_k**2*fMG_RR_prime(X))*&
!!$             second_derivative_conformal_Hubble_parameter(X) + X**2*(-1 + F_MG(X))*&
!!$             third_derivative_conformal_Hubble_parameter(X)))*Y(1) + 6*((2*derivative_conformal_Hubble_parameter(X)*(-10*X + &
!!$             10*X*F_MG(X) + 9*wavenumber_k**2*fMG_RR_prime(X)) + 3*X**2*(-1 + F_MG(X))*&
!!$             second_derivative_conformal_Hubble_parameter(X))*Y(5) + (derivative_conformal_Hubble_parameter(X)*&
!!$             (20*X - 20*X*F_MG(X) - 9*wavenumber_k**2*fMG_RR_prime(X)) - 3*X**2*(-1 + F_MG(X))*&
!!$             second_derivative_conformal_Hubble_parameter(X))*Y(6)))) + Omega_Matter(X)*((2*X**2*F_MG(X) + &
!!$             6*(wavenumber_k**2 + 10*X**4*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR(X) + X*(-2*X + &
!!$             X**2*fMG_R_prime(X) - 4*wavenumber_k**2*fMG_RR_prime(X)))*Y(3)*((-X**2 + X**2*F_MG(X) + &
!!$             4*wavenumber_k**2*fMG_RR(X))*Y(5) + (X**2 - X**2*F_MG(X) - 2*wavenumber_k**2*fMG_RR(X))*Y(6)) + &
!!$             2*X*(2*fMG_RR(X)**2*(3*wavenumber_k**4*(2*Y(5) - Y(6))*(2*X**2*&
!!$             (20*derivative_conformal_Hubble_parameter(X) + 3*X*second_derivative_conformal_Hubble_parameter(X))*Y(5) - &
!!$             derivative_conformal_Hubble_parameter(X)*Y(6)) + Y(1)*(wavenumber_k**4*((-2 + &
!!$             135*X**2)*derivative_conformal_Hubble_parameter(X) + X**3*(29*second_derivative_conformal_Hubble_parameter(X) + &
!!$             X*third_derivative_conformal_Hubble_parameter(X)))*Y(5) - &
!!$             wavenumber_k**2*derivative_conformal_Hubble_parameter(X)*(wavenumber_k**2 + &
!!$             3*X**4*derivative_conformal_Hubble_parameter(X)**2)*Y(6))) + X**3*(-1 + &
!!$             F_MG(X))*(27*wavenumber_k**2*X**2*derivative_conformal_Hubble_parameter(X)*fMG_RR_prime(X)*Y(5)*(Y(5) - &
!!$             Y(6)) + Y(1)*(wavenumber_k**2*(X**3*derivative_conformal_Hubble_parameter(X)*fMG_RR_double_prime(X) + &
!!$             fMG_RR_prime(X)*((-2 + 41*X**2)*derivative_conformal_Hubble_parameter(X) + &
!!$             2*X**3*second_derivative_conformal_Hubble_parameter(X)))*Y(5) + & 
!!$             derivative_conformal_Hubble_parameter(X)*(X - X*F_MG(X) - X**2*fMG_R_prime(X) + &
!!$             wavenumber_k**2*fMG_RR_prime(X))*Y(6))) + X**2*fMG_RR(X)*(Y(1)*(wavenumber_k**2*&
!!$             (derivative_conformal_Hubble_parameter(X)*(-4 - 95*X**2 + (4 + 95*X**2)*F_MG(X) + X*fMG_R_prime(X) + &
!!$             2*wavenumber_k**2*X**2*fMG_RR_double_prime(X) + 100*wavenumber_k**2*X*fMG_RR_prime(X)) + &
!!$             X**2*((-23*X + 23*X*F_MG(X) + 4*wavenumber_k**2*fMG_RR_prime(X))*second_derivative_conformal_Hubble_parameter(X) + &
!!$             X**2*(-1 + F_MG(X))*third_derivative_conformal_Hubble_parameter(X)))*Y(5) + &
!!$             derivative_conformal_Hubble_parameter(X)*(4*wavenumber_k**2 + &
!!$             3*X**4*derivative_conformal_Hubble_parameter(X)**2 - (4*wavenumber_k**2 + &
!!$             3*X**4*derivative_conformal_Hubble_parameter(X)**2)*F_MG(X) - 3*wavenumber_k**2*X*fMG_R_prime(X))*Y(6)) + &
!!$             3*wavenumber_k**2*(2*X*(2*derivative_conformal_Hubble_parameter(X)*(-10*X + 10*X*F_MG(X) + &
!!$             9*wavenumber_k**2*fMG_RR_prime(X)) + 3*X**2*(-1 + F_MG(X))*second_derivative_conformal_Hubble_parameter(X))*Y(5)**2 + & 
!!$             (-(derivative_conformal_Hubble_parameter(X)*(-1 - 40*X**2 + F_MG(X) + 40*X**2*F_MG(X) + &
!!$             18*wavenumber_k**2*X*fMG_RR_prime(X))) - 6*X**3*(-1 + F_MG(X))*&
!!$             second_derivative_conformal_Hubble_parameter(X))*Y(5)*Y(6) + &
!!$             derivative_conformal_Hubble_parameter(X)*(-1 + F_MG(X))*Y(6)**2)))) + &
!!$             2*wavenumber_k**2*X**2*derivative_conformal_Hubble_parameter(X)*fMG_RR(X)*Y(1)*((-X**2 + X**2*F_MG(X) + &
!!$             4*wavenumber_k**2*fMG_RR(X))*Y(5) + (X**2 - X**2*F_MG(X) - 2*wavenumber_k**2*fMG_RR(X))*Y(6))*&
!!$             derivative_Omega_Matter(X))) + 6*X**2*conformal_Hubble_parameter(X)**4*Omega_Matter(X)*Y(1)*&
!!$             (X**2*(6*H0**4*X**5*(3*X*fMG_RR_double_prime(X) + 19*fMG_RR_prime(X))*Omega_DE(X)**3*Y(2) + &
!!$             H0**2*X*Omega_DE(X)**2*(6*X**4*(3*X*fMG_RR_double_prime(X) + 19*fMG_RR_prime(X))*Omega_Matter(X)*Y(1)*(H0**2 + &
!!$             (-1 + F_MG(X))*Y(2)) + (-4*X + 2*X**2*fMG_R_prime(X) - 4*wavenumber_k**2*X*fMG_RR_double_prime(X) + &
!!$             12*wavenumber_k**2*X**3*fMG_RR_double_prime(X) - 4*wavenumber_k**2*fMG_RR_prime(X) + &
!!$             76*wavenumber_k**2*X**2*fMG_RR_prime(X) - 9*X**5*fMG_RR_double_prime(X)*Omega_Matter(X)*Y(2) - &
!!$             117*X**4*fMG_RR_prime(X)*Omega_Matter(X)*Y(2) + X*F_MG(X)*(4 + 9*X**3*(X*fMG_RR_double_prime(X) + &
!!$             13*fMG_RR_prime(X))*Omega_Matter(X)*Y(2)))*Y(5) + (-2*X**3*fMG_R_double_prime(X) - 5*X**2*fMG_R_prime(X) + &
!!$             2*wavenumber_k**2*X*fMG_RR_double_prime(X) + 9*X**5*fMG_RR_double_prime(X)*Omega_Matter(X)*Y(2) - &
!!$             9*X**5*F_MG(X)*fMG_RR_double_prime(X)*Omega_Matter(X)*Y(2) + fMG_RR_prime(X)*(2*(5*wavenumber_k**2 - &
!!$             9*X**4*derivative_conformal_Hubble_parameter(X)**2) - 117*X**4*(-1 + F_MG(X))*Omega_Matter(X)*Y(2)))*Y(6)) + &
!!$             X**2*(-1 + F_MG(X))*(2*X*F_MG(X) + X**2*fMG_R_prime(X) - 2*(X + 2*wavenumber_k**2*fMG_RR_prime(X)))*&
!!$             Omega_Matter(X)*Y(1)*(Y(5) - Y(6))*derivative_Omega_DE(X) + Omega_DE(X)*(3*H0**2*X**5*(-1 + &
!!$             F_MG(X))*Omega_Matter(X)**2*Y(1)*((6*X*fMG_RR_double_prime(X) + 38*fMG_RR_prime(X))*Y(1) + &
!!$             3*(X*fMG_RR_double_prime(X) + 13*fMG_RR_prime(X))*(Y(5) - Y(6))) + Omega_Matter(X)*(3*X*(-1 + &
!!$             F_MG(X))*(Y(5) - Y(6))*(2*X**2*(9*X*derivative_conformal_Hubble_parameter(X)*fMG_RR_prime(X)*Y(3) + &
!!$             wavenumber_k**2*(X*fMG_RR_double_prime(X) + 13*fMG_RR_prime(X))*Y(5)) + (-2*X + 2*X*F_MG(X) + &
!!$             X**2*fMG_R_prime(X) - 4*wavenumber_k**2*fMG_RR_prime(X))*Y(6)) + Y(1)*((4*X**2*F_MG(X)**2 - &
!!$             X**4*fMG_R_prime(X)**2 + 4*X*F_MG(X)*(X*(-2 + wavenumber_k**2*(-1 + 3*X**2)*fMG_RR_double_prime(X)) + &
!!$             wavenumber_k**2*(-3 + 19*X**2)*fMG_RR_prime(X)) + 4*(X**2*(1 + wavenumber_k**2*(1 - &
!!$             3*X**2)*fMG_RR_double_prime(X)) + wavenumber_k**2*X*(3 - 19*X**2)*fMG_RR_prime(X) + &
!!$             4*wavenumber_k**4*fMG_RR_prime(X)**2))*Y(5) + (X**4*fMG_R_prime(X)**2 + &
!!$             X**2*fMG_R_prime(X)*(3*X - 3*X*F_MG(X) - 2*wavenumber_k**2*fMG_RR_prime(X)) - &
!!$             2*(X**2*(-1 + F_MG(X))*(X**2*fMG_R_double_prime(X) - wavenumber_k**2*fMG_RR_double_prime(X)) + &
!!$             X*(-7*wavenumber_k**2 + 9*X**4*derivative_conformal_Hubble_parameter(X)**2)*(-1 + &
!!$             F_MG(X))*fMG_RR_prime(X) + 4*wavenumber_k**4*fMG_RR_prime(X)**2))*Y(6))) - &
!!$             X**2*(-1 + F_MG(X))*(2*X*F_MG(X) + X**2*fMG_R_prime(X) - 2*(X + &
!!$             2*wavenumber_k**2*fMG_RR_prime(X)))*Y(1)*(Y(5) - Y(6))*derivative_Omega_Matter(X))) + &
!!$             12*wavenumber_k**2*fMG_RR(X)**2*(6*H0**2*X**4*Omega_DE(X)**2*Omega_Matter(X)*Y(2)*(2*Y(1) + &
!!$             6*Y(5) - 3*Y(6)) + X*(wavenumber_k**2 + 10*X**4*derivative_conformal_Hubble_parameter(X)**2)*&
!!$             Omega_Matter(X)*Y(1)*(2*Y(5) - Y(6))*derivative_Omega_DE(X) + Omega_DE(X)*(6*H0**2*X**4*&
!!$             Omega_Matter(X)**2*Y(1)*(2*Y(1) + 6*Y(5) - 3*Y(6)) + Omega_Matter(X)*((2*Y(5) - &
!!$             Y(6))*(2*X**3*(20*derivative_conformal_Hubble_parameter(X) + &
!!$             3*X*second_derivative_conformal_Hubble_parameter(X))*Y(3) + 12*wavenumber_k**2*X**2*Y(5) + &
!!$             3*(wavenumber_k**2 + 10*X**4*derivative_conformal_Hubble_parameter(X)**2)*Y(6)) + &
!!$             Y(1)*(2*(wavenumber_k**2*(5 + 4*X**2) + 40*X**4*derivative_conformal_Hubble_parameter(X)**2)*Y(5) - &
!!$             (3*wavenumber_k**2 + 18*X**4*derivative_conformal_Hubble_parameter(X)**2 + &
!!$             4*X**5*derivative_conformal_Hubble_parameter(X)*second_derivative_conformal_Hubble_parameter(X))*Y(6))) - &
!!$             X*(wavenumber_k**2 + 10*X**4*derivative_conformal_Hubble_parameter(X)**2)*Y(1)*(2*Y(5) - &
!!$             Y(6))*derivative_Omega_Matter(X))) - 2*X*fMG_RR(X)*(-18*H0**4*X**5*Omega_DE(X)**3*Y(2) + &
!!$             3*H0**2*X*Omega_DE(X)**2*(-4*(wavenumber_k**2*(1 + X**2) + &
!!$             5*X**4*derivative_conformal_Hubble_parameter(X)**2)*Y(5) + 2*(wavenumber_k**2 + &
!!$             4*X**4*derivative_conformal_Hubble_parameter(X)**2 + 2*X**5*derivative_conformal_Hubble_parameter(X)*&
!!$             second_derivative_conformal_Hubble_parameter(X))*Y(6) - X**3*Omega_Matter(X)*(Y(1)*(6*H0**2*X + &
!!$             (-6*X + 6*X*F_MG(X) + 7*wavenumber_k**2*X*fMG_RR_double_prime(X) + 51*wavenumber_k**2*fMG_RR_prime(X))*Y(2)) + &
!!$             3*Y(2)*(2*(-3*X + 3*X*F_MG(X) + wavenumber_k**2*X*fMG_RR_double_prime(X) + &
!!$             13*wavenumber_k**2*fMG_RR_prime(X))*Y(5) - (-6*X + 6*X*F_MG(X) + wavenumber_k**2*X*fMG_RR_double_prime(X) + &
!!$             13*wavenumber_k**2*fMG_RR_prime(X))*Y(6)))) + X*Omega_Matter(X)*Y(1)*((7*wavenumber_k**2*X + &
!!$             30*X**5*derivative_conformal_Hubble_parameter(X)**2 - (7*wavenumber_k**2*X + &
!!$             30*X**5*derivative_conformal_Hubble_parameter(X)**2)*F_MG(X) - 2*wavenumber_k**2*X**2*fMG_R_prime(X) + &
!!$             8*wavenumber_k**4*fMG_RR_prime(X))*Y(5) + (-5*wavenumber_k**2*X - &
!!$             30*X**5*derivative_conformal_Hubble_parameter(X)**2 + 5*X*(wavenumber_k**2 + &
!!$             6*X**4*derivative_conformal_Hubble_parameter(X)**2)*F_MG(X) + wavenumber_k**2*X**2*fMG_R_prime(X) - &
!!$             4*wavenumber_k**4*fMG_RR_prime(X))*Y(6))*derivative_Omega_DE(X) + Omega_DE(X)*(-3*H0**2*X**4*&
!!$             Omega_Matter(X)**2*Y(1)*((-6*X + 6*X*F_MG(X) + 7*wavenumber_k**2*X*fMG_RR_double_prime(X) + &
!!$             51*wavenumber_k**2*fMG_RR_prime(X))*Y(1) + 3*(2*(-3*X + 3*X*F_MG(X) + wavenumber_k**2*X*fMG_RR_double_prime(X) + &
!!$             13*wavenumber_k**2*fMG_RR_prime(X))*Y(5) - (-6*X + 6*X*F_MG(X) + wavenumber_k**2*X*fMG_RR_double_prime(X) + &
!!$             13*wavenumber_k**2*fMG_RR_prime(X))*Y(6))) + Omega_Matter(X)*(Y(1)*((-4*(wavenumber_k**2*X*(7 + 3*X**2) + &
!!$             15*X**5*derivative_conformal_Hubble_parameter(X)**2)*F_MG(X) + 5*(-(wavenumber_k**2*X**2) + &
!!$             6*X**6*derivative_conformal_Hubble_parameter(X)**2)*fMG_R_prime(X) + 2*(X*(30*X**4*&
!!$             derivative_conformal_Hubble_parameter(X)**2 + wavenumber_k**2*(14 + 6*X**2 + &
!!$             wavenumber_k**2*(2 - 7*X**2)*fMG_RR_double_prime(X))) + (wavenumber_k**4*(20 - 51*X**2) + &
!!$             60*wavenumber_k**2*X**4*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR_prime(X)))*Y(5) + &
!!$             ((3*wavenumber_k**2*X**2 - 30*X**6*derivative_conformal_Hubble_parameter(X)**2)*fMG_R_prime(X) + &
!!$             4*F_MG(X)*(2*wavenumber_k**2*X + 6*X**5*derivative_conformal_Hubble_parameter(X)**2 + &
!!$             3*X**6*derivative_conformal_Hubble_parameter(X)*second_derivative_conformal_Hubble_parameter(X)) - &
!!$             2*((10*wavenumber_k**4 + 21*wavenumber_k**2*X**4*derivative_conformal_Hubble_parameter(X)**2)*&
!!$             fMG_RR_prime(X) + X*(12*X**4*derivative_conformal_Hubble_parameter(X)**2 + &
!!$             wavenumber_k**2*(4 - X**2*fMG_R_double_prime(X) + wavenumber_k**2*fMG_RR_double_prime(X)) + &
!!$             6*X**5*derivative_conformal_Hubble_parameter(X)*second_derivative_conformal_Hubble_parameter(X))))*Y(6)) - &
!!$             3*(4*wavenumber_k**2*X**2*(-3*X + 3*X*F_MG(X) + wavenumber_k**2*X*fMG_RR_double_prime(X) + &
!!$             13*wavenumber_k**2*fMG_RR_prime(X))*Y(5)**2 - Y(6)*(2*X**3*(derivative_conformal_Hubble_parameter(X)*&
!!$             (-20*X + 20*X*F_MG(X) + 9*wavenumber_k**2*fMG_RR_prime(X)) + 3*X**2*(-1 + F_MG(X))*&
!!$             second_derivative_conformal_Hubble_parameter(X))*Y(3) + (-5*wavenumber_k**2*X - &
!!$             30*X**5*derivative_conformal_Hubble_parameter(X)**2 + 5*X*(wavenumber_k**2 + &
!!$             6*X**4*derivative_conformal_Hubble_parameter(X)**2)*F_MG(X) + wavenumber_k**2*X**2*fMG_R_prime(X) - &
!!$             4*wavenumber_k**4*fMG_RR_prime(X))*Y(6)) + Y(5)*(2*X**3*(2*derivative_conformal_Hubble_parameter(X)*&
!!$             (-10*X + 10*X*F_MG(X) + 9*wavenumber_k**2*fMG_RR_prime(X)) + 3*X**2*(-1 + F_MG(X))*&
!!$             second_derivative_conformal_Hubble_parameter(X))*Y(3) + (-7*wavenumber_k**2*X + 12*wavenumber_k**2*X**3 - &
!!$             30*X**5*derivative_conformal_Hubble_parameter(X)**2 + (wavenumber_k**2*X*(7 - 12*X**2) + &
!!$             30*X**5*derivative_conformal_Hubble_parameter(X)**2)*F_MG(X) + 2*wavenumber_k**2*X**2*fMG_R_prime(X) - &
!!$             2*wavenumber_k**4*X**3*fMG_RR_double_prime(X) - 8*wavenumber_k**4*fMG_RR_prime(X) - &
!!$             26*wavenumber_k**4*X**2*fMG_RR_prime(X))*Y(6)))) + X*Y(1)*((-7*wavenumber_k**2*X - &
!!$             30*X**5*derivative_conformal_Hubble_parameter(X)**2 + (7*wavenumber_k**2*X + &
!!$             30*X**5*derivative_conformal_Hubble_parameter(X)**2)*F_MG(X) + 2*wavenumber_k**2*X**2*fMG_R_prime(X) - &
!!$             8*wavenumber_k**4*fMG_RR_prime(X))*Y(5) + (5*wavenumber_k**2*X + &
!!$             30*X**5*derivative_conformal_Hubble_parameter(X)**2 - 5*X*(wavenumber_k**2 + &
!!$             6*X**4*derivative_conformal_Hubble_parameter(X)**2)*F_MG(X) - wavenumber_k**2*X**2*fMG_R_prime(X) + &
!!$             4*wavenumber_k**4*fMG_RR_prime(X))*Y(6))*derivative_Omega_Matter(X)))) - &
!!$             36*X**4*conformal_Hubble_parameter(X)**5*Omega_Matter(X)*Y(1)*&
!!$             (X**3*(H0**2*X*Omega_DE(X)**2*(-18*derivative_conformal_Hubble_parameter(X)*fMG_RR_prime(X)*Y(5) + &
!!$             (X*derivative_conformal_Hubble_parameter(X)*fMG_RR_double_prime(X) + fMG_RR_prime(X)*&
!!$             (7*derivative_conformal_Hubble_parameter(X) + 2*X*second_derivative_conformal_Hubble_parameter(X)))*Y(6)) - &
!!$             9*X**2*derivative_conformal_Hubble_parameter(X)*(-1 + F_MG(X))*fMG_RR_prime(X)*Omega_Matter(X)*Y(1)*(Y(5) - &
!!$             Y(6))*derivative_Omega_DE(X) + Omega_DE(X)*(Omega_Matter(X)*(-((-1 + F_MG(X))*(X*fMG_RR_double_prime(X) + &
!!$             13*fMG_RR_prime(X))*Y(3)*(Y(5) - Y(6))) - 27*X*derivative_conformal_Hubble_parameter(X)*(-1 + F_MG(X))*&
!!$             fMG_RR_prime(X)*(Y(5) - Y(6))*Y(6) + Y(1)*(9*derivative_conformal_Hubble_parameter(X)*fMG_RR_prime(X)*&
!!$             (2*X - 2*X*F_MG(X) + X**2*fMG_R_prime(X) + 4*wavenumber_k**2*fMG_RR_prime(X))*Y(5) - &
!!$             (-(X**2*derivative_conformal_Hubble_parameter(X)*(-1 + F_MG(X))*fMG_RR_double_prime(X)) + &
!!$             18*wavenumber_k**2*derivative_conformal_Hubble_parameter(X)*fMG_RR_prime(X)**2 + &
!!$             X*fMG_RR_prime(X)*(derivative_conformal_Hubble_parameter(X)*(7 - 7*F_MG(X) + &
!!$             9*X*fMG_R_prime(X)) - 2*X*(-1 + F_MG(X))*second_derivative_conformal_Hubble_parameter(X)))*Y(6))) + &
!!$             9*X**2*derivative_conformal_Hubble_parameter(X)*(-1 + F_MG(X))*fMG_RR_prime(X)*Y(1)*(Y(5) - Y(6))*&
!!$             derivative_Omega_Matter(X))) - 2*wavenumber_k**2*fMG_RR(X)**2*(2*X**2*Omega_Matter(X)*&
!!$             (20*derivative_conformal_Hubble_parameter(X) + 3*X*second_derivative_conformal_Hubble_parameter(X))*Y(1)*&
!!$             (2*Y(5) - Y(6))*derivative_Omega_DE(X) + Omega_DE(X)*(Omega_Matter(X)*(12*Y(3)*(2*Y(5) - Y(6)) - &
!!$             X*(-6*(20*derivative_conformal_Hubble_parameter(X) + 3*X*second_derivative_conformal_Hubble_parameter(X))*&
!!$             (2*Y(5) - Y(6))*Y(6) + Y(1)*(-16*(20*derivative_conformal_Hubble_parameter(X) + &
!!$             3*X*second_derivative_conformal_Hubble_parameter(X))*Y(5) +& 
!!$             (45*derivative_conformal_Hubble_parameter(X) + 13*X*second_derivative_conformal_Hubble_parameter(X) + &
!!$             X**2*third_derivative_conformal_Hubble_parameter(X))*Y(6)))) - &
!!$             2*X**2*(20*derivative_conformal_Hubble_parameter(X) + 3*X*second_derivative_conformal_Hubble_parameter(X))*Y(1)*&
!!$             (2*Y(5) - Y(6))*derivative_Omega_Matter(X))) + X*fMG_RR(X)*(H0**2*X**2*Omega_DE(X)**2*&
!!$             (-4*(20*derivative_conformal_Hubble_parameter(X) + 3*X*second_derivative_conformal_Hubble_parameter(X))*Y(5) + &
!!$             (5*derivative_conformal_Hubble_parameter(X) + X*(7*second_derivative_conformal_Hubble_parameter(X) + &
!!$             X*third_derivative_conformal_Hubble_parameter(X)))*Y(6)) + 2*X**2*Omega_Matter(X)*Y(1)*&
!!$             ((-2*derivative_conformal_Hubble_parameter(X)*(-10*X + 10*X*F_MG(X) + 9*wavenumber_k**2*fMG_RR_prime(X)) - &
!!$             3*X**2*(-1 + F_MG(X))*second_derivative_conformal_Hubble_parameter(X))*Y(5) + &
!!$             (derivative_conformal_Hubble_parameter(X)*(-20*X + 20*X*F_MG(X) + 9*wavenumber_k**2*fMG_RR_prime(X)) + &
!!$             3*X**2*(-1 + F_MG(X))*second_derivative_conformal_Hubble_parameter(X))*Y(6))*&
!!$             derivative_Omega_DE(X) + Omega_DE(X)*(Omega_Matter(X)*(-2*Y(3)*(2*(-3*X + 3*X*F_MG(X) + &
!!$             wavenumber_k**2*X*fMG_RR_double_prime(X) + 13*wavenumber_k**2*fMG_RR_prime(X))*Y(5) - &
!!$             (-6*X + 6*X*F_MG(X) + wavenumber_k**2*X*fMG_RR_double_prime(X) + &
!!$             13*wavenumber_k**2*fMG_RR_prime(X))*Y(6)) + X*(6*Y(6)*((-2*derivative_conformal_Hubble_parameter(X)*&
!!$             (-10*X + 10*X*F_MG(X) + 9*wavenumber_k**2*fMG_RR_prime(X)) - 3*X**2*(-1 + F_MG(X))*&
!!$             second_derivative_conformal_Hubble_parameter(X))*Y(5) + (derivative_conformal_Hubble_parameter(X)*&
!!$             (-20*X + 20*X*F_MG(X) + 9*wavenumber_k**2*fMG_RR_prime(X)) + 3*X**2*(-1 + F_MG(X))*&
!!$             second_derivative_conformal_Hubble_parameter(X))*Y(6)) + Y(1)*&
!!$             (2*(4*derivative_conformal_Hubble_parameter(X)*(10*X - 10*X*F_MG(X) + 5*X**2*fMG_R_prime(X) + &
!!$             2*wavenumber_k**2*fMG_RR_prime(X)) + 3*X*(2*X - 2*X*F_MG(X) + X**2*fMG_R_prime(X) + &
!!$             4*wavenumber_k**2*fMG_RR_prime(X))*second_derivative_conformal_Hubble_parameter(X))*Y(5) + &
!!$             (derivative_conformal_Hubble_parameter(X)*(-5*X + 5*X*F_MG(X) - 40*X**2*fMG_R_prime(X) + &
!!$             2*wavenumber_k**2*X*fMG_RR_double_prime(X) - 48*wavenumber_k**2*fMG_RR_prime(X)) - &
!!$             X*((7*X - 7*X*F_MG(X) + 6*X**2*fMG_R_prime(X) + 8*wavenumber_k**2*fMG_RR_prime(X))*&
!!$             second_derivative_conformal_Hubble_parameter(X) - X**2*(-1 + F_MG(X))*&
!!$             third_derivative_conformal_Hubble_parameter(X)))*Y(6)))) + &
!!$             2*X**2*Y(1)*((2*derivative_conformal_Hubble_parameter(X)*(-10*X + 10*X*F_MG(X) + &
!!$             9*wavenumber_k**2*fMG_RR_prime(X)) + 3*X**2*(-1 + F_MG(X))*&
!!$             second_derivative_conformal_Hubble_parameter(X))*Y(5) + & 
!!$             (derivative_conformal_Hubble_parameter(X)*(20*X - 20*X*F_MG(X) - 9*wavenumber_k**2*fMG_RR_prime(X)) - &
!!$             3*X**2*(-1 + F_MG(X))*second_derivative_conformal_Hubble_parameter(X))*Y(6))*&
!!$             derivative_Omega_Matter(X)))))/(6.*H0**2*X**5*conformal_Hubble_parameter(X)**2*Omega_DE(X)**3*&
!!$             (H0**2*X**2*Omega_DE(X) + (-X**2 + X**2*F_MG(X) + 2*wavenumber_k**2*fMG_RR(X))*Omega_Matter(X)*Y(1)))
!!$
!!$     F(3) = -((conformal_Hubble_parameter(X)*Y(3) - wavenumber_k**2*Y(6))/(X*conformal_Hubble_parameter(X)))
!!$
!!$     F(4) = (2*wavenumber_k**2*X*derivative_conformal_Hubble_parameter(X)*Omega_DE(X)*&
!!$          Omega_Matter(X)*Y(1)*(3*H0**2*X**2*Omega_Matter(X)*Y(1) + 3*H0**2*X**2*Omega_DE(X)*Y(2) + &
!!$          2*wavenumber_k**2*Y(5))*(H0**2*X**2*(-X**2 + X**2*F_MG(X) + 3*(wavenumber_k**2 + &
!!$          X**4*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR(X))*Omega_DE(X) + Omega_Matter(X)*((X**4 + &
!!$          X**4*F_MG(X)**2 - (5*wavenumber_k**2*X**2 + 3*X**6*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR(X) + &
!!$          (4*wavenumber_k**4 + 6*wavenumber_k**2*X**4*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR(X)**2 + &
!!$          F_MG(X)*(-2*X**4 + (5*wavenumber_k**2*X**2 + 3*X**6*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR(X)))*Y(1) + &
!!$          3*wavenumber_k**2*fMG_RR(X)*((X**2 - X**2*F_MG(X) - 4*wavenumber_k**2*fMG_RR(X))*Y(5) + &
!!$          (-X**2 + X**2*F_MG(X) + 2*wavenumber_k**2*fMG_RR(X))*Y(6)))) + wavenumber_k**2*conformal_Hubble_parameter(X)*&
!!$          Omega_DE(X)*(18*H0**4*X**4*(1 + equation_of_state(X))*Omega_DE(X)**3*Y(6) + 3*H0**2*X**2*Omega_DE(X)**2*&
!!$          Omega_Matter(X)*Y(1)*(H0**2*X**2*(10*X**2*F_MG(X) + X*(3*X**2*fMG_R_prime(X) + 2*(-5*X + (2*wavenumber_k**2 + &
!!$          9*X**4*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR_prime(X))) + fMG_RR(X)*(22*wavenumber_k**2 + &
!!$          228*X**4*derivative_conformal_Hubble_parameter(X)**2 + 24*X**5*derivative_conformal_Hubble_parameter(X)*&
!!$          second_derivative_conformal_Hubble_parameter(X)))*Y(2) + 6*(1 + equation_of_state(X))*(-X**2 + X**2*F_MG(X) + &
!!$          2*wavenumber_k**2*fMG_RR(X))*Y(6)) + H0**2*X**2*Omega_DE(X)*Omega_Matter(X)*Y(1)*(2*wavenumber_k**2*(10*X**2*F_MG(X) + &
!!$          X*(3*X**2*fMG_R_prime(X) + 2*(-5*X + (2*wavenumber_k**2 + 9*X**4*derivative_conformal_Hubble_parameter(X)**2)*&
!!$          fMG_RR_prime(X))) + fMG_RR(X)*(22*wavenumber_k**2 + 228*X**4*derivative_conformal_Hubble_parameter(X)**2 + &
!!$          24*X**5*derivative_conformal_Hubble_parameter(X)*second_derivative_conformal_Hubble_parameter(X)))*Y(5) + &
!!$          3*Omega_Matter(X)*(Y(1)*(X**3*(3*X**2*fMG_R_prime(X) + 2*(-5*X + (2*wavenumber_k**2 + &
!!$          9*X**4*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR_prime(X)))*(H0**2 - Y(2)) + &
!!$          10*X**4*F_MG(X)**2*Y(2) + 8*wavenumber_k**2*fMG_RR(X)**2*(7*wavenumber_k**2 + &
!!$          72*X**4*derivative_conformal_Hubble_parameter(X)**2 + 6*X**5*derivative_conformal_Hubble_parameter(X)*&
!!$          second_derivative_conformal_Hubble_parameter(X))*Y(2) + 2*X**2*fMG_RR(X)*(H0**2*(11*wavenumber_k**2 + &
!!$          114*X**4*derivative_conformal_Hubble_parameter(X)**2 + 12*X**5*derivative_conformal_Hubble_parameter(X)*&
!!$          second_derivative_conformal_Hubble_parameter(X)) + (-23*wavenumber_k**2 + 4*wavenumber_k**2*X*fMG_R_prime(X) - &
!!$          6*X**3*derivative_conformal_Hubble_parameter(X)**2*(19*X - 3*wavenumber_k**2*fMG_RR_prime(X)) - &
!!$          12*X**5*derivative_conformal_Hubble_parameter(X)*second_derivative_conformal_Hubble_parameter(X))*Y(2)) + &
!!$          X**2*F_MG(X)*(10*H0**2*X**2 + (X*(-20*X + 3*X**2*fMG_R_prime(X) + 2*(2*wavenumber_k**2 + &
!!$          9*X**4*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR_prime(X)) + &
!!$          fMG_RR(X)*(46*wavenumber_k**2 + 228*X**4*derivative_conformal_Hubble_parameter(X)**2 + &
!!$          24*X**5*derivative_conformal_Hubble_parameter(X)*second_derivative_conformal_Hubble_parameter(X)))*Y(2))) + &
!!$          3*(2*X**2*F_MG(X) + 6*(wavenumber_k**2 + 10*X**4*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR(X) + &
!!$          X*(-2*X + X**2*fMG_R_prime(X) - 4*wavenumber_k**2*fMG_RR_prime(X)))*Y(2)*((-X**2 + X**2*F_MG(X) + &
!!$          4*wavenumber_k**2*fMG_RR(X))*Y(5) + (X**2 - X**2*F_MG(X) - 2*wavenumber_k**2*fMG_RR(X))*Y(6)))) + &
!!$          Omega_Matter(X)**2*Y(1)*(3*H0**2*X**2*Omega_Matter(X)*Y(1)*((10*X**4*F_MG(X)**2 - &
!!$          X**3*(3*X**2*fMG_R_prime(X) + 2*(-5*X + (2*wavenumber_k**2 + 9*X**4*derivative_conformal_Hubble_parameter(X)**2)*&
!!$          fMG_RR_prime(X))) + 8*wavenumber_k**2*fMG_RR(X)**2*(7*wavenumber_k**2 + &
!!$          72*X**4*derivative_conformal_Hubble_parameter(X)**2 + &
!!$          6*X**5*derivative_conformal_Hubble_parameter(X)*second_derivative_conformal_Hubble_parameter(X)) - &
!!$          2*X**2*fMG_RR(X)*(23*wavenumber_k**2 - 4*wavenumber_k**2*X*fMG_R_prime(X) + &
!!$          6*X**3*derivative_conformal_Hubble_parameter(X)**2*(19*X - 3*wavenumber_k**2*fMG_RR_prime(X)) + &
!!$          12*X**5*derivative_conformal_Hubble_parameter(X)*second_derivative_conformal_Hubble_parameter(X)) + &
!!$          X**2*F_MG(X)*(X*(-20*X + 3*X**2*fMG_R_prime(X) + 2*(2*wavenumber_k**2 + &
!!$          9*X**4*derivative_conformal_Hubble_parameter(X)**2)*&
!!$          fMG_RR_prime(X)) + fMG_RR(X)*(46*wavenumber_k**2 + 228*X**4*derivative_conformal_Hubble_parameter(X)**2 + &
!!$          24*X**5*derivative_conformal_Hubble_parameter(X)*second_derivative_conformal_Hubble_parameter(X))))*Y(1) + &
!!$          3*(2*X**2*F_MG(X) + 6*(wavenumber_k**2 + 10*X**4*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR(X) + &
!!$          X*(-2*X + X**2*fMG_R_prime(X) - 4*wavenumber_k**2*fMG_RR_prime(X)))*((-X**2 + X**2*F_MG(X) + &
!!$          4*wavenumber_k**2*fMG_RR(X))*Y(5) + (X**2 - X**2*F_MG(X) - 2*wavenumber_k**2*fMG_RR(X))*Y(6))) + &
!!$          2*wavenumber_k**2*((10*X**4*F_MG(X)**2 - X**3*(3*X**2*fMG_R_prime(X) + 2*(-5*X + (2*wavenumber_k**2 + &
!!$          9*X**4*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR_prime(X))) + 8*wavenumber_k**2*fMG_RR(X)**2*&
!!$          (7*wavenumber_k**2 + 72*X**4*derivative_conformal_Hubble_parameter(X)**2 + &
!!$          6*X**5*derivative_conformal_Hubble_parameter(X)*second_derivative_conformal_Hubble_parameter(X)) - &
!!$          2*X**2*fMG_RR(X)*(23*wavenumber_k**2 - 4*wavenumber_k**2*X*fMG_R_prime(X) + &
!!$          6*X**3*derivative_conformal_Hubble_parameter(X)**2*(19*X - 3*wavenumber_k**2*fMG_RR_prime(X)) + &
!!$          12*X**5*derivative_conformal_Hubble_parameter(X)*second_derivative_conformal_Hubble_parameter(X)) + &
!!$          X**2*F_MG(X)*(X*(-20*X + 3*X**2*fMG_R_prime(X) + 2*(2*wavenumber_k**2 + &
!!$          9*X**4*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR_prime(X)) + fMG_RR(X)*(46*wavenumber_k**2 + &
!!$          228*X**4*derivative_conformal_Hubble_parameter(X)**2 + 24*X**5*derivative_conformal_Hubble_parameter(X)*&
!!$          second_derivative_conformal_Hubble_parameter(X))))*Y(1)*Y(5) + &
!!$          3*(-2*X*derivative_conformal_Hubble_parameter(X)*fMG_RR(X)*Y(3) + (2*X**2*F_MG(X) + &
!!$          6*(wavenumber_k**2 + 10*X**4*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR(X) + &
!!$          X*(-2*X + X**2*fMG_R_prime(X) - 4*wavenumber_k**2*fMG_RR_prime(X)))*Y(5))*((-X**2 + X**2*F_MG(X) + &
!!$          4*wavenumber_k**2*fMG_RR(X))*Y(5) + (X**2 - X**2*F_MG(X) - 2*wavenumber_k**2*fMG_RR(X))*Y(6))))) - &
!!$          36*wavenumber_k**2*X**2*conformal_Hubble_parameter(X)**5*(12*fMG_RR(X) + X*(X*fMG_RR_double_prime(X) + &
!!$          13*fMG_RR_prime(X)))*Omega_Matter(X)*Y(1)*(-2*H0**2*X**2*Omega_DE(X)**2*Y(5) + X*Omega_Matter(X)*Y(1)*&
!!$          ((X**2 - X**2*F_MG(X) - 4*wavenumber_k**2*fMG_RR(X))*Y(5) + (-X**2 + X**2*F_MG(X) + &
!!$          2*wavenumber_k**2*fMG_RR(X))*Y(6))*derivative_Omega_DE(X) + Omega_DE(X)*(Omega_Matter(X)*(3*Y(6)*((X**2 - &
!!$          X**2*F_MG(X) - 4*wavenumber_k**2*fMG_RR(X))*Y(5) + (-X**2 + X**2*F_MG(X) + 2*wavenumber_k**2*fMG_RR(X))*Y(6)) + &
!!$          Y(1)*((2*X**2 - 2*X**2*F_MG(X) + X**3*fMG_R_prime(X) - 16*wavenumber_k**2*fMG_RR(X) + &
!!$          4*wavenumber_k**2*X*fMG_RR_prime(X))*Y(5) + (-(X**3*fMG_R_prime(X)) + 2*wavenumber_k**2*fMG_RR(X) - &
!!$          2*wavenumber_k**2*X*fMG_RR_prime(X))*Y(6))) + X*Y(1)*((-X**2 + X**2*F_MG(X) + 4*wavenumber_k**2*fMG_RR(X))*Y(5) + &
!!$          (X**2 - X**2*F_MG(X) - 2*wavenumber_k**2*fMG_RR(X))*Y(6))*derivative_Omega_Matter(X))) + &
!!$          6*conformal_Hubble_parameter(X)**2*&
!!$          (3*H0**4*X**4*(-1 + 3*equation_of_state(X))*Omega_DE(X)**4*Y(4) + 3*H0**2*X**2*Omega_DE(X)**3*Omega_Matter(X)*Y(1)*&
!!$          (H0**2*wavenumber_k**2*X**5*(X*(X*derivative_conformal_Hubble_parameter(X)*fMG_RR_double_prime(X) + &
!!$          fMG_RR_prime(X)*(41*derivative_conformal_Hubble_parameter(X) + 2*X*second_derivative_conformal_Hubble_parameter(X))) + &
!!$          fMG_RR(X)*(95*derivative_conformal_Hubble_parameter(X) + X*(23*second_derivative_conformal_Hubble_parameter(X) + &
!!$          X*third_derivative_conformal_Hubble_parameter(X))))*Y(2) + (-1 + 3*equation_of_state(X))*(-X**2 + X**2*F_MG(X) + &
!!$          2*wavenumber_k**2*fMG_RR(X))*Y(4)) + H0**2*wavenumber_k**2*X**3*Omega_DE(X)**2*Omega_Matter(X)*Y(1)*&
!!$          (6*wavenumber_k**2*X**2*fMG_RR(X)**2*Omega_Matter(X)*Y(2)*((135*derivative_conformal_Hubble_parameter(X) + &
!!$          X*(29*second_derivative_conformal_Hubble_parameter(X) + X*third_derivative_conformal_Hubble_parameter(X)))*Y(1) + &
!!$          6*(20*derivative_conformal_Hubble_parameter(X) + 3*X*second_derivative_conformal_Hubble_parameter(X))*(2*Y(5) - &
!!$          Y(6))) + X*(2*wavenumber_k**2*(X**3*derivative_conformal_Hubble_parameter(X)*fMG_RR_double_prime(X) + &
!!$          fMG_RR_prime(X)*((-2 + 41*X**2)*derivative_conformal_Hubble_parameter(X) + 2*X**3*&
!!$          second_derivative_conformal_Hubble_parameter(X)))*Y(5) + 3*X**4*Omega_Matter(X)*&
!!$          ((X*derivative_conformal_Hubble_parameter(X)*fMG_RR_double_prime(X) + fMG_RR_prime(X)*&
!!$          (41*derivative_conformal_Hubble_parameter(X) + 2*X*second_derivative_conformal_Hubble_parameter(X)))*Y(1)*&
!!$          (H0**2 + (-1 + F_MG(X))*Y(2)) + 27*derivative_conformal_Hubble_parameter(X)*(-1 + F_MG(X))*&
!!$          fMG_RR_prime(X)*Y(2)*Y(5)) + derivative_conformal_Hubble_parameter(X)*(2*X - 2*X**2*fMG_R_prime(X) + &
!!$          2*wavenumber_k**2*fMG_RR_prime(X) + 81*X**4*fMG_RR_prime(X)*Omega_Matter(X)*Y(2) - &
!!$          X*F_MG(X)*(2 + 81*X**3*fMG_RR_prime(X)*Omega_Matter(X)*Y(2)))*Y(6)) + fMG_RR(X)*(2*&
!!$          (wavenumber_k**2*((4 + 95*X**2)*derivative_conformal_Hubble_parameter(X) + &
!!$          X**3*(23*second_derivative_conformal_Hubble_parameter(X) + &
!!$          X*third_derivative_conformal_Hubble_parameter(X)))*Y(5) - derivative_conformal_Hubble_parameter(X)*&
!!$          (2*wavenumber_k**2 + 3*X**4*derivative_conformal_Hubble_parameter(X)**2)*Y(6)) + &
!!$          3*X**3*Omega_Matter(X)*(Y(1)*(derivative_conformal_Hubble_parameter(X)*(95*H0**2*X + &
!!$          (-95*X + 95*X*F_MG(X) + 2*wavenumber_k**2*X*fMG_RR_double_prime(X) + 100*wavenumber_k**2*fMG_RR_prime(X))*Y(2)) + &
!!$          X*(X**2*third_derivative_conformal_Hubble_parameter(X)*(H0**2 + (-1 + F_MG(X))*Y(2)) + &
!!$          second_derivative_conformal_Hubble_parameter(X)*(23*H0**2*X + (-23*X + 23*X*F_MG(X) + &
!!$          4*wavenumber_k**2*fMG_RR_prime(X))*Y(2)))) + 6*Y(2)*((2*derivative_conformal_Hubble_parameter(X)*&
!!$          (-10*X + 10*X*F_MG(X) + 9*wavenumber_k**2*fMG_RR_prime(X)) + 3*X**2*(-1 + F_MG(X))*&
!!$          second_derivative_conformal_Hubble_parameter(X))*Y(5) + (derivative_conformal_Hubble_parameter(X)*&
!!$          (20*X - 20*X*F_MG(X) - 9*wavenumber_k**2*fMG_RR_prime(X)) - 3*X**2*(-1 + F_MG(X))*&
!!$          second_derivative_conformal_Hubble_parameter(X))*Y(6))))) &
!!$          + 2*wavenumber_k**4*X**2*derivative_conformal_Hubble_parameter(X)*fMG_RR(X)*Omega_Matter(X)**2*Y(1)**2*&
!!$          ((X**2 - X**2*F_MG(X) - 4*wavenumber_k**2*fMG_RR(X))*Y(5) + (-X**2 + X**2*F_MG(X) + &
!!$          2*wavenumber_k**2*fMG_RR(X))*Y(6))*derivative_Omega_DE(X) + wavenumber_k**2*Omega_DE(X)*Omega_Matter(X)*Y(1)*&
!!$          (3*H0**2*X**5*Omega_Matter(X)**2*Y(1)*(X**3*(-1 + F_MG(X))*((X*derivative_conformal_Hubble_parameter(X)*&
!!$          fMG_RR_double_prime(X) + fMG_RR_prime(X)*(41*derivative_conformal_Hubble_parameter(X) + &
!!$          2*X*second_derivative_conformal_Hubble_parameter(X)))*Y(1) + 27*derivative_conformal_Hubble_parameter(X)*fMG_RR_prime(X)*&
!!$          (Y(5) - Y(6))) + 2*wavenumber_k**2*fMG_RR(X)**2*((135*derivative_conformal_Hubble_parameter(X) + &
!!$          29*X*second_derivative_conformal_Hubble_parameter(X) + X**2*third_derivative_conformal_Hubble_parameter(X))*Y(1) + &
!!$          6*(20*derivative_conformal_Hubble_parameter(X) + 3*X*second_derivative_conformal_Hubble_parameter(X))*(2*Y(5) - &
!!$          Y(6))) + X*fMG_RR(X)*((derivative_conformal_Hubble_parameter(X)*(-95*X + 95*X*F_MG(X) + &
!!$          2*wavenumber_k**2*X*fMG_RR_double_prime(X) + 100*wavenumber_k**2*fMG_RR_prime(X)) + &
!!$          X*((-23*X + 23*X*F_MG(X) + 4*wavenumber_k**2*fMG_RR_prime(X))*second_derivative_conformal_Hubble_parameter(X) + &
!!$          X**2*(-1 + F_MG(X))*third_derivative_conformal_Hubble_parameter(X)))*Y(1) + &
!!$          6*((2*derivative_conformal_Hubble_parameter(X)*(-10*X + 10*X*F_MG(X) + 9*wavenumber_k**2*fMG_RR_prime(X)) + &
!!$          3*X**2*(-1 + F_MG(X))*second_derivative_conformal_Hubble_parameter(X))*Y(5) + &
!!$          (derivative_conformal_Hubble_parameter(X)*(20*X - 20*X*F_MG(X) - 9*wavenumber_k**2*fMG_RR_prime(X)) - &
!!$          3*X**2*(-1 + F_MG(X))*second_derivative_conformal_Hubble_parameter(X))*Y(6)))) + Omega_Matter(X)*((2*X**2*F_MG(X) + &
!!$          6*(wavenumber_k**2 + 10*X**4*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR(X) + &
!!$          X*(-2*X + X**2*fMG_R_prime(X) - 4*wavenumber_k**2*fMG_RR_prime(X)))*Y(3)*((-X**2 + X**2*F_MG(X) + &
!!$          4*wavenumber_k**2*fMG_RR(X))*Y(5) + (X**2 - X**2*F_MG(X) - 2*wavenumber_k**2*fMG_RR(X))*Y(6)) + &
!!$          2*X*(2*fMG_RR(X)**2*(3*wavenumber_k**4*(2*Y(5) - Y(6))*(2*X**2*(20*derivative_conformal_Hubble_parameter(X) + &
!!$          3*X*second_derivative_conformal_Hubble_parameter(X))*Y(5) - derivative_conformal_Hubble_parameter(X)*Y(6)) + &
!!$          Y(1)*(wavenumber_k**4*((-2 + 135*X**2)*derivative_conformal_Hubble_parameter(X) + &
!!$          X**3*(29*second_derivative_conformal_Hubble_parameter(X) + X*third_derivative_conformal_Hubble_parameter(X)))*Y(5) - & 
!!$          wavenumber_k**2*derivative_conformal_Hubble_parameter(X)*(wavenumber_k**2 + &
!!$          3*X**4*derivative_conformal_Hubble_parameter(X)**2)*Y(6))) + X**3*(-1 + F_MG(X))*&
!!$          (27*wavenumber_k**2*X**2*derivative_conformal_Hubble_parameter(X)*fMG_RR_prime(X)*Y(5)*&
!!$          (Y(5) - Y(6)) + Y(1)*(wavenumber_k**2*(X**3*derivative_conformal_Hubble_parameter(X)*fMG_RR_double_prime(X) + &
!!$          fMG_RR_prime(X)*((-2 + 41*X**2)*derivative_conformal_Hubble_parameter(X) + &
!!$          2*X**3*second_derivative_conformal_Hubble_parameter(X)))*Y(5) + derivative_conformal_Hubble_parameter(X)*&
!!$          (X - X*F_MG(X) - X**2*fMG_R_prime(X) + wavenumber_k**2*fMG_RR_prime(X))*Y(6))) + &
!!$          X**2*fMG_RR(X)*(Y(1)*(wavenumber_k**2*(derivative_conformal_Hubble_parameter(X)*(-4 - 95*X**2 + (4 + &
!!$          95*X**2)*F_MG(X) + X*fMG_R_prime(X) + 2*wavenumber_k**2*X**2*fMG_RR_double_prime(X) + &
!!$          100*wavenumber_k**2*X*fMG_RR_prime(X)) + X**2*((-23*X + 23*X*F_MG(X) + &
!!$          4*wavenumber_k**2*fMG_RR_prime(X))*second_derivative_conformal_Hubble_parameter(X) + X**2*(-1 + &
!!$          F_MG(X))*third_derivative_conformal_Hubble_parameter(X)))*Y(5) + derivative_conformal_Hubble_parameter(X)*&
!!$          (4*wavenumber_k**2 + 3*X**4*derivative_conformal_Hubble_parameter(X)**2 - (4*wavenumber_k**2 + &
!!$          3*X**4*derivative_conformal_Hubble_parameter(X)**2)*F_MG(X) - 3*wavenumber_k**2*X*fMG_R_prime(X))*Y(6)) + &
!!$          3*wavenumber_k**2*(2*X*(2*derivative_conformal_Hubble_parameter(X)*(-10*X + 10*X*F_MG(X) + &
!!$          9*wavenumber_k**2*fMG_RR_prime(X)) + 3*X**2*(-1 + F_MG(X))*second_derivative_conformal_Hubble_parameter(X))*Y(5)**2&
!!$          + (-(derivative_conformal_Hubble_parameter(X)*(-1 - 40*X**2 + F_MG(X) + 40*X**2*F_MG(X) + &
!!$          18*wavenumber_k**2*X*fMG_RR_prime(X))) - 6*X**3*(-1 + &
!!$          F_MG(X))*second_derivative_conformal_Hubble_parameter(X))*Y(5)*Y(6) + & 
!!$          derivative_conformal_Hubble_parameter(X)*(-1 + F_MG(X))*Y(6)**2)))) + 2*wavenumber_k**2*X**2*&
!!$          derivative_conformal_Hubble_parameter(X)*fMG_RR(X)*Y(1)*((-X**2 + X**2*F_MG(X) + 4*wavenumber_k**2*fMG_RR(X))*&
!!$          Y(5) + (X**2 - X**2*F_MG(X) - 2*wavenumber_k**2*fMG_RR(X))*Y(6))*derivative_Omega_Matter(X))) + &
!!$          6*wavenumber_k**2*conformal_Hubble_parameter(X)**3*Omega_Matter(X)*Y(1)*(X**2*(6*H0**4*X**5*&
!!$          (3*X*fMG_RR_double_prime(X) + 19*fMG_RR_prime(X))*Omega_DE(X)**3*Y(2) + &
!!$          H0**2*X*Omega_DE(X)**2*(6*X**4*(3*X*fMG_RR_double_prime(X) + 19*fMG_RR_prime(X))*&
!!$          Omega_Matter(X)*Y(1)*(H0**2 + (-1 + F_MG(X))*Y(2)) + (-4*X + 2*X**2*fMG_R_prime(X) - &
!!$          4*wavenumber_k**2*X*fMG_RR_double_prime(X) + 12*wavenumber_k**2*X**3*fMG_RR_double_prime(X) - &
!!$          4*wavenumber_k**2*fMG_RR_prime(X) + &
!!$          76*wavenumber_k**2*X**2*fMG_RR_prime(X) - 9*X**5*fMG_RR_double_prime(X)*Omega_Matter(X)*Y(2) - &
!!$          117*X**4*fMG_RR_prime(X)*Omega_Matter(X)*Y(2) + X*F_MG(X)*(4 + 9*X**3*(X*fMG_RR_double_prime(X) + &
!!$          13*fMG_RR_prime(X))*Omega_Matter(X)*Y(2)))*Y(5) + (-2*X**3*fMG_R_double_prime(X) - 5*X**2*fMG_R_prime(X) + &
!!$          2*wavenumber_k**2*X*fMG_RR_double_prime(X) + 9*X**5*fMG_RR_double_prime(X)*Omega_Matter(X)*Y(2) - &
!!$          9*X**5*F_MG(X)*fMG_RR_double_prime(X)*Omega_Matter(X)*Y(2) + fMG_RR_prime(X)*(2*(5*wavenumber_k**2 - &
!!$          9*X**4*derivative_conformal_Hubble_parameter(X)**2) - 117*X**4*(-1 + F_MG(X))*Omega_Matter(X)*Y(2)))*Y(6)) + &
!!$          X**2*(-1 + F_MG(X))*(2*X*F_MG(X) + X**2*fMG_R_prime(X) - 2*(X + 2*wavenumber_k**2*fMG_RR_prime(X)))*Omega_Matter(X)*Y(1)*&
!!$          (Y(5) - Y(6))*derivative_Omega_DE(X) + Omega_DE(X)*(3*H0**2*X**5*(-1 + F_MG(X))*Omega_Matter(X)**2*Y(1)*&
!!$          ((6*X*fMG_RR_double_prime(X) + 38*fMG_RR_prime(X))*Y(1) + 3*(X*fMG_RR_double_prime(X) + &
!!$          13*fMG_RR_prime(X))*(Y(5) - Y(6))) + & 
!!$          Omega_Matter(X)*(3*X*(-1 + F_MG(X))*(Y(5) - Y(6))*(2*X**2*(9*X*derivative_conformal_Hubble_parameter(X)*&
!!$          fMG_RR_prime(X)*Y(3) + wavenumber_k**2*(X*fMG_RR_double_prime(X) + 13*fMG_RR_prime(X))*Y(5)) + &
!!$          (-2*X + 2*X*F_MG(X) + X**2*fMG_R_prime(X) - 4*wavenumber_k**2*fMG_RR_prime(X))*Y(6)) + Y(1)*((4*X**2*F_MG(X)**2 - &
!!$          X**4*fMG_R_prime(X)**2 + 4*X*F_MG(X)*(X*(-2 + wavenumber_k**2*(-1 + 3*X**2)*fMG_RR_double_prime(X)) + &
!!$          wavenumber_k**2*(-3 + 19*X**2)*fMG_RR_prime(X)) + 4*(X**2*(1 + wavenumber_k**2*(1 - 3*X**2)*fMG_RR_double_prime(X)) + &
!!$          wavenumber_k**2*X*(3 - 19*X**2)*fMG_RR_prime(X) + 4*wavenumber_k**4*fMG_RR_prime(X)**2))*Y(5) + &
!!$          (X**4*fMG_R_prime(X)**2 + X**2*fMG_R_prime(X)*(3*X - 3*X*F_MG(X) - 2*wavenumber_k**2*fMG_RR_prime(X)) - &
!!$          2*(X**2*(-1 + F_MG(X))*(X**2*fMG_R_double_prime(X) - wavenumber_k**2*fMG_RR_double_prime(X)) + &
!!$          X*(-7*wavenumber_k**2 + 9*X**4*derivative_conformal_Hubble_parameter(X)**2)*(-1 + F_MG(X))*fMG_RR_prime(X) + &
!!$          4*wavenumber_k**4*fMG_RR_prime(X)**2))*Y(6))) - X**2*(-1 + F_MG(X))*(2*X*F_MG(X) + X**2*fMG_R_prime(X) - &
!!$          2*(X + 2*wavenumber_k**2*fMG_RR_prime(X)))*Y(1)*(Y(5) - Y(6))*derivative_Omega_Matter(X))) + &
!!$          12*wavenumber_k**2*fMG_RR(X)**2*(6*H0**2*X**4*Omega_DE(X)**2*Omega_Matter(X)*Y(2)*(2*Y(1) + 6*Y(5) - 3*Y(6)) + &
!!$          X*(wavenumber_k**2 + 10*X**4*derivative_conformal_Hubble_parameter(X)**2)*Omega_Matter(X)*Y(1)*(2*Y(5) - &
!!$          Y(6))*derivative_Omega_DE(X) + Omega_DE(X)*(6*H0**2*X**4*Omega_Matter(X)**2*Y(1)*(2*Y(1) + 6*Y(5) - 3*Y(6)) + &
!!$          Omega_Matter(X)*((2*Y(5) - Y(6))*(2*X**3*(20*derivative_conformal_Hubble_parameter(X) + &
!!$          3*X*second_derivative_conformal_Hubble_parameter(X))*Y(3) + 12*wavenumber_k**2*X**2*Y(5) + &
!!$          3*(wavenumber_k**2 + 10*X**4*derivative_conformal_Hubble_parameter(X)**2)*Y(6)) + &
!!$          Y(1)*(2*(wavenumber_k**2*(5 + 4*X**2) + 40*X**4*derivative_conformal_Hubble_parameter(X)**2)*Y(5) - &
!!$          (3*wavenumber_k**2 + 18*X**4*derivative_conformal_Hubble_parameter(X)**2 + &
!!$          4*X**5*derivative_conformal_Hubble_parameter(X)*second_derivative_conformal_Hubble_parameter(X))*Y(6))) - &
!!$          X*(wavenumber_k**2 + 10*X**4*derivative_conformal_Hubble_parameter(X)**2)*Y(1)*(2*Y(5) - &
!!$          Y(6))*derivative_Omega_Matter(X))) - 2*X*fMG_RR(X)*(-18*H0**4*X**5*Omega_DE(X)**3*Y(2) + &
!!$          3*H0**2*X*Omega_DE(X)**2*(-4*(wavenumber_k**2*(1 + X**2) + 5*X**4*derivative_conformal_Hubble_parameter(X)**2)*Y(5) + &
!!$          2*(wavenumber_k**2 + 4*X**4*derivative_conformal_Hubble_parameter(X)**2 + &
!!$          2*X**5*derivative_conformal_Hubble_parameter(X)*second_derivative_conformal_Hubble_parameter(X))*Y(6) - &
!!$          X**3*Omega_Matter(X)*(Y(1)*(6*H0**2*X + (-6*X + 6*X*F_MG(X) + 7*wavenumber_k**2*X*fMG_RR_double_prime(X) + &
!!$          51*wavenumber_k**2*fMG_RR_prime(X))*Y(2)) + 3*Y(2)*(2*(-3*X + 3*X*F_MG(X) + wavenumber_k**2*X*fMG_RR_double_prime(X) + &
!!$          13*wavenumber_k**2*fMG_RR_prime(X))*Y(5) - (-6*X + 6*X*F_MG(X) + wavenumber_k**2*X*fMG_RR_double_prime(X) + &
!!$          13*wavenumber_k**2*fMG_RR_prime(X))*Y(6)))) + X*Omega_Matter(X)*Y(1)*((7*wavenumber_k**2*X + &
!!$          30*X**5*derivative_conformal_Hubble_parameter(X)**2 - (7*wavenumber_k**2*X + &
!!$          30*X**5*derivative_conformal_Hubble_parameter(X)**2)*F_MG(X) - 2*wavenumber_k**2*X**2*fMG_R_prime(X) + &
!!$          8*wavenumber_k**4*fMG_RR_prime(X))*Y(5) + (-5*wavenumber_k**2*X - &
!!$          30*X**5*derivative_conformal_Hubble_parameter(X)**2 + 5*X*(wavenumber_k**2 + &
!!$          6*X**4*derivative_conformal_Hubble_parameter(X)**2)*F_MG(X) + wavenumber_k**2*X**2*fMG_R_prime(X) - &
!!$          4*wavenumber_k**4*fMG_RR_prime(X))*Y(6))*derivative_Omega_DE(X) + Omega_DE(X)*(-3*H0**2*X**4*Omega_Matter(X)**2*Y(1)*&
!!$          ((-6*X + 6*X*F_MG(X) + 7*wavenumber_k**2*X*fMG_RR_double_prime(X) + 51*wavenumber_k**2*fMG_RR_prime(X))*Y(1) + &
!!$          3*(2*(-3*X + 3*X*F_MG(X) + wavenumber_k**2*X*fMG_RR_double_prime(X) + 13*wavenumber_k**2*fMG_RR_prime(X))*Y(5) - &
!!$          (-6*X + 6*X*F_MG(X) + wavenumber_k**2*X*fMG_RR_double_prime(X) + 13*wavenumber_k**2*fMG_RR_prime(X))*Y(6))) + &
!!$          Omega_Matter(X)*(Y(1)*((-4*(wavenumber_k**2*X*(7 + 3*X**2) + &
!!$          15*X**5*derivative_conformal_Hubble_parameter(X)**2)*F_MG(X) + 5*(-(wavenumber_k**2*X**2) + &
!!$          6*X**6*derivative_conformal_Hubble_parameter(X)**2)*fMG_R_prime(X) + 2*(X*(30*X**4*&
!!$          derivative_conformal_Hubble_parameter(X)**2 + wavenumber_k**2*(14 + 6*X**2 + wavenumber_k**2*(2 - &
!!$          7*X**2)*fMG_RR_double_prime(X))) + (wavenumber_k**4*(20 - 51*X**2) + 60*wavenumber_k**2*X**4*&
!!$          derivative_conformal_Hubble_parameter(X)**2)*fMG_RR_prime(X)))*Y(5) + ((3*wavenumber_k**2*X**2 - &
!!$          30*X**6*derivative_conformal_Hubble_parameter(X)**2)*fMG_R_prime(X) + 4*F_MG(X)*(2*wavenumber_k**2*X + &
!!$          6*X**5*derivative_conformal_Hubble_parameter(X)**2 + 3*X**6*derivative_conformal_Hubble_parameter(X)*&
!!$          second_derivative_conformal_Hubble_parameter(X)) - 2*((10*wavenumber_k**4 + &
!!$          21*wavenumber_k**2*X**4*derivative_conformal_Hubble_parameter(X)**2)*fMG_RR_prime(X) + &
!!$          X*(12*X**4*derivative_conformal_Hubble_parameter(X)**2 + wavenumber_k**2*(4 - X**2*fMG_R_double_prime(X) + &
!!$          wavenumber_k**2*fMG_RR_double_prime(X)) + 6*X**5*derivative_conformal_Hubble_parameter(X)*&
!!$          second_derivative_conformal_Hubble_parameter(X))))*Y(6)) - 3*(4*wavenumber_k**2*X**2*(-3*X + 3*X*F_MG(X) + & 
!!$          wavenumber_k**2*X*fMG_RR_double_prime(X) + 13*wavenumber_k**2*fMG_RR_prime(X))*Y(5)**2 - &
!!$          Y(6)*(2*X**3*(derivative_conformal_Hubble_parameter(X)*(-20*X + 20*X*F_MG(X) + &
!!$          9*wavenumber_k**2*fMG_RR_prime(X)) + 3*X**2*(-1 + F_MG(X))*second_derivative_conformal_Hubble_parameter(X))*Y(3) + &
!!$          (-5*wavenumber_k**2*X - 30*X**5*derivative_conformal_Hubble_parameter(X)**2 + 5*X*(wavenumber_k**2 + &
!!$          6*X**4*derivative_conformal_Hubble_parameter(X)**2)*F_MG(X) + wavenumber_k**2*X**2*fMG_R_prime(X) - &
!!$          4*wavenumber_k**4*fMG_RR_prime(X))*Y(6)) + Y(5)*(2*X**3*(2*derivative_conformal_Hubble_parameter(X)*&
!!$          (-10*X + 10*X*F_MG(X) + 9*wavenumber_k**2*fMG_RR_prime(X)) + 3*X**2*(-1 + F_MG(X))*&
!!$          second_derivative_conformal_Hubble_parameter(X))*Y(3) + (-7*wavenumber_k**2*X + 12*wavenumber_k**2*X**3 - &
!!$          30*X**5*derivative_conformal_Hubble_parameter(X)**2 + (wavenumber_k**2*X*(7 - 12*X**2) + &
!!$          30*X**5*derivative_conformal_Hubble_parameter(X)**2)*F_MG(X) + 2*wavenumber_k**2*X**2*fMG_R_prime(X) - &
!!$          2*wavenumber_k**4*X**3*fMG_RR_double_prime(X) - 8*wavenumber_k**4*fMG_RR_prime(X) - &
!!$          26*wavenumber_k**4*X**2*fMG_RR_prime(X))*Y(6)))) + X*Y(1)*((-7*wavenumber_k**2*X - &
!!$          30*X**5*derivative_conformal_Hubble_parameter(X)**2 + (7*wavenumber_k**2*X + &
!!$          30*X**5*derivative_conformal_Hubble_parameter(X)**2)*F_MG(X) + 2*wavenumber_k**2*X**2*fMG_R_prime(X) - &
!!$          8*wavenumber_k**4*fMG_RR_prime(X))*Y(5) + (5*wavenumber_k**2*X + &
!!$          30*X**5*derivative_conformal_Hubble_parameter(X)**2 - 5*X*(wavenumber_k**2 + &
!!$          6*X**4*derivative_conformal_Hubble_parameter(X)**2)*F_MG(X) - wavenumber_k**2*X**2*fMG_R_prime(X) + &
!!$          4*wavenumber_k**4*fMG_RR_prime(X))*Y(6))*derivative_Omega_Matter(X)))) - 36*wavenumber_k**2*X**2*&
!!$          conformal_Hubble_parameter(X)**4*Omega_Matter(X)*Y(1)*(X**3*(H0**2*X*Omega_DE(X)**2*&
!!$          (-18*derivative_conformal_Hubble_parameter(X)*fMG_RR_prime(X)*Y(5) + &
!!$          (X*derivative_conformal_Hubble_parameter(X)*fMG_RR_double_prime(X) + fMG_RR_prime(X)*&
!!$          (7*derivative_conformal_Hubble_parameter(X) + 2*X*second_derivative_conformal_Hubble_parameter(X)))*Y(6)) - &
!!$          9*X**2*derivative_conformal_Hubble_parameter(X)*(-1 + F_MG(X))*fMG_RR_prime(X)*Omega_Matter(X)*Y(1)*(Y(5) - &
!!$          Y(6))*derivative_Omega_DE(X) + Omega_DE(X)*(Omega_Matter(X)*(-((-1 + F_MG(X))*(X*fMG_RR_double_prime(X) + &
!!$          13*fMG_RR_prime(X))*Y(3)*(Y(5) - Y(6))) - 27*X*derivative_conformal_Hubble_parameter(X)*(-1 + &
!!$          F_MG(X))*fMG_RR_prime(X)*(Y(5) - Y(6))*Y(6) + Y(1)*(9*derivative_conformal_Hubble_parameter(X)*fMG_RR_prime(X)*&
!!$          (2*X - 2*X*F_MG(X) + X**2*fMG_R_prime(X) + 4*wavenumber_k**2*fMG_RR_prime(X))*Y(5) - &
!!$          (-(X**2*derivative_conformal_Hubble_parameter(X)*(-1 + F_MG(X))*fMG_RR_double_prime(X)) + &
!!$          18*wavenumber_k**2*derivative_conformal_Hubble_parameter(X)*fMG_RR_prime(X)**2 + &
!!$          X*fMG_RR_prime(X)*(derivative_conformal_Hubble_parameter(X)*(7 - 7*F_MG(X) + 9*X*fMG_R_prime(X)) - &
!!$          2*X*(-1 + F_MG(X))*second_derivative_conformal_Hubble_parameter(X)))*Y(6))) + &
!!$          9*X**2*derivative_conformal_Hubble_parameter(X)*(-1 + F_MG(X))*fMG_RR_prime(X)*Y(1)*(Y(5) - &
!!$          Y(6))*derivative_Omega_Matter(X))) - 2*wavenumber_k**2*fMG_RR(X)**2*(2*X**2*Omega_Matter(X)*&
!!$          (20*derivative_conformal_Hubble_parameter(X) + 3*X*second_derivative_conformal_Hubble_parameter(X))*Y(1)*&
!!$          (2*Y(5) - Y(6))*derivative_Omega_DE(X) + Omega_DE(X)*(Omega_Matter(X)*(12*Y(3)*(2*Y(5) - Y(6)) - &
!!$          X*(-6*(20*derivative_conformal_Hubble_parameter(X) + 3*X*second_derivative_conformal_Hubble_parameter(X))*&
!!$          (2*Y(5) - Y(6))*Y(6) + Y(1)*(-16*(20*derivative_conformal_Hubble_parameter(X) + &
!!$          3*X*second_derivative_conformal_Hubble_parameter(X))*Y(5) + (45*derivative_conformal_Hubble_parameter(X) + &
!!$          13*X*second_derivative_conformal_Hubble_parameter(X) + X**2*third_derivative_conformal_Hubble_parameter(X))*Y(6)))) - & 
!!$          2*X**2*(20*derivative_conformal_Hubble_parameter(X) + 3*X*second_derivative_conformal_Hubble_parameter(X))*Y(1)*&
!!$          (2*Y(5) - Y(6))*derivative_Omega_Matter(X))) + X*fMG_RR(X)*(H0**2*X**2*Omega_DE(X)**2*&
!!$          (-4*(20*derivative_conformal_Hubble_parameter(X) + 3*X*second_derivative_conformal_Hubble_parameter(X))*Y(5) + &
!!$          (5*derivative_conformal_Hubble_parameter(X) + X*(7*second_derivative_conformal_Hubble_parameter(X) + &
!!$          X*third_derivative_conformal_Hubble_parameter(X)))*Y(6)) + 2*X**2*Omega_Matter(X)*Y(1)*&
!!$          ((-2*derivative_conformal_Hubble_parameter(X)*(-10*X + 10*X*F_MG(X) + 9*wavenumber_k**2*fMG_RR_prime(X)) - &
!!$          3*X**2*(-1 + F_MG(X))*second_derivative_conformal_Hubble_parameter(X))*Y(5) + &
!!$          (derivative_conformal_Hubble_parameter(X)*(-20*X + 20*X*F_MG(X) + 9*wavenumber_k**2*fMG_RR_prime(X)) + &
!!$          3*X**2*(-1 + F_MG(X))*second_derivative_conformal_Hubble_parameter(X))*Y(6))*derivative_Omega_DE(X) + &
!!$          Omega_DE(X)*(Omega_Matter(X)*(-2*Y(3)*(2*(-3*X + 3*X*F_MG(X) + wavenumber_k**2*X*fMG_RR_double_prime(X) + &
!!$          13*wavenumber_k**2*fMG_RR_prime(X))*Y(5) - (-6*X + 6*X*F_MG(X) + wavenumber_k**2*X*fMG_RR_double_prime(X) + &
!!$          13*wavenumber_k**2*fMG_RR_prime(X))*Y(6)) + X*(6*Y(6)*((-2*derivative_conformal_Hubble_parameter(X)*&
!!$          (-10*X + 10*X*F_MG(X) + 9*wavenumber_k**2*fMG_RR_prime(X)) - 3*X**2*(-1 + F_MG(X))*&
!!$          second_derivative_conformal_Hubble_parameter(X))*Y(5) + (derivative_conformal_Hubble_parameter(X)*&
!!$          (-20*X + 20*X*F_MG(X) + 9*wavenumber_k**2*fMG_RR_prime(X)) + 3*X**2*(-1 + F_MG(X))*&
!!$          second_derivative_conformal_Hubble_parameter(X))*Y(6)) + Y(1)*(2*(4*derivative_conformal_Hubble_parameter(X)*&
!!$          (10*X - 10*X*F_MG(X) + 5*X**2*fMG_R_prime(X) + 2*wavenumber_k**2*fMG_RR_prime(X)) + &
!!$          3*X*(2*X - 2*X*F_MG(X) + X**2*fMG_R_prime(X) + 4*wavenumber_k**2*fMG_RR_prime(X))*&
!!$          second_derivative_conformal_Hubble_parameter(X))*Y(5) + (derivative_conformal_Hubble_parameter(X)*&
!!$          (-5*X + 5*X*F_MG(X) - 40*X**2*fMG_R_prime(X) + 2*wavenumber_k**2*X*fMG_RR_double_prime(X) - &
!!$          48*wavenumber_k**2*fMG_RR_prime(X)) - X*((7*X - 7*X*F_MG(X) + 6*X**2*fMG_R_prime(X) + &
!!$          8*wavenumber_k**2*fMG_RR_prime(X))*second_derivative_conformal_Hubble_parameter(X) - &
!!$          X**2*(-1 + F_MG(X))*third_derivative_conformal_Hubble_parameter(X)))*Y(6)))) + &
!!$          2*X**2*Y(1)*((2*derivative_conformal_Hubble_parameter(X)*(-10*X + 10*X*F_MG(X) + &
!!$          9*wavenumber_k**2*fMG_RR_prime(X)) + 3*X**2*(-1 + F_MG(X))*&
!!$          second_derivative_conformal_Hubble_parameter(X))*Y(5) + (derivative_conformal_Hubble_parameter(X)*&
!!$          (20*X - 20*X*F_MG(X) - 9*wavenumber_k**2*fMG_RR_prime(X)) - 3*X**2*(-1 + F_MG(X))*&
!!$          second_derivative_conformal_Hubble_parameter(X))*Y(6))*derivative_Omega_Matter(X)))))/&
!!$          (18.*H0**2*X**3*conformal_Hubble_parameter(X)**2*Omega_DE(X)**3*(H0**2*X**2*Omega_DE(X) + &
!!$           (-X**2 + X**2*F_MG(X) + 2*wavenumber_k**2*fMG_RR(X))*Omega_Matter(X)*Y(1)))
!!$
!!$     F(5) = -(3*H0**2*X**2*Omega_Matter(X)*Y(1) + 3*H0**2*X**2*Omega_DE(X)*Y(2) +  2*wavenumber_k**2*Y(5) &
!!$          + 6*conformal_Hubble_parameter(X)**2*Y(6))/(6.*X*conformal_Hubble_parameter(X)**2)
!!$
!!$     F(6) = -(3*H0**4*X**4*Omega_DE(X)**3*Y(2) +  6*X*conformal_Hubble_parameter(X)**2*&
!!$          derivative_Omega_DE(X)*Omega_Matter(X)*Y(1)*((-X**2 + X**2*F_MG(X) + 4*wavenumber_k**2*fMG_RR(X))*Y(5) + & 
!!$          (X**2 - X**2*F_MG(X) - 2*wavenumber_k**2*fMG_RR(X))*Y(6)) + H0**2*X**2*Omega_DE(X)**2*&
!!$          (2*((wavenumber_k**2 + 6*conformal_Hubble_parameter(X)**2)*Y(5) - 3*conformal_Hubble_parameter(X)**2*Y(6)) + &
!!$          3*Omega_Matter(X)*(Y(1)*(H0**2*X**2 + (-X**2 + X**2*F_MG(X) + 4*wavenumber_k**2*fMG_RR(X))*Y(2)) + &
!!$          3*Y(2)*((-X**2 + X**2*F_MG(X) + 4*wavenumber_k**2*fMG_RR(X))*Y(5) + (X**2 - X**2*F_MG(X) - &
!!$          2*wavenumber_k**2*fMG_RR(X))*Y(6)))) + Omega_DE(X)*(3*H0**2*X**2*Omega_Matter(X)**2*Y(1)*&
!!$          ((-X**2 + X**2*F_MG(X) + 4*wavenumber_k**2*fMG_RR(X))*Y(1) + 3*((-X**2 + X**2*F_MG(X) + &
!!$          4*wavenumber_k**2*fMG_RR(X))*Y(5) + (X**2 - X**2*F_MG(X) - 2*wavenumber_k**2*fMG_RR(X))*Y(6))) - &
!!$          2*Omega_Matter(X)*(-3*conformal_Hubble_parameter(X)*Y(3)*((-X**2 + X**2*F_MG(X) + &
!!$          4*wavenumber_k**2*fMG_RR(X))*Y(5) + (X**2 - X**2*F_MG(X) - 2*wavenumber_k**2*fMG_RR(X))*Y(6)) + &
!!$          wavenumber_k**2*Y(5)*((X**2 - X**2*F_MG(X) - 4*wavenumber_k**2*fMG_RR(X))*Y(1) + &
!!$          3*((X**2 - X**2*F_MG(X) - 4*wavenumber_k**2*fMG_RR(X))*Y(5) + (-X**2 + X**2*F_MG(X) + &
!!$          2*wavenumber_k**2*fMG_RR(X))*Y(6))) + 3*conformal_Hubble_parameter(X)**2*(3*Y(6)*((X**2 - &
!!$          X**2*F_MG(X) - 4*wavenumber_k**2*fMG_RR(X))*Y(5) + (-X**2 + X**2*F_MG(X) + 2*wavenumber_k**2*fMG_RR(X))*Y(6)) + & 
!!$          Y(1)*((2*X**2 - 2*X**2*F_MG(X) + X**3*fMG_R_prime(X) - 16*wavenumber_k**2*fMG_RR(X) + &
!!$          4*wavenumber_k**2*X*fMG_RR_prime(X))*Y(5) - (X**2 - X**2*F_MG(X) + X**3*fMG_R_prime(X) - &
!!$          4*wavenumber_k**2*fMG_RR(X) + 2*wavenumber_k**2*X*fMG_RR_prime(X))*Y(6)))) - &
!!$          6*X*conformal_Hubble_parameter(X)**2*Y(1)*((-X**2 + X**2*F_MG(X) + 4*wavenumber_k**2*fMG_RR(X))*Y(5) + &
!!$          (X**2 - X**2*F_MG(X) - 2*wavenumber_k**2*fMG_RR(X))*Y(6))*derivative_Omega_Matter(X)))/&
!!$          (6.*X*conformal_Hubble_parameter(X)**2*Omega_DE(X)*(H0**2*X**2*Omega_DE(X) + &
!!$          (-X**2 + X**2*F_MG(X) + 2*wavenumber_k**2*fMG_RR(X))*Omega_Matter(X)*Y(1)))


  Else if ( ( MG_parametrisation .eq. 'Savvas' ) .and. (approach .eq. 'GI') ) then

     ! \Delta_m_prime STARTS
     F(1) = (-g_K(X)*E_H(X)*Y(3) + 3.d0*( Omega_Matter(X)*Y(3) + Omega_DE(X)*Y(4) )/2.d0)/X
     ! \Delta_m prime ENDS

     ! \Delta_de_prime STARTS
     F(2) = (-3.d0*Y(2) - g_K(X)*E_H(X)*Y(4) - 2.d0*wavenumber_k**2/3.d0/equation_of_state(X)/&
          g_K(X)/E_H(X)/conformal_Hubble_parameter(X)**2*( Y(2) - Y(4)*X*fMG_R_prime(X)/2.d0/(1.d0 + fMG_R(X)) &
          + Omega_Matter(X)*fMG_R(X)*Y(1)/Omega_DE(X)/(1.d0 + fMG_R(X))  &
          - Omega_Matter(X)*X*fMG_R_prime(X)/Omega_DE(X)/2.d0/(1.d0 + fMG_R(X))*Y(3)   ))/X
     ! \delta_prime ENDS

     ! Theta_m_prime STARTS
     F(3) = (-E_H(X)*Y(3) + 3.d0*Y(5) - 3.d0/(1.d0 + fMG_R(X))*( -Omega_DE(X)/g_K(X)/E_H(X)*Y(2) &
          -fMG_R(X)/g_K(X)/E_H(X)*( Omega_DE(X)*(Y(2) - X*fMG_R_prime(X)*Y(4)/2.d0/fMG_R(X))&
          + Omega_Matter(X)*(Y(1) - X*fMG_R_prime(X)*Y(3)/2.d0/fMG_R(X) )   )  ))/X
     ! Theta_m_prime ENDS

     ! Theta_de_prime STARTS
     F(4) = (-3.d0*Y(2) - E_H(X)*Y(4) - 2.d0*wavenumber_k**2/3.d0/equation_of_state(X)/&
          g_K(X)/E_H(X)/conformal_Hubble_parameter(X)**2*( Y(2) - Y(4)*X*fMG_R_prime(X)/2.d0/(1.d0 + fMG_R(X))&
          + Omega_Matter(X)*fMG_R(X)*Y(1)/Omega_DE(X)/(1.d0 + fMG_R(X))  &
          - Omega_Matter(X)*X*fMG_R_prime(X)/Omega_DE(X)/2.d0/(1.d0 + fMG_R(X))*Y(3)) &
          - 3.d0/equation_of_state(X)*( (Zeta_de(X) - (2.d0*(1.d0 + fMG_R(X)) - &
          X*fMG_R_prime(X) )/3.d0/g_K(X)/X/fMG_R_prime(X) )*Y(2) &
          - Zeta_de(X)*Y(4) + Omega_Matter(X)/Omega_DE(X)*( Zeta_m(X) - &
          (2.d0*fMG_R(X) - X*fMG_R_prime(X) )/3.d0/g_K(X)/X/fMG_R_prime(X) )*Y(1) &
          - Omega_Matter(X)/Omega_DE(X)*Zeta_m(X)*Y(3)  ))/X
     ! Theta_de_prime ENDS

     ! Z_prime STARTS
     F(5) = (( Omega_Matter(X)*Y(3) + Omega_DE(X)*Y(4) )/2.d0 &
          - Y(5) + 1.d0/(1.d0 + fMG_R(X))*( -Omega_DE(X)/g_K(X)/E_H(X)*Y(2) &
          -fMG_R(X)/g_K(X)/E_H(X)*( Omega_DE(X)*(Y(2) - X*fMG_R_prime(X)*Y(4)/2.d0/fMG_R(X))&
          + Omega_Matter(X)*(Y(1) - X*fMG_R_prime(X)*Y(3)/2.d0/fMG_R(X) )   )  ))/X
     ! Z_prime ENDS

  Else if ( ( (MG_parametrisation .eq. 'HS_Basilakos') .or. (MG_parametrisation .eq. 'Starobinsky_Basilakos') )&
       .and. (approach .eq. 'EF')  ) then

     If ( abs(fMG_R(X)) .ge. switch_GR_equations) then

        ! \delta_m_prime STARTS
        F(1) = -(3*H0**2*X**2*(Omega_Matter(X)*Y(1) + Omega_DE(X)*Y(2)) + &
             2*wavenumber_k*(conformal_Hubble_parameter(X)*Y(3) + wavenumber_k*Y(5)) + &
             6*conformal_Hubble_parameter(X)**2*Y(6))/(2.*X*conformal_Hubble_parameter(X)**2)
        ! \delta_m prime ENDS

        ! \delta_prime STARTS
        F(2) = (-2*wavenumber_k*conformal_Hubble_parameter(X)*Omega_DE(X)*Y(4) + &
             (1 + equation_of_state(X))*Omega_DE(X)*(-3*H0**2*X**2*(Omega_Matter(X)*Y(1) + &
             Omega_DE(X)*Y(2)) - 2*wavenumber_k**2*Y(5)) + 6*conformal_Hubble_parameter(X)**2*&
             (-(Omega_Matter(X)*pressure_perturbation_over_density(X)*Y(1)) + &
             Omega_DE(X)*(equation_of_state(X)*Y(2) - (1 + equation_of_state(X))*Y(6))))/&
             (2.*X*conformal_Hubble_parameter(X)**2*Omega_DE(X))
        ! \delta_prime ENDS

        ! V_m_prime STARTS
        F(3) = (-(conformal_Hubble_parameter(X)*Y(3)) + wavenumber_k*Y(6))/&
             (X*conformal_Hubble_parameter(X))
        ! V_m_prime ENDS

        ! V_prime STARTS
        F(4) = (wavenumber_k*Omega_Matter(X)*(-2*anisotropic_stress(X) + 3*pressure_perturbation_over_density(X))*Y(1) + &
             3*Omega_DE(X)*(conformal_Hubble_parameter(X)*(-1 + 3*equation_of_state(X))*Y(4) + &
             wavenumber_k*(1 + equation_of_state(X))*Y(6)))/(3.*X*conformal_Hubble_parameter(X)*Omega_DE(X))
        ! V_prime ENDS

        ! \phi_prime STARTS
        F(5) = -(3*H0**2*X**2*(Omega_Matter(X)*Y(1) + Omega_DE(X)*Y(2)) + 2*wavenumber_k**2*Y(5) + &
             6*conformal_Hubble_parameter(X)**2*Y(6))/(6.*X*conformal_Hubble_parameter(X)**2)
        ! \phi_prime ENDS

        ! \psi_prime STARTS
        F(6) = (27*H0**4*X**4*anisotropic_stress(X)*Omega_Matter(X)**2*Y(1) + &
             3*H0**2*X**2*Omega_DE(X)*(-wavenumber_k**2 + 9*H0**2*X**2*anisotropic_stress(X)*&
             Omega_Matter(X))*Y(2) - 2*((wavenumber_k**4 + 6*wavenumber_k**2*&
             conformal_Hubble_parameter(X)**2)*Y(5) - 3*conformal_Hubble_parameter(X)**2*&
             (9*H0**2*X**2*anisotropic_stress(X)*Omega_Matter(X)*Y(1) + wavenumber_k**2*Y(6))) +& 
             3*H0**2*X**2*Omega_Matter(X)*(-((wavenumber_k**2 + 6*X*conformal_Hubble_parameter(X)**2*&
             derivative_anisotropic_stress(X))*Y(1)) + 6*anisotropic_stress(X)*&
             (wavenumber_k*conformal_Hubble_parameter(X)*Y(3) + wavenumber_k**2*Y(5) +& 
             3*conformal_Hubble_parameter(X)**2*Y(6))))/(6.*wavenumber_k**2*X*conformal_Hubble_parameter(X)**2) 
        ! \psi_prime ENDS

     Else

        ! \delta_m_prime STARTS
        F(1) = -(3*H0**2*X**2*Omega_Matter(X)*Y(1) + 2*wavenumber_k*(conformal_Hubble_parameter(X)*Y(3) + &
             wavenumber_k*Y(5)) + 6*conformal_Hubble_parameter(X)**2*Y(6))/(2.*X*conformal_Hubble_parameter(X)**2)
        ! \delta_m prime ENDS

        ! \delta_prime STARTS
        F(2) = 0.d0
        ! \delta_prime ENDS

        ! V_m_prime STARTS
        F(3) = (-(conformal_Hubble_parameter(X)*Y(3)) + wavenumber_k*Y(6))/(X*conformal_Hubble_parameter(X))
        ! V_m_prime ENDS

        ! V_prime STARTS
        F(4) = 0.d0
        ! V_prime ENDS

        ! \phi_prime STARTS
        F(5) = -(3*H0**2*X**2*Omega_Matter(X)*Y(1) + 2*wavenumber_k**2*Y(5) +& 
             6*conformal_Hubble_parameter(X)**2*Y(6))/(6.*X*conformal_Hubble_parameter(X)**2)
        ! \phi_prime ENDS

        ! \psi_prime STARTS
        F(6) = (-3*H0**2*X**2*Omega_Matter(X)*Y(1) - 2*(wavenumber_k**2 + 6*conformal_Hubble_parameter(X)**2)*Y(5) + &
             6*conformal_Hubble_parameter(X)**2*Y(6))/(6.*X*conformal_Hubble_parameter(X)**2)
        ! \psi_prime ENDS

     End if

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
