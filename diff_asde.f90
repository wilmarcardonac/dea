Program time_evolution

  !################
  ! REQUIRED MODULES
  !################

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
  ! MATTER AND DARK ENERGY PERTURBATIONS: \delta_cdm, \delta_de, \v_cdm, \v_de
  DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK),RPAR(number_of_parameters)      ! Y IS AN ARRAY HOLDING THE PERTURBATION VARIABLES, RPAR HOLDS PARAMETERS 
  ! IN THE SYSTEM OF DIFFERENTIAL EQUATIONS 
  EXTERNAL RHSPER,JRHSPER,SOLOUT                         ! SUBROUTINES FROM THE SOLVER
       
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

  ! --- PARAMETERS IN THE DIFFERENTIAL EQUATION                            

  Do m=1,number_of_parameters

     If (m .eq. 1) then

        RPAR(m) = w0_fld ! EQUATION OF STATE

     else if (m .eq. 2) then

        RPAR(m) = (omega_b + omega_cdm)/(H0*1.d-2*speed_of_light)**2 ! Omega_m

     else if (m .eq. 3) then

        RPAR(m) = sqrt(cs2_fld) ! c_s IS A FRACTION OF THE SPEED OF LIGHT 
     
     else if (m .eq. 4) then

        RPAR(m) = wavenumber_k ! k
        
     else if ( m .eq. 5 ) then

        RPAR(m) = e_pi !e

     else if ( m .eq. 6 ) then

        RPAR(m) = f_pi !f

     else if ( m .eq. 7 ) then
        
        RPAR(m) = g_pi !g
     
     End if

  End Do

  write(UNIT_EXE_FILE,*) 'STARTING ANALYSIS. PARAMETERS FOR CURRENT RUN ARE AS FOLLOWS: '

  write(UNIT_EXE_FILE,*) 'EQUATION OF STATE: ', RPAR(1)

  write(UNIT_EXE_FILE,*) 'MATTER DENSITY PARAMETER: ', RPAR(2)

  write(UNIT_EXE_FILE,*) 'SOUND SPEED: ', RPAR(3)

  write(UNIT_EXE_FILE,*) 'WAVENUMBER: ', RPAR(4)

  write(UNIT_EXE_FILE,*) 'DEA PARAMETER e_\pi : ', RPAR(5)

  write(UNIT_EXE_FILE,*) 'DEA PARAMETER f_\pi : ', RPAR(6)

  write(UNIT_EXE_FILE,*) 'DEA PARAMETER g_\pi : ', RPAR(7)

  If ( ceff(RPAR(3),RPAR(6))**2 .ge. 0 ) Then

     write(UNIT_EXE_FILE,*) 'EFFECTIVE SOUND SPEED SQUARED FOR CURRENT MODEL IS : ', ceff(RPAR(3),RPAR(6))**2

  Else

     write(UNIT_EXE_FILE,*) 'EFFECTIVE SOUND SPEED SQUARED FOR CURRENT MODEL IS : ', ceff(RPAR(3),RPAR(6))**2

     write(UNIT_EXE_FILE,*) 'AND WILL SURELY LEAD TO INSTABILITIES IN THE PERTURBATIONS'

  End If

  N = dimension_system_ode ! DIMENSION OF THE SYSTEM 
                                 
  IJAC = 0 ! COMPUTE THE JACOBIAN ANALYTICALLY. '0' MEANS WHAT ??                                 
        
  MLJAC = N ! JACOBIAN IS A FULL MATRIX                                         

  IMAS = 0 ! DIFFERENTIAL EQUATION IS IN EXPLICIT FORM. '0' MEANS WHAT ??                         
        
  IOUT = 1 ! OUTPUT ROUTINE IS USED DURING INTEGRATION                         
         
  X = initial_scale_factor  ! INITIAL VALUE OF THE SCALE FACTOR
             
  write(UNIT_EXE_FILE,*) 'INITIAL CONDITIONS ARE SET AT SCALE FACTOR : ', X

  write(UNIT_EXE_FILE,*) 'SYSTEM OF DIFFERENTIAL EQUATIONS IS WRITTEN IN THE CODE AS FOLLOWS: '

  write(UNIT_EXE_FILE,*) 'Y(1) IS \delta_m, MATTER DENSITY PERTURBATION '

  write(UNIT_EXE_FILE,*) 'Y(2) IS \delta_de, DARK ENERGY DENSITY PERTURBATION '

  write(UNIT_EXE_FILE,*) 'Y(3) IS \theta_m, MATTER VELOCITY PERTURBATION '

  write(UNIT_EXE_FILE,*) 'Y(4) IS \theta_de, DARK ENERGY VELOCITY PERTURBATION '


  write(UNIT_EXE_FILE,*) 'THE CONFORMAL HUBBLE PARAMETER H(a) AT ', X, ' IS : ',cHubble(X,H0,RPAR(2),RPAR(1)) 

  write(UNIT_EXE_FILE,*) 'THE WAVENUMBER K CORRESPONDING TO THE HORIZON AT INITIAL SCALE FACTOR ', X, ' IS : ',&
       cHubble(X,H0,RPAR(2),RPAR(1))/speedL

  write(UNIT_EXE_FILE,*) 'THE CONFORMAL HUBBLE PARAMETER AT THE PRESENT TIME IS : ', cHubble(1.d0,H0,RPAR(2),RPAR(1))

  write(UNIT_EXE_FILE,*) 'THE WAVENUMBER K CORRESPONDING TO THE HORIZON AT THE PRESENT TIME IS : ', &
       cHubble(1.d0,H0,RPAR(2),RPAR(1))/speedL
  
  write(UNIT_EXE_FILE,*) 'THE WAVENUMBER K FOR THE CURRENT MODE IS : ', RPAR(4)

  write(UNIT_EXE_FILE,*) 'THE MODE CROSSES THE HORIZON (IN MATTER DOMINANCE) AT SCALE FACTOR : ', H1(H0,RPAR(4),speedL,RPAR(2)) 

  write(UNIT_EXE_FILE,*) 'THE MODE CROSSES THE EFFECTIVE SOUND HORIZON (IN MATTER DOMINANCE) AT SCALE FACTOR : ', &
       H2(H0,RPAR(4),RPAR(3),RPAR(6),RPAR(2))

  write(UNIT_OUTPUT_FILE,*) '# scale_factor        \delta_m             \delta_de             \theta_m        '//trim(' ')//&
'\theta_de            \phi                 \psi'

  !######################################################################################
  ! SETTING INITIAL CONDITIONS FOR THE PERTURBATIONS USING DOMENICO AND MARTIN'S SOLUTION
  !######################################################################################

  Y(1) = dm(X,H0,RPAR(4),Bf2(RPAR(4),kp,ks,As,ns,Pi),RPAR(2)) ! MATTER DENSITY PERTURBATIONS
  
  Y(3) = thetammd(X,H0,RPAR(4),Bf2(RPAR(4),kp,ks,As,ns,Pi),RPAR(2)) ! MATTER VELOCITY PERTURBATIONS

  If (RPAR(4)/cHubble(X,H0,RPAR(2),RPAR(1))>(1/ceff(RPAR(3),RPAR(6)))) Then

     Y(2) = dd2(RPAR(1),RPAR(3),Bf2(RPAR(4),kp,ks,As,ns,Pi)) ! DARK ENERGY DENSITY PERTURBATIONS

     Y(4) = thetademdsstsh(X,H0,RPAR(4),RPAR(1),RPAR(3),Bf2(RPAR(4),kp,ks,As,ns,Pi)) ! DARK ENERGY VELOCITY PERTURBATIONS

     write(UNIT_EXE_FILE,*) 'CURRENT MODE IS SUB-SOUND-HORIZON AT THE INITIAL SCALE FACTOR'

  Else

     Y(2) = dd1(X,H0,RPAR(4),RPAR(1),Bf2(RPAR(4),kp,ks,As,ns,Pi),RPAR(2)) ! DARK ENERGY DENSITY PERTURBATIONS

     Y(4) = thetademdsltsh(X,H0,RPAR(4),RPAR(1),Bf2(RPAR(4),kp,ks,As,ns,Pi),RPAR(2)) ! DARK ENERGY VELOCITY PERTURBATIONS

     write(UNIT_EXE_FILE,*) 'CURRENT MODEL IS SUPER-SOUND-HORIZON AT THE INITIAL SCALE FACTOR'

  End If

  !################################################################################################################
  ! SOLUTIONS ON BOTH SUPER-HORIZON AND SUB-HORIZON SCALES FOR DARK ENERGY PERTURBATIONS IN MATTER DOMINATED REGIME
  !################################################################################################################

  If ( ceff(RPAR(3),RPAR(6))**2 .ge. 0) then

     write(UNIT_EXE_FILE,*) 'WRITING ANALYTICAL SOLUTIONS FOR THE CURRENT MODE'

     write(UNIT_OUTPUT_FILE2,*) '# scale_factor    \delta_{de}^{sup-hor}    v_{de}^{sup-hor}   '//trim(' ')//&
          '\delta_{de}^{sub-sound}    v_{de}^{sub-sound}'

     Do m=1,101

        Z=10.d0**(-5.d0+Real(m-1)/100.d0*5.d0)

        write(UNIT_OUTPUT_FILE2,89) Z,dde1(Z,H0,RPAR(4),Bf2(RPAR(4),kp,ks,As,ns,Pi),RPAR(2),&
             RPAR(1),RPAR(5),RPAR(6),RPAR(3),RPAR(7)),&
             -Vd1(Z,H0,RPAR(4),Bf2(RPAR(4),kp,ks,As,ns,Pi),RPAR(2),RPAR(1),RPAR(5),RPAR(6),RPAR(3),&
             RPAR(7))*H0*sqrt(RPAR(2))/((1.d0+RPAR(1))*sqrt(Z)*RPAR(4)),&
             dd4(Z,H0,RPAR(4),Bf2(RPAR(4),kp,ks,As,ns,Pi),RPAR(2),RPAR(1),RPAR(5),RPAR(6),RPAR(7),RPAR(3)),&
             -Vd4(Z,H0,RPAR(4),Bf2(RPAR(4),kp,ks,As,ns,Pi),RPAR(2),RPAR(1),RPAR(5),RPAR(6),RPAR(7),&
             RPAR(3))*H0*sqrt(RPAR(2))/(RPAR(4)*(1.d0+RPAR(1))*sqrt(Z))

89      Format(E20.10,E20.10,E20.10,E20.10,E20.10)

     End do

  Else

     write(UNIT_EXE_FILE,*) 'EFFECTIVE SQUARE SOUND SPEED IS NEGATIVE. NOT ANALYTICAL SOLUTIONS WRITTEN'

  End If

  !################################################################
  ! ENDPOINT INTEGRATION, TOLERANCE FOR SOLUTIONS, INITIAL STEPSIZE
  !################################################################

  XEND = final_scale_factor ! ENDPOINT OF INTEGRATION                                           

  RTOL=1.0D-14              ! REQUIRED TOLERANCE
  ATOL=1.0D0*RTOL           ! REQUIRED TOLERANCE 
  ITOL=0                    ! REQUIRED TOLERANCE                                     
                                               
  H=1.0D-19                 ! INITIAL STEP SIZE

  DO I=1,20                 ! SET DEFAULT VALUES                                                

     IWORK(I)=0 

     WORK(I)=0.D0 

  END DO

  !##############################################################
  ! CALL SUBROUTINE TO SOLVE THE SYSTEM OF DIFFERENTIAL EQUATIONS
  !##############################################################

  CALL RADAU5(N,RHSPER,X,Y,XEND,H,RTOL,ATOL,ITOL,JRHSPER,IJAC,MLJAC,MUJAC,RHSPER,IMAS,MLMAS,MUMAS,&
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
  DIMENSION Y(N),CONT(LRC),RPAR(7) 
  COMMON /INTERN/XOUT 

  IF (NR.EQ.1) THEN 

     WRITE (UNIT_OUTPUT_FILE,99) X,Y(1),Y(2),Y(3)/RPAR(4),Y(4)/RPAR(4),&                        
          phi(X,RPAR(4),H0,RPAR(2),RPAR(1),Y(1),Y(2),Y(3),Y(4)),&
          psi(X,RPAR(4),H0,RPAR(2),RPAR(1),Y(1),Y(2),Y(3),Y(4),RPAR(5),RPAR(6),RPAR(7))

     XOUT = 1.0D-5 !1.0D0/1001.0D0 

  ELSE 

10   CONTINUE 

     IF (X.GE.XOUT) THEN 

        ! --- CONTINUOUS OUTPUT FOR RADAU5                                      
        WRITE (UNIT_OUTPUT_FILE,99) XOUT,CONTR5(1,XOUT,CONT,LRC),                &
             CONTR5(2,XOUT,CONT,LRC),                     &    
             CONTR5(3,XOUT,CONT,LRC)/RPAR(4),                     &
             CONTR5(4,XOUT,CONT,LRC)/RPAR(4), &                     
             phi(XOUT,RPAR(4),H0,RPAR(2),RPAR(1),CONTR5(1,XOUT,CONT,LRC),&
             CONTR5(2,XOUT,CONT,LRC),CONTR5(3,XOUT,CONT,LRC),CONTR5(4,XOUT,CONT,LRC)), &
             psi(XOUT,RPAR(4),H0,RPAR(2),RPAR(1),CONTR5(1,XOUT,CONT,LRC),CONTR5(2,XOUT,CONT,LRC),&
             CONTR5(3,XOUT,CONT,LRC),CONTR5(4,XOUT,CONT,LRC),&
             RPAR(5),RPAR(6),RPAR(7))  

        XOUT = XOUT+1.0D-6              

        GOTO 10 

     END IF

  END IF

99 FORMAT(E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,E20.10)

  RETURN 

END SUBROUTINE SOLOUT                                   

subroutine RHSPER(N,X,Y,F,RPAR,IPAR)
           
  use fiducial
  use functions

  IMPLICIT REAL*8 (A-H,O-Z)
  Dimension Y(N),F(N),RPAR(7)
 
           
  F(1) = (X**(-3.d0 - 6.d0*RPAR(1))*(9.d0*H0**3*RPAR(4)**4*X**(1.d0 + 3.d0*RPAR(1))*(-2.d0*(-1.d0 + &
       RPAR(2))*RPAR(5) + RPAR(2)*X**(3.d0*RPAR(1)))*(1.d0 + RPAR(2)*(-1.d0&
       + X**(3.d0*RPAR(1))))*Y(1) - 9.d0*H0**3*(-1.d0 + RPAR(2))*RPAR(4)**2*(1.d0 + RPAR(2)*(-1.d0 + &
       X**(3.d0*RPAR(1))))*(RPAR(4)**2*(1.d0 + 2.d0*RPAR(6))*X**(1.d0 + 3.d0*RPAR(1))&
       + 2.d0*H0**2*RPAR(7)*(1.d0 + RPAR(2)*(-1.d0 + X**(3.d0*RPAR(1)))))*Y(2) + Sqrt(X**(-1.d0 - &
       3.d0*RPAR(1))*(1.d0 + RPAR(2)*(-1.d0 + X**(3.d0*RPAR(1)))))*(RPAR(4)**2*X**(1.d0&
       + 3.d0*RPAR(1))*(9.d0*H0**2*RPAR(2)*RPAR(4)**2*X**(1.d0 + 6.d0*RPAR(1)) - 2.d0*RPAR(4)**4*X**(2.d0 + &
       6.d0*RPAR(1)) + 27.d0*H0**4*(2.d0*RPAR(5) - 2.d0*RPAR(2)*RPAR(5) +&
       RPAR(2)*X**(3.d0*RPAR(1)))*(1.d0 + RPAR(2)*(-1.d0 + X**(3.d0*RPAR(1)))))*Y(3) - 9.d0*H0**2*(1.d0 + &
       RPAR(1))*(-1.d0 + RPAR(2))*(RPAR(4)**4*X**(2.d0 + 6.d0*RPAR(1)) +&
       3.d0*H0**2*RPAR(4)**2*(1.d0 + 2.d0*RPAR(6))*X**(1.d0 + 3.d0*RPAR(1))*(1.d0 + RPAR(2)*(-1.d0 + X**(3.d0*RPAR(1)))) + &
       6.d0*H0**4*RPAR(7)*(1.d0 + RPAR(2)*(-1.d0 +&
       X**(3.d0*RPAR(1))))**2)*Y(4))))/(2.d0*H0*RPAR(4)**6*(1.d0 + RPAR(2)*(-1.d0 + X**(3.d0*RPAR(1)))))

  F(2) = (X**(-4.d0 - 9.d0*RPAR(1))*(3.d0*H0*RPAR(4)**2*X**(1.d0 + 3.d0*RPAR(1))*(6.d0*H0**4*(1.d0 + RPAR(1))*(-1.d0 + &
       RPAR(2))**2*RPAR(7) + 2.d0*(RPAR(1) - &
       RPAR(3)**2)*RPAR(4)**4*X**(2.d0 + 6.d0*RPAR(1)) - 3.d0*H0**2*(1.d0 + RPAR(1))*(-1.d0 + &
       RPAR(2))*X**(3.d0*RPAR(1))*(2.d0*H0**2*RPAR(2)*RPAR(7) + RPAR(4)**2*(1.d0 &
       + 2.d0*RPAR(6))*X))*Sqrt(X**(-1.d0 - 3.d0*RPAR(1))*(1.d0 + RPAR(2)*(-1.d0 + X**(3.d0*RPAR(1)))))*Y(2) + &
       (1.d0 + RPAR(1))*(9.d0*H0**3*RPAR(4)**4*X**(2.d0 + & 
       6.d0*RPAR(1))*(-2.d0*(-1.d0 + RPAR(2))*RPAR(5) + RPAR(2)*X**(3.d0*RPAR(1)))*Sqrt(X**(-1.d0 - &
       3.d0*RPAR(1))*(1.d0 + RPAR(2)*(-1.d0 + X**(3.d0*RPAR(1)))))*Y(1)&
       + 9.d0*H0**2*RPAR(4)**2*X**(1.d0 + 3.d0*RPAR(1))*(RPAR(2)*RPAR(4)**2*X**(1.d0 + 6.d0*RPAR(1)) + &
       3.d0*H0**2*(-2.d0*(-1.d0 + RPAR(2))*RPAR(5) + RPAR(2)*X**(3.d0*RPAR(1)))*(1.d0&
       + RPAR(2)*(-1.d0 + X**(3.d0*RPAR(1)))))*Y(3) + (-2.d0*RPAR(4)**6*X**(3.d0 + 9.d0*RPAR(1)) + &
       9.d0*H0**2*RPAR(4)**4*X**(2.d0 + 6.d0*RPAR(1))*((1.d0 - RPAR(2))*(1.d0 +&
       3.d0*RPAR(1) - 2.d0*RPAR(3)**2) + 2.d0*RPAR(2)*(RPAR(1) - RPAR(3)**2)*X**(3.d0*RPAR(1))) - &
       27.d0*H0**4*(1.d0 + RPAR(1))*(-1.d0 + RPAR(2))*RPAR(4)**2*(1.d0 + & 
       2.d0*RPAR(6))*X**(1.d0 + 3.d0*RPAR(1))*(1.d0 + RPAR(2)*(-1.d0 + X**(3.d0*RPAR(1)))) - 54.d0*H0**6*(1.d0 + &
       RPAR(1))*(-1.d0 + RPAR(2))*RPAR(7)*(1.d0 + &
       RPAR(2)*(-1.d0 + X**(3.d0*RPAR(1))))**2)*Y(4))))/(2.d0*H0*RPAR(4)**6*Sqrt(X**(-1.d0 - &
       3.d0*RPAR(1))*(1.d0 + RPAR(2)*(-1.d0 + X**(3.d0*RPAR(1))))))


  F(3) = (X**(-3.d0 - 6.d0*RPAR(1))*(-3.d0*H0*RPAR(4)**4*X**(1.d0 + 3.d0*RPAR(1))*(2.d0*RPAR(5) - &
       2.d0*RPAR(2)*RPAR(5) + RPAR(2)*X**(3.d0*RPAR(1)))*Y(1) + &
       3.d0*H0*(-1.d0 + RPAR(2))*RPAR(4)**2*(RPAR(4)**2*(1.d0 + 2.d0*RPAR(6))*X**(1.d0 + 3.d0*RPAR(1)) + &
       2.d0*H0**2*RPAR(7)*(1.d0 + RPAR(2)*(-1.d0 + X**(3.d0*RPAR(1)))))*Y(2) - &
       Sqrt(X**(-1.d0 - 3.d0*RPAR(1))*(1.d0 + RPAR(2)*(-1.d0 + X**(3.d0*RPAR(1)))))*(RPAR(4)**2*X**(1.d0 + &
       6.d0*RPAR(1))*(9.d0*H0**2*RPAR(2) + 2.d0*RPAR(4)**2*X)*Y(3) + &
       18.d0*H0**4*(1.d0 + RPAR(1))*(-1.d0 + RPAR(2))**2*RPAR(7)*Y(4) - 9.d0*H0**2*(-1.d0 + &
       RPAR(2))*X**(3.d0*RPAR(1))*(2.d0*RPAR(4)**2*RPAR(5)*X*Y(3) + &
       (1.d0 + RPAR(1))*(2.d0*H0**2*RPAR(2)*RPAR(7) + RPAR(4)**2*(1.d0 + &
       2.d0*RPAR(6))*X)*Y(4)))))/(2.d0*RPAR(4)**4*Sqrt((RPAR(2) - (-1.d0 + RPAR(2))/X**(3.d0*RPAR(1)))/X))

  F(4) = -(X**(-3.d0 - 6.d0*RPAR(1))*(RPAR(4)**4*X**(1.d0 + 3.d0*RPAR(1))*(4.d0*RPAR(4)**2*RPAR(5)*X**(1.d0 + &
       3.d0*RPAR(1)) + 9.d0*H0**2*(1.d0 + RPAR(1))*(2.d0*RPAR(5)&
       - 2.d0*RPAR(2)*RPAR(5) + RPAR(2)*X**(3.d0*RPAR(1))))*Y(1) + RPAR(4)**2*(-2.d0*RPAR(4)**4*(3.d0*RPAR(3)**2 - &
       2.d0*RPAR(6))*X**(2.d0 + 6.d0*RPAR(1)) + H0**2*RPAR(4)**2*X**(1.d0&
       + 3.d0*RPAR(1))*((1.d0 - RPAR(2))*(9.d0*(1.d0 + RPAR(1))*(1.d0 + 2.d0*RPAR(6)) + 4.d0*RPAR(7)) + &
       4.d0*RPAR(2)*RPAR(7)*X**(3.d0*RPAR(1))) - 18.d0*H0**4*(1.d0 + RPAR(1))*(-1.d0&
       + RPAR(2))*RPAR(7)*(1.d0 + RPAR(2)*(-1.d0 + X**(3.d0*RPAR(1)))))*Y(2) + 3.d0*H0*Sqrt(X**(-1.d0 - &
       3.d0*RPAR(1))*(1.d0 + RPAR(2)*(-1.d0 + X**(3.d0*RPAR(1)))))*(18.d0*H0**4*(1.d0&
       + RPAR(1))**2*(-1.d0 + RPAR(2))**2*RPAR(7)*Y(4) + RPAR(4)**2*X**(1.d0 + 6.d0*RPAR(1))*((9.d0*H0**2*(1.d0 + &
       RPAR(1))*RPAR(2) + 4.d0*RPAR(4)**2*RPAR(5)*X)*Y(3) + &
       2.d0*(1.d0 + RPAR(1))*(2.d0*H0**2*RPAR(2)*RPAR(7) + RPAR(4)**2*(1.d0 - 3.d0*RPAR(3)**2 + 2.d0*RPAR(6))*X)*Y(4)) - &
       H0**2*(1.d0 + RPAR(1))*(-1.d0 + &
       RPAR(2))*X**(3.d0*RPAR(1))*(18.d0*RPAR(4)**2*RPAR(5)*X*Y(3) + (18.d0*H0**2*(1.d0 + RPAR(1))*RPAR(2)*RPAR(7) + &
       RPAR(4)**2*(9.d0*(1.d0 + RPAR(1))*(1.d0 + 2.d0*RPAR(6)) +&
       4.d0*RPAR(7))*X)*Y(4)))))/(6.d0*H0*(1.d0 + RPAR(1))*RPAR(4)**4*Sqrt((RPAR(2) - (-1.d0 + RPAR(2))/X**(3.d0*RPAR(1)))/X))

end subroutine RHSPER                                          

subroutine JRHSPER(N,X,Y,DFY,LDFY,RPAR,IPAR)

           use functions
           use fiducial

           IMPLICIT REAL*8 (A-H,O-Z)
           Dimension Y(N),DFY(LDFY,N),RPAR(7)

end subroutine JRHSPER

