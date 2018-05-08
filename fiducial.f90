Module fiducial

  Implicit none

  save

  !#####################################
  ! PARAMETERS OF THE COSMOLOGICAL MODEL
  !#####################################


  Real*8,parameter :: kp = 5.d-2
  Real*8,parameter :: Pi = 3.141592653589793d0
!  Real*8,parameter :: Mpc_to_km = 3.1d19
!  Real*8,parameter :: km_to_m = 1.d3
!  Real*8,parameter :: keV_to_eV = 1.d3
  Real*8,parameter :: km_to_inverse_GeV = 1.d3/1.97d-16
  Real*8,parameter :: s_to_inverse_GeV = 1.d0/6.58d-25
  Real*8,parameter :: speedL = 1.0d0
  Real*8,parameter :: speed_of_light = 2.9979d5 ! km/s 
  Real*8,parameter :: omega_b = 2.218d-2
  Real*8,parameter :: omega_cdm = 1.205d-1
!  Real*8,parameter :: me = 5.11d2 ! keV
!  Real*8,parameter :: alpha_G = 1.752d-45 ! Gravitational coupling constant 
  Real*8,parameter :: h_factor = 6.693d-1
  Real*8,parameter :: H0 = 1.d2*h_factor*km_to_inverse_GeV/s_to_inverse_GeV ! in Mpc^{-1}. UNITS FOR WAVENUMBER ARE ALSO Mpc^{-1
!  Real*8,parameter :: critical_density = 1.8788d-26*h_factor**2 
!  Real*8,parameter :: G = 3.d0*H0**2/8.d0/Pi/critical_density !alpha_G/me**2/keV_to_eV**2   ! eV^{-2}
  Real*8,parameter :: ks = 1.065d-2*h_factor    ! Mpc^{-1}
  Real*8,parameter :: lower_limit_ks = 1.d-1
  Real*8,parameter :: upper_limit_ks = 1.d-1
  Real*8,parameter :: wavenumber_k = 4.d2*H0 !5.d2*H0   !5.d-2 !3.d2*H0 !1.5d-2 !1.d3*H0 ! in Mpc^{-1}
  Real*8,parameter :: Omega_m = (omega_b + omega_cdm)/h_factor**2
  Real*8,parameter :: n_s = 9.619d-1
  Real*8,parameter :: A_s = 2.12424d-9
  Real*8,parameter :: m_ncdm = 6.0d-2
  Real*8,parameter :: N_ur = 2.0328d0
  Real*8,parameter :: N_ncdm = 1.d0
  Real*8,parameter :: deg_ncdm = 1.d0
  Real*8,parameter :: tau = 5.96d-2
  Real*8,parameter :: nc_bias_b0 = 1.0d0
  Real*8,parameter :: cs2_fld = 1.d-1 !3.333334d0 ! 1.d0 ; 1.d-4 ; 1.d-6 ; 3.333334d0 ; 3.3334d0 ; 4.3d0 
  Real*8,parameter :: w0_fld = -1.05d0 !-8.0d-1 
  !Real*8,parameter :: wa_fld = 0.d0 
  Real*8,parameter :: e_pi = 0.0d0
  Real*8,parameter :: f_pi = 1.d-1 ! 0.d0 ; 5.d0
  Real*8,parameter :: g_pi = 1.0d-2 ! 1.d0 ; 1.d1 g_pi
  Real*8,parameter :: b_fR = 1.d-1  ! PARAMETER IN THE F(R) PARAMETRISATION BY BASILAKOS ET AL. 

  !######################################
  ! OTHER PARAMETERS REQUIRED BY THE CODE
  !######################################

  Integer*4,parameter    :: DEA_MODEL = 1 ! 1: DEA MODEL ONLY INCLUDING e_pi; 2: DEA MODEL INCLUDING f_pi and g_pi; 3: DEA MODEL INCLUDING e_pi, f_pi, g_pi 
  Integer*4,parameter    :: number_DEA_parameters = DEA_MODEL ! NUMBER OF DEA MODEL PARAMETERS
  Integer*4,parameter    :: number_of_parameters = 7
  Integer*4,parameter    :: dimension_system_ode = 6 
  Integer*4,parameter    :: UNIT_EXE_FILE = 91 
  Integer*4,parameter    :: UNIT_OUTPUT_FILE = 92
  Integer*4,parameter    :: UNIT_OUTPUT_FILE2 = 93
  Integer*4,parameter    :: UNIT_TEST = 94
  Integer*4,parameter    :: UNIT_OUTPUT_DERIVATIVES = 95
  Integer*4,parameter    :: IMAS_RADAU = 0 ! 0: M IS THE IDENTITY 1: M IS GIVEN BY MAS SUBROUTINE

  Real*8,parameter :: initial_scale_factor = 1.d-5
  Real*8,parameter :: final_scale_factor = 1.d0
  Real*8,parameter :: switch_off_pressure_perturbation_terms = 1.d1

  Character(len=*),parameter :: MG_parametrisation = 'HS_Basilakos' !'HS_Basilakos' ! 'HS_Basilakos', 'Starobinsky_Basilakos', 'GR',... 

  !##############
  ! PATH TO FILES
  !##############

  Character(len=*),parameter :: Execution_information = './output/execution_information.txt'
  Character(len=*),parameter :: NUMERICAL_SOLUTION = './output/numerical_solution.txt'
  Character(len=*),parameter :: ANALYTICAL_SOLUTION = './output/analytical_solution.txt'

End Module fiducial
