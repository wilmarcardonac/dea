Module fiducial

  Implicit none

  save

  !#####################################
  ! PARAMETERS OF THE COSMOLOGICAL MODEL
  !#####################################

  Real*8 :: ks=0.01065D0
  Real*8 :: kp=0.05D0
  Real*8 :: As=2.0D-9
  Real*8 :: ns=0.96D0
  Real*8,parameter :: Pi=3.141592653589793D0
  Real*8,parameter :: speedL=1.0D0
  Real*8,parameter :: speed_of_light = 2.9979d5 ! km/s 
  Real*8,parameter :: wavenumber_k = 1.5d-2 ! in Mpc^{-1}

  Real*8,parameter :: omega_b = 2.218d-2
  Real*8,parameter :: omega_cdm = 1.205d-1
  Real*8,parameter :: n_s = 9.619d-1
  Real*8,parameter :: A_s = 2.12424d-9
  Real*8,parameter :: H0 = 6.693d1/speed_of_light ! in Mpc^{-1}. UNITS FOR WAVENUMBER ARE ALSO Mpc^{-1}
  Real*8,parameter :: m_ncdm = 6.0d-2
  Real*8,parameter :: N_ur = 2.0328d0
  Real*8,parameter :: N_ncdm = 1.d0
  Real*8,parameter :: deg_ncdm = 1.d0
  Real*8,parameter :: tau = 5.96d-2
  Real*8,parameter :: nc_bias_b0 = 1.0d0
  Real*8,parameter :: cs2_fld = 3.333334d0 ! 1.d0 ; 1.d-4 ; 1.d-6 ; 3.333334d0 ; 3.3334d0 ; 4.3d0 
  Real*8,parameter :: w0_fld = -8.0d-1 
  !Real*8,parameter :: wa_fld = 0.d0 
  Real*8,parameter :: e_pi = 0.0d0
  Real*8,parameter :: f_pi = 5.0d0 ! 0.d0 ; 5.d0
  Real*8,parameter :: g_pi = 1.0d0 ! 0.d0 ; 1.d0 ! THIS IS ACTUALLY log10 g_pi

  !######################################
  ! OTHER PARAMETERS REQUIRED BY THE CODE
  !######################################

  Integer*4,parameter    :: DEA_MODEL = 1 ! 1: DEA MODEL ONLY INCLUDING e_pi; 2: DEA MODEL INCLUDING f_pi and g_pi; 3: DEA MODEL INCLUDING e_pi, f_pi, g_pi 
  Integer*4,parameter    :: number_DEA_parameters = DEA_MODEL ! NUMBER OF DEA MODEL PARAMETERS
  Integer*4,parameter    :: number_of_parameters = 7
  Integer*4,parameter    :: dimension_system_ode = 4 
  Integer*4,parameter    :: UNIT_EXE_FILE = 91 
  Integer*4,parameter    :: UNIT_OUTPUT_FILE = 92
  Integer*4,parameter    :: UNIT_OUTPUT_FILE2 = 93

  Real*8,parameter :: initial_scale_factor = 1.d-5
  Real*8,parameter :: final_scale_factor = 1.d0

  !##############
  ! PATH TO FILES
  !##############

  Character(len=*),parameter :: Execution_information = './output/execution_information.txt'
  Character(len=*),parameter :: NUMERICAL_SOLUTION = './output/numerical_solution.txt'
  Character(len=*),parameter :: ANALYTICAL_SOLUTION = './output/analytical_solution.txt'

End Module fiducial
