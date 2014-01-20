module LocalVar

  real(8), allocatable :: T3(:,:)       !Particle temperature array
  real(8), allocatable :: T3S(:)
  real(8), allocatable :: rc(:)       !Coordinates of centeres
  real(8), allocatable :: rce(:)      !Coordinates of East side cells
  real(8), allocatable :: rcw(:)      !Coordinates of West side cells
  real(8), allocatable :: ap(:)       !Coeff at T(m)
  real(8), allocatable :: ce(:)       !Coeff at T(m+1)
  real(8), allocatable :: bw(:)       !Coeff at T(m-1)
  real(8), allocatable :: RHD(:,:)      !Coeff at right hand side
  real(8), allocatable :: uu(:)       !Coeff for Tridiagonal matrix
  real(8), allocatable :: vv(:)       !Coeff for Tridiagonal matrix
  real(8), allocatable :: T3AV(:)
  
  
  real(8) :: TimeCalc, time, tau, dt                    !Время расчёта, текущее время, время подачи расплава, шаг расчёта
  real(8) :: MF_Portion, MF_T, MLiq_T, M_Cold, MLiq     !Масса одной порции расплава, общая масса расплава, масса холодного столба жидкости, масса холодного столба натрия отщипляемого при возникновении новой порции
  real(8) :: dr                                    !cell size
  integer :: mesh
  integer :: PartGroupsTot, PartGroupsCur               !Total Number of Particle Groups, Current Number of Particle Groups
  integer(8) :: n_parts                    !Number of particles in each portion
  integer :: i, j, m, n, k_loop, freq                 !Loop parametres
  real(8) :: CFL
  real(8) :: Mr, Rv, Fi, PsiRo
  real(8) :: BB4, BB5
  real(8) :: VF, VR, VNa, Vrr, dz
  real(8) :: y_film
  integer :: heatmodel
  
  real(8) :: A1, A2, A3, B1, B2, B3, C1, C2, C3, RH1, RH2, RH3
  real(8) :: Znam, Chisl_W, Chisl_P, Chisl_T
  real(8) :: W_NEW, P_NEW, TLiq_New, Z_Old, Vg_Old, Vg, T_Ar_Old
  real(8) :: EKin, ACompr, rKonvers, E_F
  
  
  real(8) :: TF_init, TLiq_init, RoLiq_init, RoG_init, vsLiq, P_init
  real(8) :: L_Init, Z_Init, zetta

 
  
  real(8) :: P, Z, W
  
  real(8) :: Time_Ac 
  real(8) :: R3, QV, Q, MF, dq_part, dq_flux, dQ, dqdt


  real(8), parameter :: Thc_UO2 = 5.6D0
  real(8), parameter :: Thc_ZrO2 = 0.
  real(8), parameter :: Thc_Al2O3 = 8.4D0
  real(8), parameter :: Thc_SS = 30.D0
  real(8), parameter :: Thc_Tin = 54.D0

  real(8), parameter :: Cp_UO2 = 500.D0    !J/(kg*K)
  real(8), parameter :: Cp_ZrO2 = 0.
  real(8), parameter :: Cp_Al2O3 = 2800.D0
  real(8), parameter :: Cp_SS = 750.D0
  real(8), parameter :: Cp_Tin = 241.D0

  real(8), parameter :: Ro_UO2 = 8700.D0  !kg/m3
  real(8), parameter :: Ro_ZrO2 = 0.
  real(8), parameter :: Ro_Al2O3 = 3000.D0       
  real(8), parameter :: Ro_SS = 7000.D0
  real(8), parameter :: Ro_Tin = 6290.D0       

  real(8) :: Thc3, Ro3, Cp3
  real(8) :: ThcF, RoF, CpF

  real(8), parameter :: CpAr = 0.52D3
  real(8), parameter :: Thc_Ar = 42.D-3     !W/m/K
  real(8), parameter :: Visc_Ar = 250.D-7   !N*s/m2
  real(8), parameter :: PI = 4.d0*datan(1.d0)		!PI
  real(8), parameter :: grav = 9.81d0
  real(8), parameter :: k = 1.4                       !isentropic compression of the fission gases = 1.4
  real(8), parameter :: Rg = 8.31441D0        !J/mol/K
  real(8), parameter :: MuAr = 39.948D-3      !kg/mol
  real(8), parameter :: MuVap = 23.D-3        !kg/mol  

!!!!Global


  real(8) :: AlfaLiq, CpLiq, BettaLiq, GammaLiq
  real(8) :: AlfaVAp, CpVap, BettaVap, GammaVap
  
  real(8) :: HLiq, HVap
  real(8) :: RoLiq, RoVap
  real(8) :: vLiq, vVap
  real(8) :: dRoLiqdT, dRoVapdT
  real(8) :: dHLiqdT, dHVapdT
  real(8) :: dPsdT
  
  
!Sodium properties
  real(8) :: Tliq, TSat
  real(8) :: H1, H2, Ro1, Ro2, v1, v2
  real(8) :: dPdT, dRo1dT, dRo2dT, Cp1, Cp2, dH1dT, dH2dT
  real(8) :: Visc1, Visc2, Thc1, Thc2, Sig1, Vs1
  real(8) :: Gr1, Pr1
  
  
  real(8) :: h_t
  real(8) :: q_chf, T_chf, q_min, T_min, n_exp, q_tb
  
  integer :: PHASE
  integer :: DirectContact, Coolant, transient, isOpenTop


  
  
!Specific for Phase2
real(8) :: X, X_init, dXdt                                       ! Void fraction
real(8) :: TVap


!Geometry
real(8) :: Fz, r_ves, dHydro

real(8) :: T_Wall_Init

!Gas parametres

real(8) :: RgMuAr, CvAr, T_Ar, MAr, PExt, Vg_Init, Ro_Ar, Pr_Ar, Re_Ar
real(8) :: P_ext, MG
real(8) :: CoefCT, CT, Alfa_Wall
real(8) :: Znam_T_Ar, Chisl_T_Ar

end module LocalVar

    
! ===========================================================================
    MODULE RH_CONST
! ===========================================================================

   REAL(8),    PARAMETER :: PI=3.141592653589793D0
   REAL(8),    PARAMETER :: GRAVITY=9.80665D0
   REAL(8),    PARAMETER :: P_MIN =613.6D0    
   REAL(8),    PARAMETER :: SIGR=5.67D-8
   !REAL(8),    PARAMETER :: EV=12.7D0
   REAL(8),    PARAMETER :: SIN_45G = 0.7071067811865475244D0
   REAL(8),    PARAMETER :: MVAPOR=0.018D0
   REAL(8),    PARAMETER :: R_UNIV=461.526D0
   REAL(8),    PARAMETER :: T_CRIT=647.096D0


   REAL(8),	PARAMETER :: R_CONST =8.31
   REAL(8),    PARAMETER :: P_MIN_SODIUM =1.80D-4
   REAL(8),    PARAMETER :: P_CRIT_SODIUM=25.64D6
   REAL(8),    PARAMETER :: MVAPOR_SODIUM=0.02299D0
   REAL(8),    PARAMETER :: R_UNIV_SODIUM=277.926
   REAL(8),    PARAMETER :: T_CRIT_SODIUM=2503.7
   REAL(8),    PARAMETER :: lambda_ev_cond=2.0d0
   REAL(8),    parameter :: tc_na = 2503.7D0
   REAL(8),    parameter :: roc_na = 219.D0
   REAL(8),    parameter :: tm_na = 371.0D0
   REAL(8),    parameter :: MU_NA = 22.99D-3 
! kachulin 28.01.10 end

    END MODULE RH_CONST