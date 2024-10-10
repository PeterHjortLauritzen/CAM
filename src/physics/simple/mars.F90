module mars
!------------------------------------------------------------------------------------
!
! Purpose: Implement simple mars physics
!
!====================================================================================
  !
  ! The only modules that are permitted
  !--------------------------------------
  use shr_kind_mod,        only: r8=>shr_kind_r8, cl=>shr_kind_cl
  use shr_const_mod, only: pi => shr_const_pi
  use spmd_utils,    only: masterproc, mpicom
  use spmd_utils,    only: mpicom, mstrid=>masterprocid,   &
                           mpi_real8

  ! Set all Global values and routine to private by default
  ! and then explicitly set their exposure
  !---------------------------------------------------------
  implicit none
  private
  save

  public:: mars_set_const
  public:: mars_condensate_NONE
  public:: mars_condensate_USER
  public:: mars_gw_drag_NONE
  public:: mars_gw_drag_USER
  public:: mars_radiation_exort
  public:: mars_radiation_NONE
  public:: mars_radiation_USER



  ! Global Tuning Parameters:
  !   T0 and E0  are the temperature and saturation vapor pressure used
  !   to calculate qsat values, the saturation value for Q (kg/kg)
  !--------------------------------------------------------------------
  real(r8):: T0
  real(r8):: E0
  real(r8):: Erad
  real(r8):: Wind_min
  real(r8):: Z0
  real(r8):: Ri_c
  real(r8):: Karman
  real(r8):: Fb
  real(r8):: Rs0
  real(r8):: DeltaS
  real(r8):: Tau_eqtr
  real(r8):: Tau_pole
  real(r8):: LinFrac
  real(r8):: Boltz
  real(r8):: C_ocn

  ! Private data
  !----------------------
  real(r8),private :: gravit    ! g: gravitational acceleration (m/s2)
  real(r8),private :: cappa     ! Rd/cp
  real(r8),private :: rair      ! Rd: dry air gas constant (J/K/kg)
  real(r8),private :: cpair     ! cp: specific heat of dry air (J/K/kg)
  real(r8),private :: latvap    ! L: latent heat of vaporization (J/kg)
  real(r8),private :: rh2o      ! Rv: water vapor gas constant (J/K/kg)
  real(r8),private :: epsilo    ! Rd/Rv: ratio of h2o to dry air molecular weights
  real(r8),private :: rhoh2o    ! density of liquid water (kg/m3)
  real(r8),private :: zvir      ! (rh2o/rair) - 1, needed for virtual temperature
  real(r8),private :: ps0       ! Base state surface pressure (Pa)

  real(r8),private :: latvap_div_cpair ! latvap/cpair
  real(r8),private :: latvap_div_rh2o  ! latvap/rh2o

  real(r8),private,allocatable:: etamid(:) ! hybrid coordinate - midpoints


contains
  !=======================================================================
  subroutine mars_set_const(I_gravit,I_cappa   ,I_rair    ,I_cpair  ,I_latvap  , &
                                I_rh2o  ,I_epsilo  ,I_rhoh2o  ,I_zvir   ,I_ps0     , &
                                I_etamid,I_T0      ,I_E0      ,I_Erad   ,I_Wind_min, &
                                I_Z0    ,I_Ri_c    ,I_Karman  ,I_Fb     ,I_Rs0     , &
                                I_DeltaS,I_Tau_eqtr,I_Tau_pole,I_LinFrac,I_Boltz   , &
                                I_Cocn                                               )
    !
    ! mars_set_const: Set parameters and constants for the Mars
    !                     Model fomulation. Optional inputs can be provided
    !                     to over-ride the model defaults.
    !=====================================================================

  use cam_abortutils,      only: handle_allocate_error

    !
    ! Passed Variables
    !-------------------
    real(r8),intent(in):: I_gravit
    real(r8),intent(in):: I_cappa
    real(r8),intent(in):: I_rair
    real(r8),intent(in):: I_cpair
    real(r8),intent(in):: I_latvap
    real(r8),intent(in):: I_rh2o
    real(r8),intent(in):: I_epsilo
    real(r8),intent(in):: I_rhoh2o
    real(r8),intent(in):: I_zvir
    real(r8),intent(in):: I_ps0
    real(r8),intent(in):: I_etamid(:)

    real(r8),intent(in) :: I_T0
    real(r8),intent(in) :: I_E0
    real(r8),intent(in) :: I_Erad
    real(r8),intent(in) :: I_Wind_min
    real(r8),intent(in) :: I_Z0
    real(r8),intent(in) :: I_Ri_c
    real(r8),intent(in) :: I_Karman
    real(r8),intent(in) :: I_Fb
    real(r8),intent(in) :: I_Rs0
    real(r8),intent(in) :: I_DeltaS
    real(r8),intent(in) :: I_Tau_eqtr
    real(r8),intent(in) :: I_Tau_pole
    real(r8),intent(in) :: I_LinFrac
    real(r8),intent(in) :: I_Boltz
    real(r8),intent(in) :: I_Cocn

    integer :: ierr

    ! Set global constants for later use
    !------------------------------------
    gravit   = I_gravit
    cappa    = I_cappa
    rair     = I_rair
    cpair    = I_cpair
    latvap   = I_latvap
    rh2o     = I_rh2o
    epsilo   = I_epsilo
    rhoh2o   = I_rhoh2o
    zvir     = I_zvir
    ps0      = I_ps0
    T0       = I_T0
    E0       = I_E0
    Erad     = I_Erad
    Wind_min = I_Wind_min
    Z0       = I_Z0
    Ri_c     = I_Ri_c
    Karman   = I_Karman
    Fb       = I_Fb
    Rs0      = I_Rs0
    DeltaS   = I_DeltaS
    Tau_eqtr = I_Tau_eqtr
    Tau_pole = I_Tau_pole
    LinFrac  = I_LinFrac
    Boltz    = I_Boltz
    C_ocn    = I_Cocn

    latvap_div_cpair = latvap/cpair
    latvap_div_rh2o  = latvap/rh2o

    ! allocate space and set the level information
    !----------------------------------------------
    allocate(etamid(size(I_etamid)),stat=ierr)
    if (ierr /= 0) then
      call handle_allocate_error(ierr, 'mars_set_const', 'etamid')
    end if

    etamid = I_etamid

  end subroutine mars_set_const
  !=======================================================================


  !=======================================================================
  subroutine mars_condensate_NONE(ncol,pver,pmid,T,qv)
    !
    ! Precip_process: Implement NO large-scale condensation/precipitation
    !=======================================================================
    !
    ! Passed Variables
    !---------------------
    integer ,intent(in)   :: ncol              ! number of columns
    integer ,intent(in)   :: pver              ! number of vertical levels
    real(r8),intent(in)   :: pmid  (ncol,pver) ! mid-point pressure (Pa)
    real(r8),intent(inout):: T     (ncol,pver) ! temperature (K)
    real(r8),intent(inout):: qv    (ncol,pver) ! specific humidity Q (kg/kg)
    !
    ! Local Values
    !-------------


  end subroutine mars_condensate_NONE
  !=======================================================================

  !=======================================================================
  subroutine mars_condensate_USER(ncol,pver,dtime,pmid,pdel,T,qv)
    !
    ! mars_condensate_USER: This routine is a stub which users can use
    !                           to develop and test their own large scale
    !                           condensation scheme
    !=======================================================================
    !
    ! Passed Variables
    !---------------------
    integer ,intent(in)   :: ncol              ! number of columns
    integer ,intent(in)   :: pver              ! number of vertical levels
    real(r8),intent(in)   :: dtime             ! time step (s)
    real(r8),intent(in)   :: pmid  (ncol,pver) ! mid-point pressure (Pa)
    real(r8),intent(in)   :: pdel  (ncol,pver) ! layer thickness (Pa)
    real(r8),intent(inout):: T     (ncol,pver) ! temperature (K)
    real(r8),intent(inout):: qv    (ncol,pver) ! specific humidity Q (kg/kg)

    ! Local Values
    !-------------

  end subroutine mars_condensate_USER
  !=======================================================================


  !=======================================================================
  subroutine mars_gw_drag_NONE( &
       )
    !
    !
    ! mars_gw_drag: TBD
    !
    !==========================================================================
    !
    ! Passed Variables
    !------------------

    !
    ! Local Values
    !---------------


  end subroutine mars_gw_drag_NONE

  !=======================================================================
  subroutine mars_gw_drag_USER(    &
       )
    !
    ! mars_gw_drag_USER: This routine is a stub which users can use
    !                    to develop and test their own GW_DRAG scheme
    !==========================================================================
    !
    ! Passed Variables
    !------------------

    !
    ! Local Values
    !---------------


  end subroutine mars_gw_drag_USER

  !=======================================================================
!  subroutine mars_radiation_exort(ncol,pver,dtime,clat,pint,pmid,  &
!                                Psfc,Tsfc,Qsfc,T,qv,dtdt_rad, &
!                                Fsolar,Fup_s,Fdown_s,Fup_toa,Fdown_toa)
  subroutine mars_radiation_exort()
    !
    ! The exo-rt radiation parameterization
    !
    ! mars_radiation: This is an implementation of the exo-rt radiation
    !                     scheme used in the Mars model.
    !==========================================================================
    !
    ! Passed Variables
    !-------------------

    !
    ! Local Values
    !-------------
  end subroutine mars_radiation_exort
  !=======================================================================

  subroutine mars_radiation_USER(ncol,pver,dtime,clat,pint,pmid,  &
                                     Psfc,Tsfc,Qsfc,T,qv,dtdt_rad, &
                                     Fsolar,Fup_s,Fdown_s,Fup_toa,Fdown_toa)
    !
    ! mars_radiation_USER: This routine is a stub which users can use
    !                          to develop and test their own radiation scheme
    !==========================================================================
    !
    ! Passed Variables
    !-------------------
    integer ,intent(in)   :: ncol                  ! number of columns
    integer ,intent(in)   :: pver                  ! number of vertical levels
    real(r8),intent(in)   :: dtime                 ! time step (s)
    real(r8),intent(in)   :: clat    (ncol)        ! latitude
    real(r8),intent(in)   :: pint    (ncol,pver+1) ! mid-point pressure (Pa)
    real(r8),intent(in)   :: pmid    (ncol,pver)   ! mid-point pressure (Pa)
    real(r8),intent(in)   :: Psfc    (ncol)        ! surface pressure
    real(r8),intent(in)   :: Tsfc    (ncol)        ! surface temperature (K)
    real(r8),intent(in)   :: Qsfc    (ncol)
    real(r8),intent(inout):: T       (ncol,pver)   ! temperature (K)
    real(r8),intent(in)   :: qv      (ncol,pver)   ! Q (kg/kg)
    real(r8),intent(out)  :: dtdt_rad(ncol,pver)   ! temperature tendency in K/s from relaxation
    real(r8),intent(out)  :: Fsolar  (ncol)        !
    real(r8),intent(out)  :: Fup_s   (ncol)        !
    real(r8),intent(out)  :: Fdown_s (ncol)        !
    real(r8),intent(out)  :: Fup_toa  (ncol)       !
    real(r8),intent(out)  :: Fdown_toa(ncol)      !
    !
    ! Local Values
    !-------------

    dtdt_rad = 0._r8
    Fsolar   = 0._r8
    Fup_s    = 0._r8
    Fdown_s  = 0._r8
    Fup_toa    = 0._r8
    Fdown_toa  = 0._r8


  end subroutine mars_radiation_USER

  !=======================================================================

  subroutine mars_radiation_NONE(ncol,pver,dtime,clat,pint,pmid,  &
                                     Psfc,Tsfc,Qsfc,T,qv,dtdt_rad, &
                                     Fsolar,Fup_s,Fdown_s,Fup_toa,Fdown_toa)
    !
    ! mars_radiation_NONE: This routine is a stub for no radiation forcing
    !==========================================================================
    !
    ! Passed Variables
    !-------------------
    integer ,intent(in)   :: ncol                  ! number of columns
    integer ,intent(in)   :: pver                  ! number of vertical levels
    real(r8),intent(in)   :: dtime                 ! time step (s)
    real(r8),intent(in)   :: clat    (ncol)        ! latitude
    real(r8),intent(in)   :: pint    (ncol,pver+1) ! mid-point pressure (Pa)
    real(r8),intent(in)   :: pmid    (ncol,pver)   ! mid-point pressure (Pa)
    real(r8),intent(in)   :: Psfc    (ncol)        ! surface pressure
    real(r8),intent(in)   :: Tsfc    (ncol)        ! surface temperature (K)
    real(r8),intent(in)   :: Qsfc    (ncol)
    real(r8),intent(inout):: T       (ncol,pver)   ! temperature (K)
    real(r8),intent(in)   :: qv      (ncol,pver)   ! Q (kg/kg)
    real(r8),intent(out)  :: dtdt_rad(ncol,pver)   ! temperature tendency in K/s from relaxation
    real(r8),intent(out)  :: Fsolar  (ncol)        !
    real(r8),intent(out)  :: Fup_s   (ncol)        !
    real(r8),intent(out)  :: Fdown_s (ncol)        !
    real(r8),intent(out)  :: Fup_toa  (ncol)       !
    real(r8),intent(out)  :: Fdown_toa(ncol)       !
    !
    ! Local Values
    !-------------

    dtdt_rad = 0._r8
    Fsolar   = 0._r8
    Fup_s    = 0._r8
    Fdown_s  = 0._r8
    Fup_toa    = 0._r8
    Fdown_toa  = 0._r8


  end subroutine mars_radiation_NONE
  !=======================================================================

end module mars
