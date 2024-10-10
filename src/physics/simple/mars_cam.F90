module mars_cam
!-----------------------------------------------------------------------
!
! Purpose: Implement Mars idealized forcings
!
!============================================================================
  ! Useful modules
  !-------------------
  use shr_kind_mod,   only: r8 => shr_kind_r8, cl=>shr_kind_cl
  use shr_const_mod,  only: pi => shr_const_pi
  use physconst,      only: gravit, cappa, rair, cpair, latvap, rh2o, epsilo, rhoh2o, zvir,pstd
  use phys_grid,      only: get_rlat_all_p, get_rlon_all_p
  use ppgrid,         only: pcols, pver, pverp, begchunk, endchunk
  use constituents,   only: pcnst
  use physics_buffer, only: dtype_r8, pbuf_add_field, physics_buffer_desc, pbuf_old_tim_idx, &
                            dyn_time_lvls, pbuf_get_field
  use camsrfexch,     only: cam_in_t,cam_out_t
  use cam_history,    only: outfld
  use time_manager,   only: is_first_step, get_nstep, get_step_size, get_curr_calday, get_curr_date
  use cam_abortutils, only: endrun
  use spmd_utils,     only: masterproc
  use cam_logfile,    only: iulog
  use hycoef,         only: ps0, etamid
  use spmd_utils,     only: mpicom, mstrid=>masterprocid, mpi_real8

  use pio             ,only: file_desc_t, var_desc_t, io_desc_t, pio_double, pio_def_var
  use pio             ,only: pio_write_darray, pio_read_darray, pio_inq_varid
  use cam_grid_support,only: cam_grid_id, cam_grid_dimensions, cam_grid_get_decomp
  use shr_const_mod,   only: SHR_CONST_STEBOL, SHR_CONST_REARTH, SHR_CONST_KARMAN, SHR_CONST_TKTRIP
  use cam_control_mod, only: lambm0, obliqr, eccen, mvelpp
  use exoplanet_mod,   only: do_exo_rt_clearsky, exo_rad_step, do_exo_rt_spectral, exo_diurnal, &
                              exo_n2mmr, exo_h2mmr,exo_c2h6mmr,exo_ndays !exo_o2mmr, exo_o3mmr?
  use radgrid,         only: ntot_wavlnrng
  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure.
  !---------------------------------------------------------
  implicit none
  private
  save

  public :: mars_register
  public :: mars_readnl
  public :: mars_init
  public :: mars_condensate_tend
  public :: mars_gw_drag_tend
  public :: mars_restart_init
  public :: mars_restart_write
  public :: mars_restart_read
  public :: mars_radiation_do
!!$  public :: rad_out_t
  private :: mars_surface_init

  ! Gravity Wave Drag Configuatons
  !------------------
  integer,parameter :: GW_DRAG_NONE     = 0           ! Implementation of Mars Gravity Wave drag
  integer,parameter :: GW_DRAG_USER     = 1           ! Optional call for user defined Gravity Wave drag

  ! Tags to identify optional model formulations
  !------------------------------------------------
  integer,parameter :: CONDENSATE_NONE       = 0  ! No Condensation, PRECL=0
  integer,parameter :: CONDENSATE_USER       = 1  ! Optional user defined Condensation scheme

  integer,parameter :: RADIATION_NONE        = 0  ! No Radiation
  integer,parameter :: RADIATION_EXORT       = 1  ! Mars exo-rt
  integer,parameter :: RADIATION_USER        = 2  ! Optional user defined Radiation scheme

  ! Options selecting which CONDENSATE GW_DRAG, RADIATION, etc.. formulations to use.
  !---------------------------------------------------------------------------------
  integer,parameter :: GW_DRAG_OPT      = GW_DRAG_NONE
  integer,parameter :: CONDENSATE_OPT   = CONDENSATE_NONE
  integer,parameter :: RADIATION_OPT    = RADIATION_EXORT

  ! Global Constants
  !---------------------
  real(r8),parameter :: mars_T0     = SHR_CONST_TKTRIP  ! Reference Temperature for E0
  real(r8),parameter :: mars_E0     = 610.78_r8     ! Saturation Vapor pressure @ T0
  real(r8),parameter :: mars_Rs0    = 1360.0_r8     ! Solar Constant
  real(r8),parameter :: mars_Erad   = SHR_CONST_REARTH  ! Earth Radius
  real(r8),parameter :: mars_Karman = SHR_CONST_KARMAN  ! Von Karman constant
  real(r8),parameter :: mars_Boltz  = SHR_CONST_STEBOL  ! Stefan-Boltzmann constant

  ! Global values for Surface Temp, surface fluxes, and radiative heating
  !----------------------------------------------------------------------
  type(var_desc_t)    :: Tsurf_desc      ! Vardesc for restarts
  type(var_desc_t)    :: Qsurf_desc      ! Vardesc for restarts
  real(r8),allocatable :: Tsurf (:,:)     ! Surface Temp
  real(r8),allocatable :: Qsurf (:,:)     ! Surface Q
  real(r8),allocatable :: Fsolar(:,:)     ! Net Solar Heating
  real(r8),allocatable :: Fup   (:,:)     ! Upward Longwave flux
  real(r8),allocatable :: Fdown (:,:)     ! Downward Longwave flux
  real(r8),allocatable :: Fup_toa  (:,:)  ! Upward Longwave flux at TOA
  real(r8),allocatable :: Fdown_toa(:,:)  ! Downward Longwave flux at TOA
  real(r8),allocatable :: SHflux(:,:)     ! Sensible Heat flux
  real(r8),allocatable :: LHflux(:,:)     ! Latent Heat Flux
  real(r8),allocatable :: TUflux(:,:)     ! U stress momentum flux
  real(r8),allocatable :: TVflux(:,:)     ! V stress momentum flux
  real(r8),allocatable :: Cd    (:,:)     ! Surface Drag
  real(r8),allocatable :: clat  (:,:)     ! latitudes(radians) for columns
  real(r8),allocatable :: Fnet  (:,:)     ! Net Radiative Surface Heating
  real(r8),allocatable :: Fnet_toa(:,:)   ! Net Radiative Surface Heating at TOA

  real(r8), parameter :: unset_r8 = huge(1.0_r8)

  ! Global Tuning values
  !------------------------
  real(r8) :: mars_Wind_min   = unset_r8      ! Minimum wind threshold
  real(r8) :: mars_Z0         = unset_r8      ! Roughness Length
  real(r8) :: mars_Ri_c       = unset_r8      ! Crit. Richardson # for stable mixing
  real(r8) :: mars_Fb         = unset_r8      ! Surface layer Fraction
  real(r8) :: mars_Albedo     = unset_r8      ! Mars Albedo
  real(r8) :: mars_DeltaS     = unset_r8      ! Lat variation of shortwave radiation
  real(r8) :: mars_Tau_eqtr   = unset_r8      ! Longwave optical depth at Equator
  real(r8) :: mars_Tau_pole   = unset_r8      ! Longwave optical depth at poles.
  real(r8) :: mars_LinFrac    = unset_r8      ! Stratosphere Linear optical depth param
  real(r8) :: mars_C0         = unset_r8      ! Ocean mixed layer heat capacity
  real(r8) :: mars_WetDryCoef = unset_r8      ! E0 Scale factor to control moisture
  real(r8) :: mars_Tmin       = unset_r8      ! IC: Minimum sst (K)
  real(r8) :: mars_Tdlt       = unset_r8      ! IC: eq-polar difference sst (K)
  real(r8) :: mars_Twidth     = unset_r8      ! IC: Latitudinal width parameter for sst (degrees latitude)


contains
  !==============================================================================
  subroutine mars_register()
    !
    ! mars_register:
    !=====================================================================

  end subroutine mars_register
  !==============================================================================


  !==============================================================================
  subroutine mars_readnl(nlfile)
    !
    ! mars_readnl: Read in parameters controlling Mars parameterizations.
    !=====================================================================
    use namelist_utils,only: find_group_name
    use units         ,only: getunit, freeunit
    ! Passed Variables
    !------------------
    character(len=*),intent(in):: nlfile
    !
    ! Local Values
    !--------------
    integer:: ierr,unitn

    character(len=*), parameter :: sub = 'mars_readnl'

    namelist /mars_nl/ mars_Wind_min, mars_Z0        , mars_Ri_c   , &
                           mars_Fb      , mars_Albedo    , mars_DeltaS , &
                           mars_Tau_eqtr, mars_Tau_pole  , mars_LinFrac, &
                           mars_C0      , mars_WetDryCoef, mars_Tmin   , &
                           mars_Tdlt    , mars_Twidth

   ! Read in namelist values
    !-------------------------
    if(masterproc) then
      unitn = getunit()
      open(unitn,file=trim(nlfile),status='old')
      call find_group_name(unitn,'mars_nl',status=ierr)
      if(ierr == 0) then
        read(unitn,mars_nl,iostat=ierr)
        if(ierr /= 0) then
          call endrun(sub//': ERROR reading namelist')
        endif
      endif
      close(unitn)
      call freeunit(unitn)
    endif

    ! Broadcast namelist values
    !---------------------------
    call mpi_bcast(mars_Wind_min  , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: mars_Wind_min")
    call mpi_bcast(mars_Z0        , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: mars_Z0")
    call mpi_bcast(mars_Ri_c      , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: mars_Ri_c")
    call mpi_bcast(mars_Fb        , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: mars_Fb")
    call mpi_bcast(mars_Albedo    , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: mars_Albedo")
    call mpi_bcast(mars_DeltaS    , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: mars_DeltaS")
    call mpi_bcast(mars_Tau_eqtr  , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: mars_Tau_eqtr")
    call mpi_bcast(mars_Tau_pole  , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: mars_Tau_pole")
    call mpi_bcast(mars_LinFrac   , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: mars_LinFrac")
    call mpi_bcast(mars_C0        , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: mars_C0")
    call mpi_bcast(mars_Tmin      , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: mars_Tmin")
    call mpi_bcast(mars_Tdlt      , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: mars_Tdlt")
    call mpi_bcast(mars_Twidth    , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: mars_Twidth")
    call mpi_bcast(mars_WetDryCoef, 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: mars_WetDryCoef")

  end subroutine mars_readnl
  !==============================================================================


  !==============================================================================
  subroutine mars_init(phys_state,pbuf2d)
    !
    ! mars_init: allocate space for global arrays and initialize values.
    !                Add variables to history outputs
    !=====================================================================
    use physics_types, only: physics_state
    use error_messages,only: alloc_err
    use cam_history,   only: addfld, add_default,horiz_only
    use phys_grid,     only: get_ncols_p, get_rlat_p
    use mars,          only: mars_set_const

    !
    ! Passed Variables
    !------------------
    type(physics_state)      ,pointer:: phys_state(:)
    type(physics_buffer_desc),pointer:: pbuf2d    (:,:)
    !
    ! Local Values
    !---------------
    integer :: istat,lchnk,icol,ncol
    real(r8):: adjusted_E0
    real(r8):: mars_Rs

    ! Initialize constants in Mars module
    !------------------------------------------
    adjusted_E0 = mars_WetDryCoef*mars_E0
    mars_Rs = mars_Rs0*(1._r8-mars_Albedo)
    call mars_set_const(gravit,cappa,rair,cpair,latvap,rh2o,epsilo,rhoh2o,zvir,ps0,etamid,     &
                            mars_T0     ,adjusted_E0    ,mars_Erad    ,mars_Wind_min,  &
                            mars_Z0     ,mars_Ri_c  ,mars_Karman  ,mars_Fb      ,  &
                            mars_Rs     ,mars_DeltaS,mars_Tau_eqtr,mars_Tau_pole,  &
                            mars_LinFrac,mars_Boltz ,mars_C0                           )

    ! Add values for history output
    !---------------------------------
    call addfld('QRL'   ,(/'lev' /),'A','K/s'    ,'Longwave heating rate for gray atmosphere'       )
    call addfld('QRS'   ,(/'lev' /),'A','K/s'    ,'Solar heating rate for gray atmosphere'          )
    call addfld('DTCOND',(/'lev' /),'A','K/s'    ,'T tendency - gray atmosphere moist process'      )
    call addfld('DQCOND',(/'lev' /),'A','kg/kg/s','Q tendency - gray atmosphere moist process'      )
    call addfld('KVH'   ,(/'ilev'/),'A','m2/s'   ,'Vertical diffusion diffusivities (heat/moisture)')
    call addfld('KVM'   ,(/'ilev'/),'A','m2/s'   ,'Vertical diffusion diffusivities (momentum)'     )
    call addfld('VSE'   ,(/'lev' /),'A','K'      ,'VSE: (Tv + gZ/Cp)'                               )
    call addfld('Zm'    ,(/'lev' /),'A','m'      ,'Geopotential height'                             )
    call addfld('Rf'    ,(/'lev' /),'A','1'      ,'Bulk Richardson number (Mars et al 2006, eq 16) / Ri_c' )
    call addfld('DTV'   ,(/'lev' /),'A','K/s'    ,'T tendency due to vertical diffusion'            )
    call addfld('DUV'   ,(/'lev' /),'A','m/s2'   ,'U tendency due to vertical diffusion'            )
    call addfld('DVV'   ,(/'lev' /),'A','m/s2'   ,'V tendency due to vertical diffusion'            )
    call addfld('VD01'  ,(/'lev' /),'A','kg/kg/s','Q tendency (vertical diffusion)'                 )
    call addfld('PRECL' ,horiz_only,'A','m/s'    ,'Large-scale precipitation rate'                  )
    call addfld('PRECC' ,horiz_only,'A','m/s'    ,'Convective precipitation rate'                   )
    call addfld('Tsurf ',horiz_only,'I','K'      ,'Surface Temperature'                             )
    call addfld('Qsurf ',horiz_only,'I','kg/kg'  ,'Surface Water Vapor'                             )
    call addfld('Cdrag' ,horiz_only,'A','1'      ,'Surface Drag Coefficient'                        )
    call addfld('Zgw_drag'  ,horiz_only,'I','m'      ,'GW_DRAG Height'                                      )
    call addfld('SWflux',horiz_only,'I','W/m2'   ,'SW Solar Flux'                                   )
    call addfld('LUflux',horiz_only,'I','W/m2'   ,'LW Upward Radiative Flux at Surface'             )
    call addfld('LDflux',horiz_only,'I','W/m2'   ,'LW Downward Radiative Flux at Surface'           )
    call addfld('LWflux',horiz_only,'I','W/m2'   ,'LW Net Radiative Flux at Surface'                )
    call addfld('LUflux_TOA',horiz_only,'I','W/m2'   ,'LW Upward Radiative Flux at TOA'             )
    call addfld('LDflux_TOA',horiz_only,'I','W/m2'   ,'LW Downward Radiative Flux at TOA'           )
    call addfld('LWflux_TOA',horiz_only,'I','W/m2'   ,'LW Net Radiative Flux at TOA'                )
    call addfld('SHflux',horiz_only,'I','W/m2'   , 'Sensible Heat Flux'        )
    call addfld('LHflux',horiz_only,'I','W/m2'   , 'Latent Heat Flux'          )
    call addfld('TauU'  ,horiz_only,'I','N/m2'   , 'U Surface Stress'          )
    call addfld('TauV'  ,horiz_only,'I','N/m2'   , 'V Surface Stress'          )

    call add_default('QRL'   ,1,' ')
    call add_default('QRS'   ,1,' ')
    call add_default('DTCOND',1,' ')
    call add_default('DQCOND',1,' ')
    call add_default('KVH'   ,1,' ')
    call add_default('KVM'   ,1,' ')
    call add_default('VSE'   ,1,' ')
    call add_default('Zm'    ,1,' ')
    call add_default('Rf'    ,1,' ')
    call add_default('DTV'   ,1,' ')
    call add_default('DUV'   ,1,' ')
    call add_default('DVV'   ,1,' ')
    call add_default('VD01'  ,1,' ')
    call add_default('PRECC' ,1,' ')
    call add_default('PRECL' ,1,' ')
    call add_default('Tsurf' ,1,' ')
    call add_default('Qsurf' ,1,' ')
    call add_default('Cdrag' ,1,' ')
    call add_default('Zgw_drag'  ,1,' ')
    call add_default('SWflux',1,' ')
    call add_default('LUflux',1,' ')
    call add_default('LDflux',1,' ')
    call add_default('LWflux',1,' ')
    call add_default('LUflux_TOA',1,' ')
    call add_default('LDflux_TOA',1,' ')
    call add_default('LWflux_TOA',1,' ')
    call add_default('SHflux',1,' ')
    call add_default('LHflux',1,' ')
    call add_default('TauU'  ,1,' ')
    call add_default('TauV'  ,1,' ')

    call addfld_spectral_intervals()

    ! Allocate Global arrays
    !-------------------------
    allocate(Fsolar(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Mars INIT','Fsolar',pcols*(endchunk-begchunk+1))
    allocate(Fup   (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Mars INIT','Fup'   ,pcols*(endchunk-begchunk+1))
    allocate(Fdown (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Mars INIT','Fdown' ,pcols*(endchunk-begchunk+1))
    allocate(Fup_toa   (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Mars INIT','Fup_toa'   ,pcols*(endchunk-begchunk+1))
    allocate(Fdown_toa (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Mars INIT','Fdown_toa' ,pcols*(endchunk-begchunk+1))
    allocate(Fnet(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Mars INIT','Fnet',pcols*(endchunk-begchunk+1))
    allocate(Fnet_toa(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Mars INIT','Fnet_toa'  ,pcols*(endchunk-begchunk+1))

    allocate(SHflux(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Mars INIT','SHflux',pcols*(endchunk-begchunk+1))
    allocate(LHflux(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Mars INIT','LHflux',pcols*(endchunk-begchunk+1))
    allocate(TUflux(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Mars INIT','TUflux',pcols*(endchunk-begchunk+1))
    allocate(TVflux(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Mars INIT','TVflux',pcols*(endchunk-begchunk+1))
    allocate(Cd    (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Mars INIT','Cd'    ,pcols*(endchunk-begchunk+1))
    allocate(clat  (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Mars INIT','clat'  ,pcols*(endchunk-begchunk+1))

    ! Initialize time indices and latitudes
    !----------------------------------------
    do lchnk = begchunk,endchunk
      ncol = get_ncols_p(lchnk)
      do icol = 1,ncol
        clat(icol,lchnk) = get_rlat_p(lchnk,icol)
      end do
    end do

    ! At first model step, initialize some values
    !-----------------------------------------------
!!$    if(is_first_step()) then
!!$      ! Initialize physics buffer values
!!$      !----------------------------------
!!$
!!$      ! Allocate Surface fields
!!$      !-------------------------
!!$      allocate(Tsurf (pcols,begchunk:endchunk),stat=istat)
!!$      call alloc_err(istat,'Mars INIT','Tsurf' ,pcols*(endchunk-begchunk+1))
!!$      allocate(Qsurf (pcols,begchunk:endchunk)  ,stat=istat)
!!$      call alloc_err(istat,'Mars INIT','Qsurf' ,pcols*(endchunk-begchunk+1))
!!$
!!$      ! Initialize Surface temperatures and Q
!!$      !-----------------------------------------------------------------------
!!$      do lchnk = begchunk,endchunk
!!$        ncol = get_ncols_p(lchnk)
!!$
!!$        ! Set to reference values for initialization
!!$        !------------------------------------------------------------
!!$        phys_state(lchnk)%ps(:ncol) = ps0
!!$
!!$        call mars_surface_init(ncol,         clat(:ncol,lchnk), &
!!$                                 phys_state(lchnk)%ps(:ncol),       &
!!$                                                Tsurf(:ncol,lchnk), &
!!$                                                Qsurf(:ncol,lchnk)  )
!!$      end do
!!$    endif

    ! Initialize radiation and flux values to 0.0
    !---------------------------------------------------------------------------
    do lchnk = begchunk,endchunk
      Fsolar(:,lchnk) = 0._r8
      Fup   (:,lchnk) = 0._r8
      Fdown (:,lchnk) = 0._r8
      Fup_toa   (:,lchnk) = 0._r8
      Fdown_toa (:,lchnk) = 0._r8
      SHflux(:,lchnk) = 0._r8
      LHflux(:,lchnk) = 0._r8
      TUflux(:,lchnk) = 0._r8
      TVflux(:,lchnk) = 0._r8
      Cd    (:,lchnk) = 0._r8
      Fnet  (:,lchnk) = 0._r8
    end do

    ! Informational Output
    !----------------------
    if(masterproc) then
      write(iulog,*) ' '
      write(iulog,*) '-----------------------------------------------------------'
      write(iulog,*) '  MARS MODULE INITIALIZED WITH THE FOLLOWING SETTINGS: '
      write(iulog,*) '-----------------------------------------------------------'
      write(iulog,*) 'MARS: gravit='    , gravit
      write(iulog,*) 'MARS: cappa='     , cappa
      write(iulog,*) 'MARS: rair ='     , rair
      write(iulog,*) 'MARS: cpair='     , cpair
      write(iulog,*) 'MARS: latvap='    , latvap
      write(iulog,*) 'MARS: rh2o='      , rh2o
      write(iulog,*) 'MARS: epsilo='    , epsilo
      write(iulog,*) 'MARS: rhoh2o='    , rhoh2o
      write(iulog,*) 'MARS: zvir='      , zvir
      write(iulog,*) 'MARS: ps0='       , ps0
      write(iulog,*) 'MARS: etamid='    , etamid
      write(iulog,*) 'MARS: T0='        , mars_T0
      write(iulog,*) 'MARS: E0='        , mars_E0
      write(iulog,*) 'MARS: Erad='      , mars_Erad
      write(iulog,*) 'MARS: Wind_min='  , mars_Wind_min
      write(iulog,*) 'MARS: Z0='        , mars_Z0
      write(iulog,*) 'MARS: Ri_c='      , mars_Ri_c
      write(iulog,*) 'MARS: Karman='    , mars_Karman
      write(iulog,*) 'MARS: Fb='        , mars_Fb
      write(iulog,*) 'MARS: Rs0='       , mars_Rs0
      write(iulog,*) 'MARS: Albedo='    , mars_Albedo
      write(iulog,*) 'MARS: Rs='        , mars_Rs
      write(iulog,*) 'MARS: DeltaS='    , mars_DeltaS
      write(iulog,*) 'MARS: Tau_eqtr='  , mars_Tau_eqtr
      write(iulog,*) 'MARS: Tau_pole='  , mars_Tau_pole
      write(iulog,*) 'MARS: LinFrac='   , mars_LinFrac
      write(iulog,*) 'MARS: Boltz='     , mars_Boltz
      write(iulog,*) 'MARS: C0='        , mars_C0
      write(iulog,*) 'MARS: Tmin='      , mars_Tmin
      write(iulog,*) 'MARS: Tdlt='      , mars_Tdlt
      write(iulog,*) 'MARS: Twidth='    , mars_Twidth
      write(iulog,*) 'MARS: WetDryCoef=', mars_WetDryCoef
      write(iulog,*) ' '
    endif
contains
    subroutine addfld_spectral_intervals

!------------------------------------------------------------------------
!
! Purpose:  Add spectral output to master field list
!
!------------------------------------------------------------------------
!
    use cam_history,     only: addfld, horiz_only

    ! shortwave fieldsa
    call addfld ('FUS_int01     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 01')
    call addfld ('FDS_int01     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 01')
    call addfld ('FUSC_int01     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 01')
    call addfld ('FDSC_int01     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 01')
    call addfld ('FUS_toa_int01     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 01, toa')
    call addfld ('FDS_toa_int01     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 01, toa')
    call addfld ('FUSC_toa_int01     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 01, toa')
    call addfld ('FDSC_toa_int01     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 01, toa')

    call addfld ('FUS_int02     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 02')
    call addfld ('FDS_int02     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 02')
    call addfld ('FUSC_int02     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 02')
    call addfld ('FDSC_int02     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 02')
    call addfld ('FUS_toa_int02     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 02, toa')
    call addfld ('FDS_toa_int02     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 02, toa')
    call addfld ('FUSC_toa_int02     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 02, toa')
    call addfld ('FDSC_toa_int02     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 02, toa')

    call addfld ('FUS_int03     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 03')
    call addfld ('FDS_int03     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 03')
    call addfld ('FUSC_int03     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 03')
    call addfld ('FDSC_int03     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 03')
    call addfld ('FUS_toa_int03     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 03, toa')
    call addfld ('FDS_toa_int03     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 03, toa')
    call addfld ('FUSC_toa_int03     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 03, toa')
    call addfld ('FDSC_toa_int03     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 03, toa')

    call addfld ('FUS_int04     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 04')
    call addfld ('FDS_int04     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 04')
    call addfld ('FUSC_int04     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 04')
    call addfld ('FDSC_int04     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 04')
    call addfld ('FUS_toa_int04     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 04, toa')
    call addfld ('FDS_toa_int04     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 04, toa')
    call addfld ('FUSC_toa_int04     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 04, toa')
    call addfld ('FDSC_toa_int04     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 04, toa')

    call addfld ('FUS_int05     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 05')
    call addfld ('FDS_int05     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 05')
    call addfld ('FUSC_int05     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 05')
    call addfld ('FDSC_int05     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 05')
    call addfld ('FUS_toa_int05     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 05, toa')
    call addfld ('FDS_toa_int05     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 05, toa')
    call addfld ('FUSC_toa_int05     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 05, toa')
    call addfld ('FDSC_toa_int05     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 05, toa')

    call addfld ('FUS_int06     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 06')
    call addfld ('FDS_int06     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 06')
    call addfld ('FUSC_int06     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 06')
    call addfld ('FDSC_int06     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 06')
    call addfld ('FUS_toa_int06     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 06, toa')
    call addfld ('FDS_toa_int06     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 06, toa')
    call addfld ('FUSC_toa_int06     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 06, toa')
    call addfld ('FDSC_toa_int06     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 06, toa')

    call addfld ('FUS_int07     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 07')
    call addfld ('FDS_int07     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 07')
    call addfld ('FUSC_int07     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 07')
    call addfld ('FDSC_int07     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 07')
    call addfld ('FUS_toa_int07     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 07, toa')
    call addfld ('FDS_toa_int07     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 07, toa')
    call addfld ('FUSC_toa_int07     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 07, toa')
    call addfld ('FDSC_toa_int07     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 07, toa')

    call addfld ('FUS_int08     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 08')
    call addfld ('FDS_int08     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 08')
    call addfld ('FUSC_int08     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 08')
    call addfld ('FDSC_int08     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 08')
    call addfld ('FUS_toa_int08     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 08, toa')
    call addfld ('FDS_toa_int08     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 08, toa')
    call addfld ('FUSC_toa_int08     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 08, toa')
    call addfld ('FDSC_toa_int08     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 08, toa')

    call addfld ('FUS_int09     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 09')
    call addfld ('FDS_int09     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 09')
    call addfld ('FUSC_int09     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 09')
    call addfld ('FDSC_int09     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 09')
    call addfld ('FUS_toa_int09     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 09, toa')
    call addfld ('FDS_toa_int09     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 09, toa')
    call addfld ('FUSC_toa_int09     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 09, toa')
    call addfld ('FDSC_toa_int09     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 09, toa')

    call addfld ('FUS_int10     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 10')
    call addfld ('FDS_int10     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 10')
    call addfld ('FUSC_int10     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 10')
    call addfld ('FDSC_int10     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 10')
    call addfld ('FUS_toa_int10     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 10, toa')
    call addfld ('FDS_toa_int10     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 10, toa')
    call addfld ('FUSC_toa_int10     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 10, toa')
    call addfld ('FDSC_toa_int10     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 10, toa')

    call addfld ('FUS_int11     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 11')
    call addfld ('FDS_int11     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 11')
    call addfld ('FUSC_int11     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 11')
    call addfld ('FDSC_int11     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 11')
    call addfld ('FUS_toa_int11     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 11, toa')
    call addfld ('FDS_toa_int11     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 11, toa')
    call addfld ('FUSC_toa_int11     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 11, toa')
    call addfld ('FDSC_toa_int11     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 11, toa')

    call addfld ('FUS_int12     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 12')
    call addfld ('FDS_int12     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 12')
    call addfld ('FUSC_int12     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 12')
    call addfld ('FDSC_int12     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 12')
    call addfld ('FUS_toa_int12     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 12, toa')
    call addfld ('FDS_toa_int12     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 12, toa')
    call addfld ('FUSC_toa_int12     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 12, toa')
    call addfld ('FDSC_toa_int12     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 12, toa')

    call addfld ('FUS_int13     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 13')
    call addfld ('FDS_int13     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 13')
    call addfld ('FUSC_int13     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 13')
    call addfld ('FDSC_int13     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 13')
    call addfld ('FUS_toa_int13     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 13, toa')
    call addfld ('FDS_toa_int13     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 13, toa')
    call addfld ('FUSC_toa_int13     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 13, toa')
    call addfld ('FDSC_toa_int13     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 13, toa')

    call addfld ('FUS_int14     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 14')
    call addfld ('FDS_int14     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 14')
    call addfld ('FUSC_int14     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 14')
    call addfld ('FDSC_int14     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 14')
    call addfld ('FUS_toa_int14     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 14, toa')
    call addfld ('FDS_toa_int14     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 14, toa')
    call addfld ('FUSC_toa_int14     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 14, toa')
    call addfld ('FDSC_toa_int14     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 14, toa')

    call addfld ('FUS_int15     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 15')
    call addfld ('FDS_int15     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 15')
    call addfld ('FUSC_int15     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 15')
    call addfld ('FDSC_int15     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 15')
    call addfld ('FUS_toa_int15     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 15, toa')
    call addfld ('FDS_toa_int15     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 15, toa')
    call addfld ('FUSC_toa_int15     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 15, toa')
    call addfld ('FDSC_toa_int15     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 15, toa')

    call addfld ('FUS_int16     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 16')
    call addfld ('FDS_int16     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 16')
    call addfld ('FUSC_int16     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 16')
    call addfld ('FDSC_int16     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 16')
    call addfld ('FUS_toa_int16     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 16, toa')
    call addfld ('FDS_toa_int16     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 16, toa')
    call addfld ('FUSC_toa_int16     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 16, toa')
    call addfld ('FDSC_toa_int16     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 16, toa')

    call addfld ('FUS_int17     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 17')
    call addfld ('FDS_int17     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 17')
    call addfld ('FUSC_int17     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 17')
    call addfld ('FDSC_int17     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 17')
    call addfld ('FUS_toa_int17     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 17, toa')
    call addfld ('FDS_toa_int17     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 17, toa')
    call addfld ('FUSC_toa_int17     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 17, toa')
    call addfld ('FDSC_toa_int17     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 17, toa')

    call addfld ('FUS_int18     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 18')
    call addfld ('FDS_int18     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 18')
    call addfld ('FUSC_int18     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 18')
    call addfld ('FDSC_int18     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 18')
    call addfld ('FUS_toa_int18     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 18, toa')
    call addfld ('FDS_toa_int18     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 18, toa')
    call addfld ('FUSC_toa_int18     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 18, toa')
    call addfld ('FDSC_toa_int18     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 18, toa')

    call addfld ('FUS_int19     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 19')
    call addfld ('FDS_int19     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 19')
    call addfld ('FUSC_int19     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 19')
    call addfld ('FDSC_int19     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 19')
    call addfld ('FUS_toa_int19     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 19, toa')
    call addfld ('FDS_toa_int19     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 19, toa')
    call addfld ('FUSC_toa_int19     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 19, toa')
    call addfld ('FDSC_toa_int19     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 19, toa')

    call addfld ('FUS_int20     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 20')
    call addfld ('FDS_int20     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 20')
    call addfld ('FUSC_int20     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 20')
    call addfld ('FDSC_int20     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 20')
    call addfld ('FUS_toa_int20     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 20, toa')
    call addfld ('FDS_toa_int20     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 20, toa')
    call addfld ('FUSC_toa_int20     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 20, toa')
    call addfld ('FDSC_toa_int20     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 20, toa')

    call addfld ('FUS_int21     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 21')
    call addfld ('FDS_int21     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 21')
    call addfld ('FUSC_int21     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 21')
    call addfld ('FDSC_int21     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 21')
    call addfld ('FUS_toa_int21     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 21, toa')
    call addfld ('FDS_toa_int21     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 21, toa')
    call addfld ('FUSC_toa_int21     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 21, toa')
    call addfld ('FDSC_toa_int21     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 21, toa')

    call addfld ('FUS_int22     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 22')
    call addfld ('FDS_int22     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 22')
    call addfld ('FUSC_int22     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 22')
    call addfld ('FDSC_int22     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 22')
    call addfld ('FUS_toa_int22     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 22, toa')
    call addfld ('FDS_toa_int22     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 22, toa')
    call addfld ('FUSC_toa_int22     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 22, toa')
    call addfld ('FDSC_toa_int22     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 22, toa')

    call addfld ('FUS_int23     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 23')
    call addfld ('FDS_int23     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 23')
    call addfld ('FUSC_int23     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 23')
    call addfld ('FDSC_int23     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 23')
    call addfld ('FUS_toa_int23     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 23, toa')
    call addfld ('FDS_toa_int23     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 23, toa')
    call addfld ('FUSC_toa_int23     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 23, toa')
    call addfld ('FDSC_toa_int23     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 23, toa')

    call addfld ('FUS_int24     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 24')
    call addfld ('FDS_int24     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 24')
    call addfld ('FUSC_int24     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 24')
    call addfld ('FDSC_int24     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 24')
    call addfld ('FUS_toa_int24     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 24, toa')
    call addfld ('FDS_toa_int24     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 24, toa')
    call addfld ('FUSC_toa_int24     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 24, toa')
    call addfld ('FDSC_toa_int24     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 24, toa')

    call addfld ('FUS_int25     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 25')
    call addfld ('FDS_int25     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 25')
    call addfld ('FUSC_int25     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 25')
    call addfld ('FDSC_int25     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 25')
    call addfld ('FUS_toa_int25     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 25, toa')
    call addfld ('FDS_toa_int25     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 25, toa')
    call addfld ('FUSC_toa_int25     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 25, toa')
    call addfld ('FDSC_toa_int25     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 25, toa')

    call addfld ('FUS_int26     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 26')
    call addfld ('FDS_int26     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 26')
    call addfld ('FUSC_int26     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 26')
    call addfld ('FDSC_int26     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 26')
    call addfld ('FUS_toa_int26     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 26, toa')
    call addfld ('FDS_toa_int26     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 26, toa')
    call addfld ('FUSC_toa_int26     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 26, toa')
    call addfld ('FDSC_toa_int26     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 26, toa')

    call addfld ('FUS_int27     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 27')
    call addfld ('FDS_int27     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 27')
    call addfld ('FUSC_int27     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 27')
    call addfld ('FDSC_int27     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 27')
    call addfld ('FUS_toa_int27     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 27, toa')
    call addfld ('FDS_toa_int27     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 27, toa')
    call addfld ('FUSC_toa_int27     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 27, toa')
    call addfld ('FDSC_toa_int27     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 27, toa')

    call addfld ('FUS_int28     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 28')
    call addfld ('FDS_int28     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 28')
    call addfld ('FUSC_int28     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 28')
    call addfld ('FDSC_int28     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 28')
    call addfld ('FUS_toa_int28     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 28, toa')
    call addfld ('FDS_toa_int28     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 28, toa')
    call addfld ('FUSC_toa_int28     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 28, toa')
    call addfld ('FDSC_toa_int28     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 28, toa')

    call addfld ('FUS_int29     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 29')
    call addfld ('FDS_int29     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 29')
    call addfld ('FUSC_int29     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 29')
    call addfld ('FDSC_int29     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 29')
    call addfld ('FUS_toa_int29     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 29, toa')
    call addfld ('FDS_toa_int29     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 29, toa')
    call addfld ('FUSC_toa_int29     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 29, toa')
    call addfld ('FDSC_toa_int29     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 29, toa')

    call addfld ('FUS_int30     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 30')
    call addfld ('FDS_int30     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 30')
    call addfld ('FUSC_int30     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 30')
    call addfld ('FDSC_int30     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 30')
    call addfld ('FUS_toa_int30     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 30, toa')
    call addfld ('FDS_toa_int30     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 30, toa')
    call addfld ('FUSC_toa_int30     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 30, toa')
    call addfld ('FDSC_toa_int30     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 30, toa')

    call addfld ('FUS_int31     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 31')
    call addfld ('FDS_int31     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 31')
    call addfld ('FUSC_int31     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 31')
    call addfld ('FDSC_int31     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 31')
    call addfld ('FUS_toa_int31     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 31, toa')
    call addfld ('FDS_toa_int31     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 31, toa')
    call addfld ('FUSC_toa_int31     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 31, toa')
    call addfld ('FDSC_toa_int31     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 31, toa')

    call addfld ('FUS_int32     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 32')
    call addfld ('FDS_int32     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 32')
    call addfld ('FUSC_int32     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 32')
    call addfld ('FDSC_int32     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 32')
    call addfld ('FUS_toa_int32     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 32, toa')
    call addfld ('FDS_toa_int32     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 32, toa')
    call addfld ('FUSC_toa_int32     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 32, toa')
    call addfld ('FDSC_toa_int32     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 32, toa')

    call addfld ('FUS_int33     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 33')
    call addfld ('FDS_int33     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 33')
    call addfld ('FUSC_int33     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 33')
    call addfld ('FDSC_int33     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 33')
    call addfld ('FUS_toa_int33     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 33, toa')
    call addfld ('FDS_toa_int33     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 33, toa')
    call addfld ('FUSC_toa_int33     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 33, toa')
    call addfld ('FDSC_toa_int33     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 33, toa')

    call addfld ('FUS_int34     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 34')
    call addfld ('FDS_int34     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 34')
    call addfld ('FUSC_int34     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 34')
    call addfld ('FDSC_int34     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 34')
    call addfld ('FUS_toa_int34     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 34, toa')
    call addfld ('FDS_toa_int34     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 34, toa')
    call addfld ('FUSC_toa_int34     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 34, toa')
    call addfld ('FDSC_toa_int34     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 34, toa')

    call addfld ('FUS_int35     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 35')
    call addfld ('FDS_int35     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 35')
    call addfld ('FUSC_int35     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 35')
    call addfld ('FDSC_int35     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 35')
    call addfld ('FUS_toa_int35     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 35, toa')
    call addfld ('FDS_toa_int35     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 35, toa')
    call addfld ('FUSC_toa_int35     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 35, toa')
    call addfld ('FDSC_toa_int35     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 35, toa')

    call addfld ('FUS_int36     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 36')
    call addfld ('FDS_int36     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 36')
    call addfld ('FUSC_int36     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 36')
    call addfld ('FDSC_int36     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 36')
    call addfld ('FUS_toa_int36     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 36, toa')
    call addfld ('FDS_toa_int36     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 36, toa')
    call addfld ('FUSC_toa_int36     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 36, toa')
    call addfld ('FDSC_toa_int36     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 36, toa')

    call addfld ('FUS_int37     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 37')
    call addfld ('FDS_int37     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 37')
    call addfld ('FUSC_int37     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 37')
    call addfld ('FDSC_int37     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 37')
    call addfld ('FUS_toa_int37     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 37, toa')
    call addfld ('FDS_toa_int37     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 37, toa')
    call addfld ('FUSC_toa_int37     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 37, toa')
    call addfld ('FDSC_toa_int37     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 37, toa')

    call addfld ('FUS_int38     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 38')
    call addfld ('FDS_int38     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 38')
    call addfld ('FUSC_int38     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 38')
    call addfld ('FDSC_int38     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 38')
    call addfld ('FUS_toa_int38     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 38, toa')
    call addfld ('FDS_toa_int38     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 38, toa')
    call addfld ('FUSC_toa_int38     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 38, toa')
    call addfld ('FDSC_toa_int38     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 38, toa')

    call addfld ('FUS_int39     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 39')
    call addfld ('FDS_int39     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 39')
    call addfld ('FUSC_int39     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 39')
    call addfld ('FDSC_int39     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 39')
    call addfld ('FUS_toa_int39     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 39, toa')
    call addfld ('FDS_toa_int39     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 39, toa')
    call addfld ('FUSC_toa_int39     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 39, toa')
    call addfld ('FDSC_toa_int39     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 39, toa')

    call addfld ('FUS_int40     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 40')
    call addfld ('FDS_int40     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 40')
    call addfld ('FUSC_int40     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 40')
    call addfld ('FDSC_int40     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 40')
    call addfld ('FUS_toa_int40     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 40, toa')
    call addfld ('FDS_toa_int40     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 40, toa')
    call addfld ('FUSC_toa_int40     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 40, toa')
    call addfld ('FDSC_toa_int40     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 40, toa')

    call addfld ('FUS_int41     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 41')
    call addfld ('FDS_int41     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 41')
    call addfld ('FUSC_int41     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 41')
    call addfld ('FDSC_int41     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 41')
    call addfld ('FUS_toa_int41     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 41, toa')
    call addfld ('FDS_toa_int41     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 41, toa')
    call addfld ('FUSC_toa_int41     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 41, toa')
    call addfld ('FDSC_toa_int41     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 41, toa')

    call addfld ('FUS_int42     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 42')
    call addfld ('FDS_int42     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 42')
    call addfld ('FUSC_int42     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 42')
    call addfld ('FDSC_int42     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 42')
    call addfld ('FUS_toa_int42     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 42, toa')
    call addfld ('FDS_toa_int42     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 42, toa')
    call addfld ('FUSC_toa_int42     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 42, toa')
    call addfld ('FDSC_toa_int42     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 42, toa')

    call addfld ('FUS_int43     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 43')
    call addfld ('FDS_int43     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 43')
    call addfld ('FUSC_int43     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 43')
    call addfld ('FDSC_int43     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 43')
    call addfld ('FUS_toa_int43     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 43, toa')
    call addfld ('FDS_toa_int43     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 43, toa')
    call addfld ('FUSC_toa_int43     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 43, toa')
    call addfld ('FDSC_toa_int43     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 43, toa')

    call addfld ('FUS_int44     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 44')
    call addfld ('FDS_int44     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 44')
    call addfld ('FUSC_int44     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 44')
    call addfld ('FDSC_int44     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 44')
    call addfld ('FUS_toa_int44     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 44, toa')
    call addfld ('FDS_toa_int44     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 44, toa')
    call addfld ('FUSC_toa_int44     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 44, toa')
    call addfld ('FDSC_toa_int44     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 44, toa')

    call addfld ('FUS_int45     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 45')
    call addfld ('FDS_int45     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 45')
    call addfld ('FUSC_int45     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 45')
    call addfld ('FDSC_int45     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 45')
    call addfld ('FUS_toa_int45     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 45, toa')
    call addfld ('FDS_toa_int45     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 45, toa')
    call addfld ('FUSC_toa_int45     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 45, toa')
    call addfld ('FDSC_toa_int45     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 45, toa')

    call addfld ('FUS_int46     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 46')
    call addfld ('FDS_int46     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 46')
    call addfld ('FUSC_int46     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 46')
    call addfld ('FDSC_int46     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 46')
    call addfld ('FUS_toa_int46     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 46, toa')
    call addfld ('FDS_toa_int46     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 46, toa')
    call addfld ('FUSC_toa_int46     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 46, toa')
    call addfld ('FDSC_toa_int46     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 46, toa')

    call addfld ('FUS_int47     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 47')
    call addfld ('FDS_int47     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 47')
    call addfld ('FUSC_int47     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 47')
    call addfld ('FDSC_int47     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 47')
    call addfld ('FUS_toa_int47     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 47, toa')
    call addfld ('FDS_toa_int47     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 47, toa')
    call addfld ('FUSC_toa_int47     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 47, toa')
    call addfld ('FDSC_toa_int47     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 47, toa')

    call addfld ('FUS_int48     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 48')
    call addfld ('FDS_int48     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 48')
    call addfld ('FUSC_int48     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 48')
    call addfld ('FDSC_int48     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 48')
    call addfld ('FUS_toa_int48     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 48, toa')
    call addfld ('FDS_toa_int48     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 48, toa')
    call addfld ('FUSC_toa_int48     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 48, toa')
    call addfld ('FDSC_toa_int48     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 48, toa')

    call addfld ('FUS_int49     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 49')
    call addfld ('FDS_int49     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 49')
    call addfld ('FUSC_int49     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 49')
    call addfld ('FDSC_int49     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 49')
    call addfld ('FUS_toa_int49     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 49, toa')
    call addfld ('FDS_toa_int49     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 49, toa')
    call addfld ('FUSC_toa_int49     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 49, toa')
    call addfld ('FDSC_toa_int49     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 49, toa')

    call addfld ('FUS_int50     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 50')
    call addfld ('FDS_int50     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 50')
    call addfld ('FUSC_int50     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 50')
    call addfld ('FDSC_int50     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 50')
    call addfld ('FUS_toa_int50     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 50, toa')
    call addfld ('FDS_toa_int50     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 50, toa')
    call addfld ('FUSC_toa_int50     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 50, toa')
    call addfld ('FDSC_toa_int50     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 50, toa')

    call addfld ('FUS_int51     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 51')
    call addfld ('FDS_int51     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 51')
    call addfld ('FUSC_int51     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 51')
    call addfld ('FDSC_int51     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 51')
    call addfld ('FUS_toa_int51     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 51, toa')
    call addfld ('FDS_toa_int51     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 51, toa')
    call addfld ('FUSC_toa_int51     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 51, toa')
    call addfld ('FDSC_toa_int51     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 51, toa')

    call addfld ('FUS_int52     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 52')
    call addfld ('FDS_int52     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 52')
    call addfld ('FUSC_int52     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 52')
    call addfld ('FDSC_int52     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 52')
    call addfld ('FUS_toa_int52     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 52, toa')
    call addfld ('FDS_toa_int52     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 52, toa')
    call addfld ('FUSC_toa_int52     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 52, toa')
    call addfld ('FDSC_toa_int52     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 52, toa')

    call addfld ('FUS_int53     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 53')
    call addfld ('FDS_int53     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 53')
    call addfld ('FUSC_int53     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 53')
    call addfld ('FDSC_int53     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 53')
    call addfld ('FUS_toa_int53     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 53, toa')
    call addfld ('FDS_toa_int53     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 53, toa')
    call addfld ('FUSC_toa_int53     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 53, toa')
    call addfld ('FDSC_toa_int53     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 53, toa')

    call addfld ('FUS_int54     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 54')
    call addfld ('FDS_int54     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 54')
    call addfld ('FUSC_int54     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 54')
    call addfld ('FDSC_int54     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 54')
    call addfld ('FUS_toa_int54     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 54, toa')
    call addfld ('FDS_toa_int54     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 54, toa')
    call addfld ('FUSC_toa_int54     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 54, toa')
    call addfld ('FDSC_toa_int54     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 54, toa')

    call addfld ('FUS_int55     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 55')
    call addfld ('FDS_int55     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 55')
    call addfld ('FUSC_int55     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 55')
    call addfld ('FDSC_int55     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 55')
    call addfld ('FUS_toa_int55     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 55, toa')
    call addfld ('FDS_toa_int55     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 55, toa')
    call addfld ('FUSC_toa_int55     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 55, toa')
    call addfld ('FDSC_toa_int55     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 55, toa')

    call addfld ('FUS_int56     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 56')
    call addfld ('FDS_int56     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 56')
    call addfld ('FUSC_int56     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 56')
    call addfld ('FDSC_int56     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 56')
    call addfld ('FUS_toa_int56     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 56, toa')
    call addfld ('FDS_toa_int56     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 56, toa')
    call addfld ('FUSC_toa_int56     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 56, toa')
    call addfld ('FDSC_toa_int56     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 56, toa')

    call addfld ('FUS_int57     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 57')
    call addfld ('FDS_int57     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 57')
    call addfld ('FUSC_int57     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 57')
    call addfld ('FDSC_int57     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 57')
    call addfld ('FUS_toa_int57     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 57, toa')
    call addfld ('FDS_toa_int57     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 57, toa')
    call addfld ('FUSC_toa_int57     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 57, toa')
    call addfld ('FDSC_toa_int57     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 57, toa')

    call addfld ('FUS_int58     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 58')
    call addfld ('FDS_int58     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 58')
    call addfld ('FUSC_int58     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 58')
    call addfld ('FDSC_int58     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 58')
    call addfld ('FUS_toa_int58     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 58, toa')
    call addfld ('FDS_toa_int58     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 58, toa')
    call addfld ('FUSC_toa_int58     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 58, toa')
    call addfld ('FDSC_toa_int58     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 58, toa')

    call addfld ('FUS_int59     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 59')
    call addfld ('FDS_int59     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 59')
    call addfld ('FUSC_int59     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 59')
    call addfld ('FDSC_int59     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 59')
    call addfld ('FUS_toa_int59     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 59, toa')
    call addfld ('FDS_toa_int59     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 59, toa')
    call addfld ('FUSC_toa_int59     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 59, toa')
    call addfld ('FDSC_toa_int59     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 59, toa')

    call addfld ('FUS_int60     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 60')
    call addfld ('FDS_int60     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 60')
    call addfld ('FUSC_int60     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 60')
    call addfld ('FDSC_int60     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 60')
    call addfld ('FUS_toa_int60     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 60, toa')
    call addfld ('FDS_toa_int60     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 60, toa')
    call addfld ('FUSC_toa_int60     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 60, toa')
    call addfld ('FDSC_toa_int60     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 60, toa')

    call addfld ('FUS_int61     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 61')
    call addfld ('FDS_int61     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 61')
    call addfld ('FUSC_int61     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 61')
    call addfld ('FDSC_int61     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 61')
    call addfld ('FUS_toa_int61     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 61, toa')
    call addfld ('FDS_toa_int61     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 61, toa')
    call addfld ('FUSC_toa_int61     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 61, toa')
    call addfld ('FDSC_toa_int61     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 61, toa')

    call addfld ('FUS_int62     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 62')
    call addfld ('FDS_int62     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 62')
    call addfld ('FUSC_int62     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 62')
    call addfld ('FDSC_int62     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 62')
    call addfld ('FUS_toa_int62     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 62, toa')
    call addfld ('FDS_toa_int62     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 62, toa')
    call addfld ('FUSC_toa_int62     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 62, toa')
    call addfld ('FDSC_toa_int62     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 62, toa')

    call addfld ('FUS_int63     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 63')
    call addfld ('FDS_int63     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 63')
    call addfld ('FUSC_int63     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 63')
    call addfld ('FDSC_int63     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 63')
    call addfld ('FUS_toa_int63     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 63, toa')
    call addfld ('FDS_toa_int63     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 63, toa')
    call addfld ('FUSC_toa_int63     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 63, toa')
    call addfld ('FDSC_toa_int63     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 63, toa')

    call addfld ('FUS_int64     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 64')
    call addfld ('FDS_int64     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 64')
    call addfld ('FUSC_int64     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 64')
    call addfld ('FDSC_int64     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 64')
    call addfld ('FUS_toa_int64     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 64, toa')
    call addfld ('FDS_toa_int64     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 64, toa')
    call addfld ('FUSC_toa_int64     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 64, toa')
    call addfld ('FDSC_toa_int64     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 64, toa')

    call addfld ('FUS_int65     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 65')
    call addfld ('FDS_int65     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 65')
    call addfld ('FUSC_int65     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 65')
    call addfld ('FDSC_int65     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 65')
    call addfld ('FUS_toa_int65     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 65, toa')
    call addfld ('FDS_toa_int65     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 65, toa')
    call addfld ('FUSC_toa_int65     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 65, toa')
    call addfld ('FDSC_toa_int65     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 65, toa')

    call addfld ('FUS_int66     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 66')
    call addfld ('FDS_int66     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 66')
    call addfld ('FUSC_int66     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 66')
    call addfld ('FDSC_int66     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 66')
    call addfld ('FUS_toa_int66     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 65, toa')
    call addfld ('FDS_toa_int66     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 65, toa')
    call addfld ('FUSC_toa_int66     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 66, toa')
    call addfld ('FDSC_toa_int66     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 66, toa')

    call addfld ('FUS_int67     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 67')
    call addfld ('FDS_int67     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 67')
    call addfld ('FUSC_int67     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 67')
    call addfld ('FDSC_int67     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 67')
    call addfld ('FUS_toa_int67     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 67, toa')
    call addfld ('FDS_toa_int67     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 67, toa')
    call addfld ('FUSC_toa_int67     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 67, toa')
    call addfld ('FDSC_toa_int67     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 67, toa')

    call addfld ('FUS_int68     ',(/'ilev'/),'A','W/m2    ','Shortwave upward flux, spectral interval 68')
    call addfld ('FDS_int68     ',(/'ilev'/),'A','W/m2    ','Shortwave downward flux, spectral interval 68')
    call addfld ('FUSC_int68     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 68')
    call addfld ('FDSC_int68     ',(/'ilev'/),'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 68')
    call addfld ('FUS_toa_int68     ',horiz_only,'A','W/m2    ','Shortwave upward flux, spectral interval 68, toa')
    call addfld ('FDS_toa_int68     ',horiz_only,'A','W/m2    ','Shortwave downward flux, spectral interval 68, toa')
    call addfld ('FUSC_toa_int68     ',horiz_only,'A','W/m2    ','Shortwave clear-sky upward flux, spectral interval 68, toa')
    call addfld ('FDSC_toa_int68     ',horiz_only,'A','W/m2    ','Shortwave clear-sky downward flux, spectral interval 68, toa')

    ! long wave fields
    call addfld ('FUL_int01     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 01')
    call addfld ('FDL_int01     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 01')
    call addfld ('FULC_int01    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 01')
    call addfld ('FDLC_int01    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 01')
    call addfld ('FUL_toa_int01     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 01, toa')
    call addfld ('FDL_toa_int01     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 01, toa')
    call addfld ('FULC_toa_int01    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 01, toa')
    call addfld ('FDLC_toa_int01    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 01, toa')

    call addfld ('FUL_int02     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 02')
    call addfld ('FDL_int02     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 02')
    call addfld ('FULC_int02    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 02')
    call addfld ('FDLC_int02    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 02')
    call addfld ('FUL_toa_int02     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 02, toa')
    call addfld ('FDL_toa_int02     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 02, toa')
    call addfld ('FULC_toa_int02    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 02, toa')
    call addfld ('FDLC_toa_int02    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 02, toa')

    call addfld ('FUL_int03     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 03')
    call addfld ('FDL_int03     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 03')
    call addfld ('FULC_int03    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 03')
    call addfld ('FDLC_int03    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 03')
    call addfld ('FUL_toa_int03     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 03, toa')
    call addfld ('FDL_toa_int03     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 03, toa')
    call addfld ('FULC_toa_int03    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 03, toa')
    call addfld ('FDLC_toa_int03    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 03, toa')

    call addfld ('FUL_int04     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 04')
    call addfld ('FDL_int04     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 04')
    call addfld ('FULC_int04    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 04')
    call addfld ('FDLC_int04    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 04')
    call addfld ('FUL_toa_int04     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 04, toa')
    call addfld ('FDL_toa_int04     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 04, toa')
    call addfld ('FULC_toa_int04    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 04, toa')
    call addfld ('FDLC_toa_int04    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 04, toa')

    call addfld ('FUL_int05     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 05')
    call addfld ('FDL_int05     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 05')
    call addfld ('FULC_int05    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 05')
    call addfld ('FDLC_int05    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 05')
    call addfld ('FUL_toa_int05     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 05, toa')
    call addfld ('FDL_toa_int05     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 05, toa')
    call addfld ('FULC_toa_int05    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 05, toa')
    call addfld ('FDLC_toa_int05    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 05, toa')

    call addfld ('FUL_int06     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 06')
    call addfld ('FDL_int06     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 06')
    call addfld ('FULC_int06    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 06')
    call addfld ('FDLC_int06    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 06')
    call addfld ('FUL_toa_int06     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 06, toa')
    call addfld ('FDL_toa_int06     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 06, toa')
    call addfld ('FULC_toa_int06    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 06, toa')
    call addfld ('FDLC_toa_int06    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 06, toa')

    call addfld ('FUL_int07     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 07')
    call addfld ('FDL_int07     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 07')
    call addfld ('FULC_int07    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 07')
    call addfld ('FDLC_int07    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 07')
    call addfld ('FUL_toa_int07     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 07, toa')
    call addfld ('FDL_toa_int07     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 07, toa')
    call addfld ('FULC_toa_int07    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 07, toa')
    call addfld ('FDLC_toa_int07    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 07, toa')

    call addfld ('FUL_int08     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 08')
    call addfld ('FDL_int08     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 08')
    call addfld ('FULC_int08    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 08')
    call addfld ('FDLC_int08    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 08')
    call addfld ('FUL_toa_int08     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 08, toa')
    call addfld ('FDL_toa_int08     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 08, toa')
    call addfld ('FULC_toa_int08    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 08, toa')
    call addfld ('FDLC_toa_int08    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 08, toa')

    call addfld ('FUL_int09     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 09')
    call addfld ('FDL_int09     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 09')
    call addfld ('FULC_int09    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 09')
    call addfld ('FDLC_int09    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 09')
    call addfld ('FUL_toa_int09     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 09, toa')
    call addfld ('FDL_toa_int09     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 09, toa')
    call addfld ('FULC_toa_int09    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 09, toa')
    call addfld ('FDLC_toa_int09    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 09, toa')

    call addfld ('FUL_int10     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 10')
    call addfld ('FDL_int10     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 10')
    call addfld ('FULC_int10    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 10')
    call addfld ('FDLC_int10    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 10')
    call addfld ('FUL_toa_int10     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 10, toa')
    call addfld ('FDL_toa_int10     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 10, toa')
    call addfld ('FULC_toa_int10    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 10, toa')
    call addfld ('FDLC_toa_int10    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 10, toa')

    call addfld ('FUL_int11     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 11')
    call addfld ('FDL_int11     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 11')
    call addfld ('FULC_int11    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 11')
    call addfld ('FDLC_int11    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 11')
    call addfld ('FUL_toa_int11     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 11, toa')
    call addfld ('FDL_toa_int11     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 11, toa')
    call addfld ('FULC_toa_int11    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 11, toa')
    call addfld ('FDLC_toa_int11    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 11, toa')

    call addfld ('FUL_int12     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 12')
    call addfld ('FDL_int12     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 12')
    call addfld ('FULC_int12    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 12')
    call addfld ('FDLC_int12    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 12')
    call addfld ('FUL_toa_int12     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 12, toa')
    call addfld ('FDL_toa_int12     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 12, toa')
    call addfld ('FULC_toa_int12    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 12, toa')
    call addfld ('FDLC_toa_int12    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 12, toa')

    call addfld ('FUL_int13     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 13')
    call addfld ('FDL_int13     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 13')
    call addfld ('FULC_int13    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 13')
    call addfld ('FDLC_int13    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 13')
    call addfld ('FUL_toa_int13     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 13, toa')
    call addfld ('FDL_toa_int13     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 13, toa')
    call addfld ('FULC_toa_int13    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 13, toa')
    call addfld ('FDLC_toa_int13    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 13, toa')

    call addfld ('FUL_int14     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 14')
    call addfld ('FDL_int14     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 14')
    call addfld ('FULC_int14    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 14')
    call addfld ('FDLC_int14    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 14')
    call addfld ('FUL_toa_int14     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 13, toa')
    call addfld ('FDL_toa_int14     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 13, toa')
    call addfld ('FULC_toa_int14    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 13, toa')
    call addfld ('FDLC_toa_int14    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 13, toa')

    call addfld ('FUL_int15     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 15')
    call addfld ('FDL_int15     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 15')
    call addfld ('FULC_int15    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 15')
    call addfld ('FDLC_int15    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 15')
    call addfld ('FUL_toa_int15     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 15, toa')
    call addfld ('FDL_toa_int15     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 15, toa')
    call addfld ('FULC_toa_int15    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 15, toa')
    call addfld ('FDLC_toa_int15    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 15, toa')

    call addfld ('FUL_int16     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 16')
    call addfld ('FDL_int16     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 16')
    call addfld ('FULC_int16    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 16')
    call addfld ('FDLC_int16    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 16')
    call addfld ('FUL_toa_int16     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 16, toa')
    call addfld ('FDL_toa_int16     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 16, toa')
    call addfld ('FULC_toa_int16    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 16, toa')
    call addfld ('FDLC_toa_int16    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 16, toa')

    call addfld ('FUL_int17     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 17')
    call addfld ('FDL_int17     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 17')
    call addfld ('FULC_int17    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 17')
    call addfld ('FDLC_int17    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 17')
    call addfld ('FUL_toa_int17     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 17, toa')
    call addfld ('FDL_toa_int17     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 17, toa')
    call addfld ('FULC_toa_int17    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 17, toa')
    call addfld ('FDLC_toa_int17    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 17, toa')

    call addfld ('FUL_int18     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 18')
    call addfld ('FDL_int18     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 18')
    call addfld ('FULC_int18    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 18')
    call addfld ('FDLC_int18    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 18')
    call addfld ('FUL_toa_int18     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 18, toa')
    call addfld ('FDL_toa_int18     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 18, toa')
    call addfld ('FULC_toa_int18    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 18, toa')
    call addfld ('FDLC_toa_int18    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 18, toa')

    call addfld ('FUL_int19     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 19')
    call addfld ('FDL_int19     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 19')
    call addfld ('FULC_int19    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 19')
    call addfld ('FDLC_int19    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 19')
    call addfld ('FUL_toa_int19     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 19, toa')
    call addfld ('FDL_toa_int19     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 19, toa')
    call addfld ('FULC_toa_int19    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 19, toa')
    call addfld ('FDLC_toa_int19    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 19, toa')

    call addfld ('FUL_int20     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 20')
    call addfld ('FDL_int20     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 20')
    call addfld ('FULC_int20    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 20')
    call addfld ('FDLC_int20    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 20')
    call addfld ('FUL_toa_int20     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 20, toa')
    call addfld ('FDL_toa_int20     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 20, toa')
    call addfld ('FULC_toa_int20    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 20, toa')
    call addfld ('FDLC_toa_int20    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 20, toa')

    call addfld ('FUL_int21     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 21')
    call addfld ('FDL_int21     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 21')
    call addfld ('FULC_int21    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 21')
    call addfld ('FDLC_int21    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 21')
    call addfld ('FUL_toa_int21     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 21, toa')
    call addfld ('FDL_toa_int21     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 21, toa')
    call addfld ('FULC_toa_int21    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 21, toa')
    call addfld ('FDLC_toa_int21    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 21, toa')

    call addfld ('FUL_int22     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 22')
    call addfld ('FDL_int22     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 22')
    call addfld ('FULC_int22    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 22')
    call addfld ('FDLC_int22    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 22')
    call addfld ('FUL_toa_int22     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 22, toa')
    call addfld ('FDL_toa_int22     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 22, toa')
    call addfld ('FULC_toa_int22    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 22, toa')
    call addfld ('FDLC_toa_int22    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 22 toa')

    call addfld ('FUL_int23     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 23')
    call addfld ('FDL_int23     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 23')
    call addfld ('FULC_int23    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 23')
    call addfld ('FDLC_int23    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 23')
    call addfld ('FUL_toa_int23     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 23, toa')
    call addfld ('FDL_toa_int23     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 23, toa')
    call addfld ('FULC_toa_int23    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 23, toa')
    call addfld ('FDLC_toa_int23    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 23, toa')

    call addfld ('FUL_int24     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 24')
    call addfld ('FDL_int24     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 24')
    call addfld ('FULC_int24    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 24')
    call addfld ('FDLC_int24    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 24')
    call addfld ('FUL_toa_int24     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 24, toa')
    call addfld ('FDL_toa_int24     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 24, toa')
    call addfld ('FULC_toa_int24    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 24, toa')
    call addfld ('FDLC_toa_int24    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 24, toa')

    call addfld ('FUL_int25     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 25')
    call addfld ('FDL_int25     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 25')
    call addfld ('FULC_int25    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 25')
    call addfld ('FDLC_int25    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 25')
    call addfld ('FUL_toa_int25     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 25, toa')
    call addfld ('FDL_toa_int25     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 25, toa')
    call addfld ('FULC_toa_int25    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 25, toa')
    call addfld ('FDLC_toa_int25    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 25, toa')

    call addfld ('FUL_int26     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 26')
    call addfld ('FDL_int26     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 26')
    call addfld ('FULC_int26    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 26')
    call addfld ('FDLC_int26    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 26')
    call addfld ('FUL_toa_int26     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 26, toa')
    call addfld ('FDL_toa_int26     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 26, toa')
    call addfld ('FULC_toa_int26    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 26, toa')
    call addfld ('FDLC_toa_int26    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 26, toa')

    call addfld ('FUL_int27     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 27')
    call addfld ('FDL_int27     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 27')
    call addfld ('FULC_int27    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 27')
    call addfld ('FDLC_int27    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 27')
    call addfld ('FUL_toa_int27     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 27, toa')
    call addfld ('FDL_toa_int27     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 27, toa')
    call addfld ('FULC_toa_int27    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 27, toa')
    call addfld ('FDLC_toa_int27    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 27, toa')

    call addfld ('FUL_int28     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 28')
    call addfld ('FDL_int28     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 28')
    call addfld ('FULC_int28    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 28')
    call addfld ('FDLC_int28    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 28')
    call addfld ('FUL_toa_int28     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 28, toa')
    call addfld ('FDL_toa_int28     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 28, toa')
    call addfld ('FULC_toa_int28    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 28, toa')
    call addfld ('FDLC_toa_int28    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 28, toa')

    call addfld ('FUL_int29     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 29')
    call addfld ('FDL_int29     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 29')
    call addfld ('FULC_int29    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 29')
    call addfld ('FDLC_int29    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 29')
    call addfld ('FUL_toa_int29     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 29, toa')
    call addfld ('FDL_toa_int29     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 29, toa')
    call addfld ('FULC_toa_int29    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 29, toa')
    call addfld ('FDLC_toa_int29    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 29, toa')

    call addfld ('FUL_int30     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 30')
    call addfld ('FDL_int30     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 30')
    call addfld ('FULC_int30    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 30')
    call addfld ('FDLC_int30    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 30')
    call addfld ('FUL_toa_int30     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 30, toa')
    call addfld ('FDL_toa_int30     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 30, toa')
    call addfld ('FULC_toa_int30    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 30, toa')
    call addfld ('FDLC_toa_int30    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 30, toa')

    call addfld ('FUL_int31     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 31')
    call addfld ('FDL_int31     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 31')
    call addfld ('FULC_int31    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 31')
    call addfld ('FDLC_int31    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 31')
    call addfld ('FUL_toa_int31     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 31, toa')
    call addfld ('FDL_toa_int31     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 31, toa')
    call addfld ('FULC_toa_int31    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 31, toa')
    call addfld ('FDLC_toa_int31    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 31, toa')

    call addfld ('FUL_int32     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 32')
    call addfld ('FDL_int32     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 32')
    call addfld ('FULC_int32    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 32')
    call addfld ('FDLC_int32    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 32')
    call addfld ('FUL_toa_int32     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 32, toa')
    call addfld ('FDL_toa_int32     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 32, toa')
    call addfld ('FULC_toa_int32    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 32, toa')
    call addfld ('FDLC_toa_int32    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 32, toa')

    call addfld ('FUL_int33     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 33')
    call addfld ('FDL_int33     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 33')
    call addfld ('FULC_int33    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 33')
    call addfld ('FDLC_int33    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 33')
    call addfld ('FUL_toa_int33     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 33, toa')
    call addfld ('FDL_toa_int33     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 33, toa')
    call addfld ('FULC_toa_int33    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 33, toa')
    call addfld ('FDLC_toa_int33    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 33, toa')

    call addfld ('FUL_int34     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 34')
    call addfld ('FDL_int34     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 34')
    call addfld ('FULC_int34    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 34')
    call addfld ('FDLC_int34    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 34')
    call addfld ('FUL_toa_int34     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 34, toa')
    call addfld ('FDL_toa_int34     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 34, toa')
    call addfld ('FULC_toa_int34    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 34, toa')
    call addfld ('FDLC_toa_int34    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 34, toa')

    call addfld ('FUL_int35     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 35')
    call addfld ('FDL_int35     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 35')
    call addfld ('FULC_int35    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 35')
    call addfld ('FDLC_int35    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 35')
    call addfld ('FUL_toa_int35     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 35, toa')
    call addfld ('FDL_toa_int35     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 35, toa')
    call addfld ('FULC_toa_int35    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 35, toa')
    call addfld ('FDLC_toa_int35    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 35, toa')

    call addfld ('FUL_int36     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 36')
    call addfld ('FDL_int36     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 36')
    call addfld ('FULC_int36    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 36')
    call addfld ('FDLC_int36    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 36')
    call addfld ('FUL_toa_int36     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 36, toa')
    call addfld ('FDL_toa_int36     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 36, toa')
    call addfld ('FULC_toa_int36    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 36, toa')
    call addfld ('FDLC_toa_int36    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 36, toa')

    call addfld ('FUL_int37     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 37')
    call addfld ('FDL_int37     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 37')
    call addfld ('FULC_int37    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 37')
    call addfld ('FDLC_int37    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 37')
    call addfld ('FUL_toa_int37     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 37, toa')
    call addfld ('FDL_toa_int37     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 37, toa')
    call addfld ('FULC_toa_int37    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 37, toa')
    call addfld ('FDLC_toa_int37    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 37, toa')

    call addfld ('FUL_int38     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 38')
    call addfld ('FDL_int38     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 38')
    call addfld ('FULC_int38    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 38')
    call addfld ('FDLC_int38    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 38')
    call addfld ('FUL_toa_int38     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 38, toa')
    call addfld ('FDL_toa_int38     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 38, toa')
    call addfld ('FULC_toa_int38    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 38, toa')
    call addfld ('FDLC_toa_int38    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 38, toa')

    call addfld ('FUL_int39     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 39')
    call addfld ('FDL_int39     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 39')
    call addfld ('FULC_int39    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 39')
    call addfld ('FDLC_int39    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 39')
    call addfld ('FUL_toa_int39     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 39, toa')
    call addfld ('FDL_toa_int39     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 39, toa')
    call addfld ('FULC_toa_int39    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 39, toa')
    call addfld ('FDLC_toa_int39    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 39, toa')

    call addfld ('FUL_int40     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 40')
    call addfld ('FDL_int40     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 40')
    call addfld ('FULC_int40    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 40')
    call addfld ('FDLC_int40    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 40')
    call addfld ('FUL_toa_int40     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 40, toa')
    call addfld ('FDL_toa_int40     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 40, toa')
    call addfld ('FULC_toa_int40    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 40, toa')
    call addfld ('FDLC_toa_int40    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 40, toa')

    call addfld ('FUL_int41     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 41')
    call addfld ('FDL_int41     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 41')
    call addfld ('FULC_int41    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 41')
    call addfld ('FDLC_int41    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 41')
    call addfld ('FUL_toa_int41     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 41, toa')
    call addfld ('FDL_toa_int41     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 41, toa')
    call addfld ('FULC_toa_int41    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 41, toa')
    call addfld ('FDLC_toa_int41    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 41, toa')

    call addfld ('FUL_int42     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 42')
    call addfld ('FDL_int42     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 42')
    call addfld ('FULC_int42    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 42')
    call addfld ('FDLC_int42    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 42')
    call addfld ('FUL_toa_int42     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 42, toa')
    call addfld ('FDL_toa_int42     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 42, toa')
    call addfld ('FULC_toa_int42    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 42, toa')
    call addfld ('FDLC_toa_int42    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 42, toa')

    call addfld ('FUL_int43     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 43')
    call addfld ('FDL_int43     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 43')
    call addfld ('FULC_int43    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 43')
    call addfld ('FDLC_int43    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 43')
    call addfld ('FUL_toa_int43     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 43, toa')
    call addfld ('FDL_toa_int43     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 43, toa')
    call addfld ('FULC_toa_int43    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 43, toa')
    call addfld ('FDLC_toa_int43    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 43, toa')

    call addfld ('FUL_int44     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 44')
    call addfld ('FDL_int44     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 44')
    call addfld ('FULC_int44    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 44')
    call addfld ('FDLC_int44    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 44')
    call addfld ('FUL_toa_int44     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 44, toa')
    call addfld ('FDL_toa_int44     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 44, toa')
    call addfld ('FULC_toa_int44    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 44, toa')
    call addfld ('FDLC_toa_int44    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 44, toa')

    call addfld ('FUL_int45     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 45')
    call addfld ('FDL_int45     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 45')
    call addfld ('FULC_int45    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 45')
    call addfld ('FDLC_int45    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 45')
    call addfld ('FUL_toa_int45     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 45, toa')
    call addfld ('FDL_toa_int45     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 45, toa')
    call addfld ('FULC_toa_int45    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 45, toa')
    call addfld ('FDLC_toa_int45    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 45, toa')

    call addfld ('FUL_int46     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 46')
    call addfld ('FDL_int46     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 46')
    call addfld ('FULC_int46    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 46')
    call addfld ('FDLC_int46    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 46')
    call addfld ('FUL_toa_int46     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 46, toa')
    call addfld ('FDL_toa_int46     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 46, toa')
    call addfld ('FULC_toa_int46    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 46, toa')
    call addfld ('FDLC_toa_int46    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 46, toa')

    call addfld ('FUL_int47     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 47')
    call addfld ('FDL_int47     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 47')
    call addfld ('FULC_int47    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 47')
    call addfld ('FDLC_int47    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 47')
    call addfld ('FUL_toa_int47     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 48, toa')
    call addfld ('FDL_toa_int47     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 48, toa')
    call addfld ('FULC_toa_int47    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 48, toa')
    call addfld ('FDLC_toa_int47    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 48, toa')

    call addfld ('FUL_int48     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 48')
    call addfld ('FDL_int48     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 48')
    call addfld ('FULC_int48    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 48')
    call addfld ('FDLC_int48    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 48')
    call addfld ('FUL_toa_int48     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 48, toa')
    call addfld ('FDL_toa_int48     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 48, toa')
    call addfld ('FULC_toa_int48    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 48, toa')
    call addfld ('FDLC_toa_int48    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 48, toa')

    call addfld ('FUL_int49     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 49')
    call addfld ('FDL_int49     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 49')
    call addfld ('FULC_int49    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 49')
    call addfld ('FDLC_int49    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 49')
    call addfld ('FUL_toa_int49     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 49, toa')
    call addfld ('FDL_toa_int49     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 49, toa')
    call addfld ('FULC_toa_int49    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 49, toa')
    call addfld ('FDLC_toa_int49    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 49, toa')

    call addfld ('FUL_int50     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 50')
    call addfld ('FDL_int50     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 50')
    call addfld ('FULC_int50    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 50')
    call addfld ('FDLC_int50    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 50')
    call addfld ('FUL_toa_int50     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 50, toa')
    call addfld ('FDL_toa_int50     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 50, toa')
    call addfld ('FULC_toa_int50    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 50, toa')
    call addfld ('FDLC_toa_int50    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 50, toa')

    call addfld ('FUL_int51     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 51')
    call addfld ('FDL_int51     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 51')
    call addfld ('FULC_int51    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 51')
    call addfld ('FDLC_int51    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 51')
    call addfld ('FUL_toa_int51     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 51, toa')
    call addfld ('FDL_toa_int51     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 51, toa')
    call addfld ('FULC_toa_int51    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 51, toa')
    call addfld ('FDLC_toa_int51    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 51, toa')

    call addfld ('FUL_int52     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 52')
    call addfld ('FDL_int52     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 52')
    call addfld ('FULC_int52    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 52')
    call addfld ('FDLC_int52    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 52')
    call addfld ('FUL_toa_int52     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 52, toa')
    call addfld ('FDL_toa_int52     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 52, toa')
    call addfld ('FULC_toa_int52    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 52, toa')
    call addfld ('FDLC_toa_int52    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 52, toa')

    call addfld ('FUL_int53     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 53')
    call addfld ('FDL_int53     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 53')
    call addfld ('FULC_int53    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 53')
    call addfld ('FDLC_int53    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 53')
    call addfld ('FUL_toa_int53     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 53, toa')
    call addfld ('FDL_toa_int53     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 53, toa')
    call addfld ('FULC_toa_int53    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 53, toa')
    call addfld ('FDLC_toa_int53    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 53, toa')

    call addfld ('FUL_int54     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 54')
    call addfld ('FDL_int54     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 54')
    call addfld ('FULC_int54    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 54')
    call addfld ('FDLC_int54    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 54')
    call addfld ('FUL_toa_int54     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 54, toa')
    call addfld ('FDL_toa_int54     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 54, toa')
    call addfld ('FULC_toa_int54    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 54, toa')
    call addfld ('FDLC_toa_int54    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 54, toa')

    call addfld ('FUL_int55     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 55')
    call addfld ('FDL_int55     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 55')
    call addfld ('FULC_int55    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 55')
    call addfld ('FDLC_int55    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 55')
    call addfld ('FUL_toa_int55     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 55, toa')
    call addfld ('FDL_toa_int55     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 55, toa')
    call addfld ('FULC_toa_int55    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 55, toa')
    call addfld ('FDLC_toa_int55    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 55, toa')

    call addfld ('FUL_int56     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 56')
    call addfld ('FDL_int56     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 56')
    call addfld ('FULC_int56    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 56')
    call addfld ('FDLC_int56    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 56')
    call addfld ('FUL_toa_int56     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 56, toa')
    call addfld ('FDL_toa_int56     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 56, toa')
    call addfld ('FULC_toa_int56    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 56, toa')
    call addfld ('FDLC_toa_int56    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 56, toa')

    call addfld ('FUL_int57     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 57')
    call addfld ('FDL_int57     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 57')
    call addfld ('FULC_int57    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 57')
    call addfld ('FDLC_int57    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 57')
    call addfld ('FUL_toa_int57     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 57, toa')
    call addfld ('FDL_toa_int57     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 57, toa')
    call addfld ('FULC_toa_int57    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 57, toa')
    call addfld ('FDLC_toa_int57    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 57, toa')

    call addfld ('FUL_int58     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 58')
    call addfld ('FDL_int58     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 58')
    call addfld ('FULC_int58    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 58')
    call addfld ('FDLC_int58    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 58')
    call addfld ('FUL_toa_int58     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 58, toa')
    call addfld ('FDL_toa_int58     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 58, toa')
    call addfld ('FULC_toa_int58    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 58, toa')
    call addfld ('FDLC_toa_int58    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 58, toa')

    call addfld ('FUL_int59     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 59')
    call addfld ('FDL_int59     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 59')
    call addfld ('FULC_int59    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 59')
    call addfld ('FDLC_int59    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 59')
    call addfld ('FUL_toa_int59     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 59, toa')
    call addfld ('FDL_toa_int59     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 59, toa')
    call addfld ('FULC_toa_int59    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 59, toa')
    call addfld ('FDLC_toa_int59    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 59, toa')

    call addfld ('FUL_int60     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 60')
    call addfld ('FDL_int60     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 60')
    call addfld ('FULC_int60    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 60')
    call addfld ('FDLC_int60    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 60')
    call addfld ('FUL_toa_int60     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 60, toa')
    call addfld ('FDL_toa_int60     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 60, toa')
    call addfld ('FULC_toa_int60    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 60, toa')
    call addfld ('FDLC_toa_int60    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 60, toa')

    call addfld ('FUL_int61     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 61')
    call addfld ('FDL_int61     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 61')
    call addfld ('FULC_int61    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 61')
    call addfld ('FDLC_int61    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 61')
    call addfld ('FUL_toa_int61     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 61, toa')
    call addfld ('FDL_toa_int61     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 61, toa')
    call addfld ('FULC_toa_int61    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 61, toa')
    call addfld ('FDLC_toa_int61    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 61, toa')

    call addfld ('FUL_int62     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 62')
    call addfld ('FDL_int62     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 62')
    call addfld ('FULC_int62    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 62')
    call addfld ('FDLC_int62    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 62')
    call addfld ('FUL_toa_int62     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 62, toa')
    call addfld ('FDL_toa_int62     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 62, toa')
    call addfld ('FULC_toa_int62    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 62, toa')
    call addfld ('FDLC_toa_int62    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 62, toa')

    call addfld ('FUL_int63     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 63')
    call addfld ('FDL_int63     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 63')
    call addfld ('FULC_int63    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 63')
    call addfld ('FDLC_int63    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 63')
    call addfld ('FUL_toa_int63     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 63, toa')
    call addfld ('FDL_toa_int63     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 63, toa')
    call addfld ('FULC_toa_int63    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 63, toa')
    call addfld ('FDLC_toa_int63    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 63, toa')

    call addfld ('FUL_int64     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 64')
    call addfld ('FDL_int64     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 64')
    call addfld ('FULC_int64    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 64')
    call addfld ('FDLC_int64    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 64')
    call addfld ('FUL_toa_int64     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 64, toa')
    call addfld ('FDL_toa_int64     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 64, toa')
    call addfld ('FULC_toa_int64    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 64, toa')
    call addfld ('FDLC_toa_int64    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 64, toa')

    call addfld ('FUL_int65     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 65')
    call addfld ('FDL_int65     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 65')
    call addfld ('FULC_int65    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 65')
    call addfld ('FDLC_int65    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 65')
    call addfld ('FUL_toa_int65     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 65, toa')
    call addfld ('FDL_toa_int65     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 65, toa')
    call addfld ('FULC_toa_int65    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 65, toa')
    call addfld ('FDLC_toa_int65    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 65, toa')

    call addfld ('FUL_int66     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 66')
    call addfld ('FDL_int66     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 66')
    call addfld ('FULC_int66    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 66')
    call addfld ('FDLC_int66    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 66')
    call addfld ('FUL_toa_int66     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 66, toa')
    call addfld ('FDL_toa_int66     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 66, toa')
    call addfld ('FULC_toa_int66    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 66, toa')
    call addfld ('FDLC_toa_int66    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 66, toa')

    call addfld ('FUL_int67     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 67')
    call addfld ('FDL_int67     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 67')
    call addfld ('FULC_int67    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 67')
    call addfld ('FDLC_int67    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 67')
    call addfld ('FUL_toa_int67     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 67, toa')
    call addfld ('FDL_toa_int67     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 67, toa')
    call addfld ('FULC_toa_int67    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 67, toa')
    call addfld ('FDLC_toa_int67    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 67, toa')

    call addfld ('FUL_int68     ',(/'ilev'/),'A','W/m2    ','Longwave upward flux, spectral interval 68')
    call addfld ('FDL_int68     ',(/'ilev'/),'A','W/m2    ','Longwave downward flux, spectral interval 68')
    call addfld ('FULC_int68    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 68')
    call addfld ('FDLC_int68    ',(/'ilev'/),'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 68')
    call addfld ('FUL_toa_int68     ',horiz_only,'A','W/m2    ','Longwave upward flux, spectral interval 68, toa')
    call addfld ('FDL_toa_int68     ',horiz_only,'A','W/m2    ','Longwave downward flux, spectral interval 68, toa')
    call addfld ('FULC_toa_int68    ',horiz_only,'A','W/m2    ','Longwave clear-sky upward flux, spectral interval 68, toa')
    call addfld ('FDLC_toa_int68    ',horiz_only,'A','W/m2    ','Longwave clear-sky downward flux, spectral interval 68, toa')

  end subroutine addfld_spectral_intervals

  end subroutine mars_init
  !==============================================================================

  !==============================================================================
  subroutine mars_condensate_tend(state, ptend, ztodt, pbuf)
    !
    ! mars_condensate_tend: Run the selected process for condensation and sublimation of CO2
    !=====================================================================
    use physics_types,only: physics_state, physics_ptend
    use physics_types,only: physics_ptend_init
    use mars,         only: mars_condensate_NONE,mars_condensate_USER
    !
    ! Passed Variables
    !------------------
    type(physics_state)      ,intent(inout):: state
    real(r8)                 ,intent(in)   :: ztodt
    type(physics_ptend)      ,intent(out)  :: ptend
    type(physics_buffer_desc),pointer      :: pbuf(:)
    !
    ! Local Values
    !-----------------
    real(r8)        :: dtcond(state%ncol, pver) ! Temperature tendency due to condensation
    real(r8)        :: dqcond(state%ncol, pver) ! Q tendency due to condensation
    real(r8)        :: T     (state%ncol, pver) ! T temporary
    real(r8)        :: qco2    (state%ncol, pver) ! Q temporary
    logical         :: lq(pcnst)                ! Calc tendencies?
    integer         :: lchnk                    ! chunk identifier
    integer         :: ncol                     ! number of atmospheric columns
    integer         :: k

    ! Set local copies of values
    !---------------------------------
    lchnk       = state%lchnk
    ncol        = state%ncol
    T (:ncol,:) = state%T(:ncol,:)
    qco2(:ncol,:) = state%Q(:ncol,:,1)

    ! initialize individual parameterization tendencies
    !---------------------------------------------------
    lq    = .false.
    lq(1) = .true.
    call physics_ptend_init(ptend, state%psetcols, 'Mars condensate', &
                                ls=.true., lu=.true., lv=.true., lq=lq)

    ! Get values from the physics buffer
    !------------------------------------

    ! Initialize values for condensate tendencies
    !---------------------------------------------
    do k = 1, pver
      dtcond(:ncol,k) = state%T(:ncol,k)
      dqcond(:ncol,k) = state%q(:ncol,k,1)
    end do

    ! Call the Selected condensation routine  ~~DEVO style~~
    !--------------------------------------------------------
    if(CONDENSATE_OPT == CONDENSATE_NONE) then
      call mars_condensate_NONE(ncol,pver,state%pmid(:ncol,:), &
                                                             T(:ncol,:), &
                                                            qco2(:ncol,:) )
    elseif(CONDENSATE_OPT == CONDENSATE_USER) then
      call mars_condensate_USER(ncol,pver,ztodt,state%pmid(:ncol,:), &
                                                    state%pdel(:ncol,:), &
                                                             T(:ncol,:), &
                                                          qco2(:ncol,:) )
    else
      ! ERROR: Unknown CONDENSATE_OPT value
      !-------------------------------------
      write(iulog,*) 'ERROR: unknown CONDENSATE_OPT=',CONDENSATE_OPT
      call endrun('mars_condensate_tend() CONDENSATE_OPT ERROR')
    endif

    ! Back out temperature and specific humidity
    ! tendencies from updated fields
    !--------------------------------------------
    do k = 1, pver
      ptend%s(:ncol,k)   = (T (:,k)-state%T(:ncol,k)  )/ztodt*cpair
      ptend%q(:ncol,k,1) = (qco2(:,k)-state%q(:ncol,k,1))/ztodt
    end do

    ! Output condensate tendencies
    !------------------------------
    do k = 1, pver
      dtcond(:ncol,k) = (T (:ncol,k) - dtcond(:ncol,k))/ztodt
      dqcond(:ncol,k) = (qco2(:ncol,k) - dqcond(:ncol,k))/ztodt
    end do

  end subroutine mars_condensate_tend
  !==============================================================================

  subroutine mars_gw_drag_tend(state, ptend, ztodt, cam_in)
    !
    ! mars_gw_drag_tend: Run the selected GW_DRAG process.
    !=========================================================================
    use physics_types,only: physics_state, physics_ptend
    use physics_types,only: physics_ptend_init
    use phys_grid,    only: get_rlat_all_p
    use mars,         only: mars_gw_drag_NONE,mars_gw_drag_USER
    !
    ! Passed Variables
    !-------------------
    type(physics_state),intent(in)   :: state
    real(r8),           intent(in)   :: ztodt
    type(physics_ptend),intent(out)  :: ptend
    type(cam_in_t),     intent(inout):: cam_in
    !
    ! Local Values
    !----------------
    real(r8) :: T         (state%ncol,pver)   ! T temporary
    real(r8) :: qv        (state%ncol,pver)   ! Q temporary (specific humidity)
    real(r8) :: U         (state%ncol,pver)   ! U temporary
    real(r8) :: V         (state%ncol,pver)   ! V temporary
    real(r8) :: dqdt_vdiff(state%ncol,pver)   ! GW_DRAG Q vertical diffusion tend kg/kg/s
    real(r8) :: dtdt_vdiff(state%ncol,pver)   ! GW_DRAG T vertical diffusion tend  K/s
    real(r8) :: dudt_vdiff(state%ncol,pver)   ! GW_DRAG U vertical diffusion tend  m/s/s
    real(r8) :: dvdt_vdiff(state%ncol,pver)   ! GW_DRAG V vertical diffusion tend  m/s/s
    real(r8) :: Km        (state%ncol,pverp)  ! Eddy diffusivity at layer interfaces (m2/s)
    real(r8) :: Ke        (state%ncol,pverp)  ! Eddy diffusivity at layer interfaces (m2/s)
    real(r8) :: VSE       (state%ncol,pver)   ! Dry Static Energy divided by Cp (K)
    real(r8) :: Zm        (state%ncol,pver)   !
    real(r8) :: Zi        (state%ncol,pver)   !
    real(r8) :: Z_gw_drag     (state%ncol)        !
    real(r8) :: Rf        (state%ncol,pver)   !
    real(r8) :: Tsfc      (state%ncol)        ! Surface T
    real(r8) :: Qsfc      (state%ncol)        ! Surface Q (saturated)
    real(r8) :: Cdrag     (state%ncol)        ! Cdrag coef from surface calculation

    logical  :: lq        (pcnst)             ! Calc tendencies?
    real(r8) :: dTs       (state%ncol)
    real(r8) :: dUa       (state%ncol,pver)
    real(r8) :: dVa       (state%ncol,pver)
    real(r8) :: dTa       (state%ncol,pver)
    real(r8) :: dQa       (state%ncol,pver)
    integer  :: lchnk                        ! chunk identifier
    integer  :: ncol                         ! number of atmospheric columns
    integer  :: kk                           ! loop index

    ! Set local copies of values
    !---------------------------------
    lchnk              = state%lchnk
    ncol               = state%ncol
    Zm  (:ncol,:)      = state%zm  (:ncol,:)
    Zi  (:ncol,1:pver) = state%zi  (:ncol,1:pver)
    T   (:ncol,:)      = state%T   (:ncol,:)
    U   (:ncol,:)      = state%U   (:ncol,:)
    V   (:ncol,:)      = state%V   (:ncol,:)
    qv  (:ncol,:)      = state%Q   (:ncol,:,1)

    ! Initialize individual parameterization tendencies
    !-----------------------------------------------------
    lq    = .false.
    lq(1) = .true.
    call physics_ptend_init(ptend,state%psetcols,'Mars gw_drag_tend',        &
                                       ls=.true., lu=.true., lv=.true., lq=lq)

    ! Call the Selected GW_DRAG routine
    !--------------------------------------------------------
    Tsfc(:ncol) = Tsurf(:ncol,lchnk)
    Qsfc(:ncol) = Qsurf(:ncol,lchnk)
    if(GW_DRAG_OPT == GW_DRAG_NONE) then
      ! Call Mars GW_DRAG scheme
      !--------------------------------------------------
      call mars_gw_drag_NONE()
    elseif(GW_DRAG_OPT == GW_DRAG_USER) then
      ! Call USER implemented routine in mars module
      !--------------------------------------------------
      call mars_gw_drag_USER()
    else
      ! ERROR: Unknown GW_DRAG_OPT value
      !-------------------------------------
      write(iulog,*) 'ERROR: unknown GW_DRAG_OPT=',GW_DRAG_OPT
      call endrun('mars_gw_drag_tend() GW_DRAG_OPT ERROR')
    endif
    Tsurf(:ncol,lchnk) = Tsfc (:ncol)
    Qsurf(:ncol,lchnk) = Qsfc (:ncol)
    Cd   (:ncol,lchnk) = Cdrag(:ncol)

    ! Back out tendencies from updated fields
    !-----------------------------------------
    do kk = 1, pver
      ptend%s(:ncol,kk  ) = (T (:,kk)-state%T(:ncol,kk  ))/ztodt*cpair
      ptend%u(:ncol,kk  ) = (U (:,kk)-state%U(:ncol,kk  ))/ztodt
      ptend%v(:ncol,kk  ) = (V (:,kk)-state%V(:ncol,kk  ))/ztodt
      ptend%q(:ncol,kk,1) = (qv(:,kk)-state%q(:ncol,kk,1))/ztodt
    end do

    ! Archive diagnostic fields
    !----------------------------
    call outfld('Tsurf' ,Tsurf(:ncol,lchnk) ,ncol,lchnk)
    call outfld('Qsurf' ,Qsurf(:ncol,lchnk) ,ncol,lchnk)
    call outfld('Cdrag' ,Cd   (:ncol,lchnk) ,ncol,lchnk)
    call outfld('Zgw_drag'  ,Z_gw_drag              ,ncol,lchnk) !
    call outfld('KVH'   ,Ke                 ,ncol,lchnk) ! Eddy diffusivity (heat and moisture,m2/s)
    call outfld('KVM'   ,Km                 ,ncol,lchnk) ! Eddy diffusivity (momentum, m2/s)
    call outfld('VSE'   ,VSE                ,ncol,lchnk) ! Virtual Dry Static Energy divided by Cp (K)
    call outfld('Zm'    ,Zm                 ,ncol,lchnk) !
    call outfld('Rf'    ,Rf                 ,ncol,lchnk) !
    call outfld('DTV'   ,dtdt_vdiff         ,ncol,lchnk) ! GW_DRAG + surface flux T tendency (K/s)
    call outfld('DUV'   ,dudt_vdiff         ,ncol,lchnk) ! GW_DRAG u tendency (m/s2)
    call outfld('DVV'   ,dvdt_vdiff         ,ncol,lchnk) ! GW_DRAG v tendency (m/s2)
    call outfld('VD01'  ,dqdt_vdiff         ,ncol,lchnk) ! GW_DRAG + surface flux Q tendency (kg/kg/s)
    call outfld('SHflux',SHflux(:ncol,lchnk),ncol,lchnk) ! Sensible Heat Flux
    call outfld('LHflux',LHflux(:ncol,lchnk),ncol,lchnk) ! Latent Heat Flux
    call outfld('TauU'  ,TUflux(:ncol,lchnk),ncol,lchnk) ! U Surface Stress
    call outfld('TauV'  ,TVflux(:ncol,lchnk),ncol,lchnk) ! V Surface Stress

end subroutine mars_gw_drag_tend
!============================================================================

subroutine mars_surface_init(ncol, clat, PS, Tsfc, Qsfc)
  !
  !
  !==========================================================================
  !
  ! Passed variables
  !--------------------
  integer ,intent(in) :: ncol
  real(r8),intent(in) :: clat (ncol)
  real(r8),intent(in) :: PS   (ncol)
  real(r8),intent(out):: Tsfc(ncol)
  real(r8),intent(out):: Qsfc(ncol)
  !
  ! Local values
  !--------------
  integer :: ii
  real(r8):: T_width

end subroutine mars_surface_init
!=======================================================================


!=======================================================================
subroutine mars_restart_init(File,hdimids,hdimcnt)
  !
  ! mars_restart_init:
  !==========================================================================
  !
  ! Passed variables
  !--------------------
  type(file_desc_t),intent(inout):: File
  integer          ,intent(in)   :: hdimcnt
  integer          ,intent(in)   :: hdimids(1:hdimcnt)
  !
  ! Local values
  !--------------
  integer:: ierr

  ierr = pio_def_var(File,'Mars_Tsfc',pio_double, hdimids, Tsurf_desc)
  if (ierr /= 0) then
     call endrun('mars_restart_init: ERROR defining Mars_Tsfc')
  end if

  ierr = pio_def_var(File,'Mars_Qsfc',pio_double, hdimids, Qsurf_desc)
  if (ierr /= 0) then
     call endrun('mars_restart_init: ERROR defining Mars_Qsfc')
  end if

end subroutine mars_restart_init
!=======================================================================


!=======================================================================
subroutine mars_restart_write(File)
  !
  ! mars_restart_write:
  !==========================================================================
  !
  ! Passed variables
  !--------------------
  type(file_desc_t),intent(inout):: File
  !
  ! Local values
  !--------------
  type(io_desc_t),pointer:: iodesc
  integer:: dims(3),gdims(3),nhdims
  integer:: physgrid
  integer:: ierr

  ! Get the iodesc for write calls
  !---------------------------------
  dims(1) = pcols
  dims(2) = endchunk - begchunk + 1
  physgrid = cam_grid_id('physgrid')
  call cam_grid_dimensions(physgrid, gdims(1:2), nhdims)
  call cam_grid_get_decomp(physgrid,  dims(1:2), gdims(1:nhdims), pio_double, iodesc)

  ! Write Surface values
  !---------------------
  call pio_write_darray(File, Tsurf_desc, iodesc, Tsurf, ierr)
  if (ierr /= 0) then
     call endrun('mars_restart_write: ERROR writing Tsurf')
  end if

  call pio_write_darray(File, Qsurf_desc, iodesc, Qsurf, ierr)
  if (ierr /= 0) then
     call endrun('mars_restart_write: ERROR writing Qsurf')
  end if

end subroutine mars_restart_write
!=======================================================================


!=======================================================================
subroutine mars_restart_read(File)
  !
  ! mars_restart_read:
  !==========================================================================
  use error_messages,only: alloc_err
  !
  ! Passed variables
  !--------------------
  type(file_desc_t),intent(inout):: File
  !
  ! Local values
  !--------------
  type( io_desc_t),pointer:: iodesc
  type(var_desc_t)        :: vardesc
  integer:: dims(3),gdims(3),nhdims
  integer:: physgrid
  integer:: ierr

  ! Allocate space for the restart fields
  !-----------------------------------------
  allocate(Tsurf (pcols,begchunk:endchunk),stat=ierr)
  call alloc_err(ierr,'Mars RESTART','Tsurf' ,pcols*(endchunk-begchunk+1))
  allocate(Qsurf (pcols,begchunk:endchunk)  ,stat=ierr)
  call alloc_err(ierr,'Mars RESTART','Qsurf' ,pcols*(endchunk-begchunk+1))

  ! Get the iodesc for read calls
  !---------------------------------
  dims(1) = pcols
  dims(2) = endchunk - begchunk + 1
  physgrid = cam_grid_id('physgrid')
  call cam_grid_dimensions(physgrid, gdims(1:2), nhdims)
  call cam_grid_get_decomp(physgrid,  dims(1:2), gdims(1:nhdims), pio_double, iodesc)

  ! Read Surface values
  !---------------------
  ierr = pio_inq_varid(File,'Mars_Tsfc',vardesc)
  if (ierr /= 0) then
     call endrun('mars_restart_read: ERROR PIO unable to find variable Mars_Tsfc')
  end if

  call pio_read_darray(File, vardesc, iodesc, Tsurf, ierr)
  if (ierr /= 0) then
     call endrun('mars_restart_read: ERROR PIO unable to read variable Tsurf')
  end if

  ierr = pio_inq_varid(File,'Mars_Qsfc',vardesc)
  if (ierr /= 0) then
     call endrun('mars_restart_read: ERROR PIO unable to find variable Mars_Qsfc')
  end if

  call pio_read_darray(File, vardesc, iodesc, Qsurf, ierr)
  if (ierr /= 0) then
     call endrun('mars_restart_read: ERROR PIO unable to read variable Qsurf')
  end if

end subroutine mars_restart_read
!============================================================================

  function mars_radiation_do(timestep)

!------------------------------------------------------------------------
!
! Purpose:  Returns true if the exo_rt is done this timestep
!
!------------------------------------------------------------------------

    integer, intent(in), optional :: timestep
    logical :: mars_radiation_do

!------------------------------------------------------------------------
!
! Local Variables
!
    integer :: nstep
!------------------------------------------------------------------------
!
! Start Code
!
  if (present(timestep)) then
      nstep = timestep
   else
      nstep = get_nstep()
   end if

!   write(*,*)  "mars_radiation_do ", nstep, exo_rad_step, mod(nstep-1,exo_rad_step)
!    mars_radiation_do = nstep == 0  .or.  exo_rad_step == 1                     &
!                       .or. (mod(nstep-1,exo_rad_step) == 0  .and.  nstep /= 1)

    mars_radiation_do = nstep == 0  .or.  exo_rad_step == 1                     &
                       .or. (mod(nstep,exo_rad_step) == 0  .and.  nstep /= 1)

  end function mars_radiation_do
  real(r8) function mars_radiation_nextsw_cday()

    !-----------------------------------------------------------------------
    ! Purpose: Returns calendar day of next sw radiation calculation
    !          This is used to ensure surface albedo claculations are sync'd
    !          with mars radiation calls.
    !-----------------------------------------------------------------------

    use time_manager, only: get_curr_calday, get_nstep, get_step_size


    ! Local variables
    integer :: nstep      ! timestep counter
    logical :: dosw       ! true => do shosrtwave calc
    integer :: offset     ! offset for calendar day calculation
    integer :: dTime      ! integer timestep size
    real(r8):: calday     ! calendar day of
    !-----------------------------------------------------------------------

    mars_radiation_nextsw_cday = -1._r8
    dosw   = .false.
    nstep  = get_nstep()
    dtime  = get_step_size()
    offset = 0
    do while (.not. dosw)
       nstep = nstep + 1
       offset = offset + dtime
       if (mars_radiation_do(nstep)) then
          mars_radiation_nextsw_cday = get_curr_calday(offset=offset)
          dosw = .true.
       end if
    end do
    if(mars_radiation_nextsw_cday == -1._r8) then
       call endrun('error in mars_radiation_nextsw_cday')
    end if

  end function mars_radiation_nextsw_cday

!============================================================================

end module mars_cam
