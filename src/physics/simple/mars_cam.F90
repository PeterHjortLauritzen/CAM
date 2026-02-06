module mars_cam
!-----------------------------------------------------------------------
!
! Purpose: Implement Mars idealized forcings
!
!============================================================================
  ! Useful modules
  !-------------------
  use cam_abortutils, only: endrun
  use cam_grid_support,only: cam_grid_id, cam_grid_dimensions, cam_grid_get_decomp
  use cam_history,         only: addfld, add_default, outfld, register_vector_field
  use cam_history_support, only: fillvalue
  use cam_logfile,    only: iulog
  use camsrfexch,     only: cam_in_t,cam_out_t
  use constituents,   only: pcnst
  use co2_cycle,      only: co2_transport
  use physconst,      only: cpair
  use physics_buffer, only: dtype_r8, pbuf_add_field, physics_buffer_desc
  use physics_types,  only: physics_ptend_init, physics_ptend_dealloc
  use physics_types,  only: physics_state, physics_ptend
  use pio,            only: file_desc_t, io_desc_t, var_desc_t, &
                            pio_double, pio_int, pio_noerr, &
                            pio_seterrorhandling, pio_bcast_error, &
                            pio_inq_varid, &
                            pio_def_var, pio_def_dim, &
                            pio_put_var, pio_get_var
  use ppgrid,         only: pcols, pver, pverp, begchunk, endchunk
  use radiation,      only: radiation_define_restart, radiation_write_restart, &
                            radiation_read_restart
  use shr_const_mod,   only: SHR_CONST_STEBOL, SHR_CONST_REARTH, SHR_CONST_KARMAN, SHR_CONST_TKTRIP
  use shr_const_mod,  only: pi => shr_const_pi
  use shr_kind_mod,   only: r8 => shr_kind_r8, cl=>shr_kind_cl
  use held_suarez_cam,  only: held_suarez_init, held_suarez_tend
  use spmd_utils,       only: masterproc
  use string_utils,     only: to_upper
  use subcol_utils,     only: is_subcol_on

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
  public :: mars_radiative_tend
  public :: mars_restart_init
  public :: mars_restart_write
  public :: mars_restart_read
  public :: mars_grow_topo
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

! NEW: Topography Ramp Parameter
  real(r8), public, protected :: mars_topo_ramp_days = 0.0_r8 ! Days to ramp up topo (0 = off)

  integer :: &
         ixcldliq = -1,      &! cloud liquid amount index
         ixcldice = -1        ! cloud ice amount index

!
! Private data
!

    type(var_desc_t) :: flwds_desc, &
         solld_desc, co2prog_desc, co2diag_desc, sols_desc, soll_desc, &
         solsd_desc

    type(var_desc_t) :: bcphidry_desc, bcphodry_desc, ocphidry_desc, ocphodry_desc, &
       dstdry1_desc, dstdry2_desc, dstdry3_desc, dstdry4_desc

    type(var_desc_t) :: cflx_desc(pcnst)

    type(var_desc_t) :: wsx_desc
    type(var_desc_t) :: wsy_desc
    type(var_desc_t) :: shf_desc


contains
  !==============================================================================
  subroutine mars_register()
    !
    ! mars_register:
    !=====================================================================
    use physconst,      only: cpair, mwh2o
    use constituents,   only: cnst_add
    use physics_buffer, only: pbuf_add_field, dtype_r8
    use radiation,      only: radiation_register
    use radheat,        only: radheat_register

    integer                :: mm

    call radiation_register()
    call radheat_register

    call cnst_add('CO2', 44._r8, 800._r8, 1.e-12_r8, mm, fixed_ubc=.false., &
         longname='CO2', readiv=.true., is_convtran1=.true.)

    call cnst_add('N2', 28._r8, 800._r8, 1.e-12_r8, mm, fixed_ubc=.false., &
         longname='N2', readiv=.true., is_convtran1=.true.)

    ! Add liquid constituents
    call cnst_add('CLDLIQ', mwh2o, cpair, 0._r8, ixcldliq,                    &
         longname='Grid box averaged cloud liquid amount', is_convtran1=.true.)
    ! Add ice constituents
    call cnst_add('CLDICE', mwh2o, cpair, 0._r8, ixcldice, &
      longname='Grid box averaged cloud ice amount', is_convtran1=.true.)

! Register Gravity Wave Diagnostics (so they appear in history files)
      call addfld('UTEND_GWDTOT', (/'lev'/), 'A', 'm/s2', 'Zonal wind tendency by all GWs')
      call addfld('VTEND_GWDTOT', (/'lev'/), 'A', 'm/s2', 'Meridional wind tendency by all GWs')
      call register_vector_field('UTEND_GWDTOT', 'VTEND_GWDTOT')

  end subroutine mars_register
  !==============================================================================


  !==============================================================================
  subroutine mars_readnl(nlfile)
    !
    ! mars_readnl: Read in parameters controlling Mars parameterizations.
    !=====================================================================
    use namelist_utils,only: find_group_name
    use units         ,only: getunit, freeunit
    use spmd_utils,    only: masterproc, masterprocid, mpi_logical, mpicom, mpi_real8, &
                             mpi_character
    use radiation,     only: radiation_readnl
    use rad_constituents,    only: rad_cnst_readnl
    ! Input Parameters
    !------------------
    character(len=*),intent(in):: nlfile
    !
    ! Local Values
    !--------------
    integer:: ierr,unitn,i
    logical  :: adv

    character(len=*), parameter :: sub = 'mars_readnl'

    namelist /mars_nl/ mars_topo_ramp_days

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

    call mpi_bcast(mars_topo_ramp_days, 1, mpi_real8, masterprocid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: mars_topo_ramp_days")

   if (masterproc) then
      write(iulog,*) 'MARS Simple Physics namelist parameters:'
      write(iulog,10) mars_topo_ramp_days
   end if

10 format(' ars_topo_ramp_days    :                         ',f8.1/)

  end subroutine mars_readnl
  !==============================================================================


  !==============================================================================
  subroutine mars_init(phys_state,pbuf2d)
    !
    ! mars_init: allocate space for global arrays and initialize values.
    !                Add variables to history outputs
    !=====================================================================
    use physics_types, only: physics_state
    use radiation,     only: radiation_init
    use radheat,     only: radheat_init
    use rad_constituents,    only: rad_cnst_init
    use ref_pres,           only: pref_mid
    use constituents,  only: apcnst, bpcnst, cnst_name

    !
    ! Input Parameters
    !------------------
    type(physics_state)      ,pointer:: phys_state(:)
    type(physics_buffer_desc),pointer:: pbuf2d    (:,:)
    !
    ! Local Values
    !---------------

    ! Until radiation is working we will run it as a diagnostic calculation and use the
    ! heating rates generated from held suarez

    call held_suarez_init()

    !initialize radiation scheme
    call radiation_init(pbuf2d)
    call radheat_init(pref_mid)
    call rad_cnst_init()

!  These needed for cldliq/ice diagnostics - added for vdiff
    call addfld(apcnst(ixcldliq), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldliq))//' after physics'  )
   call addfld(apcnst(ixcldice), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldice))//' after physics'  )
   call addfld(bpcnst(ixcldliq), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldliq))//' before physics' )
   call addfld(bpcnst(ixcldice), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldice))//' before physics' )

   if (ixcldliq > 0) then
      call addfld (cnst_name(ixcldliq),(/ 'lev' /), 'A', 'kg/kg',' Cloud Liquid '      )
   end if
   if (ixcldice > 0) then
      call addfld (cnst_name(ixcldice),(/ 'lev' /), 'A', 'kg/kg',' Cloud Ice ')
   end if


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
    ! Input Parameters
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
    use mars,         only: mars_gw_drag_NONE,mars_gw_drag_USER
    !
    ! Input Parameters
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

subroutine mars_surface_init(cam_in,state)
  !
  !
  !==========================================================================
  !
  ! Passed variables
  !--------------------
  type(cam_in_t),     intent(inout)          :: cam_in
  type(physics_state),intent(in)             :: state
  !
  ! Local values
  !--------------
  integer ncol


  ncol=cam_in%ncol
  write(6,*)'max cam_in%ts(:)=',maxval(cam_in%ts(:ncol))
  write(6,*)'max cam_in%landfrac(:)=',maxval(cam_in%landfrac(:ncol))
  write(6,*)'max cam_in%asdir(:) should be .17 =',maxval(cam_in%asdir(:ncol))
  write(6,*)'max cam_in%asdif(:) should be .17 =',maxval(cam_in%asdif(:ncol))
  write(6,*)'max cam_in%aldir(:) should be .35 =',maxval(cam_in%aldir(:ncol))
  write(6,*)'max cam_in%aldif(:) should be .35 =',maxval(cam_in%aldif(:ncol))
  cam_in%ts(:)= state%t(:,pver)
  cam_in%landfrac(:)= 1.0_r8
  !asdir, asdif (Visible)	0.15 to 0.20	Mars is dark in the visible (absorbs blue/green).
  !aldir, aldif (Near-IR)	0.30 to 0.40	Mars is bright in the infrared (reflects heat).
  !Broadband Average	~0.25	Matches global observation.
  ! Visible Band (Darker)
  cam_in%asdir(:) = 0.17_r8
  cam_in%asdif(:) = 0.17_r8

  ! Near-IR Band (Brighter)
  cam_in%aldir(:) = 0.35_r8
  cam_in%aldif(:) = 0.35_r8
end subroutine mars_surface_init
!=======================================================================

subroutine mars_radiative_tend( state, ptend, pbuf, cam_out, cam_in, net_flx)

  !-----------------------------------------------------------------------
  !
  ! Driver for radiation computation.
  !
  ! Revision history:
  !-----------------------------------------------------------------------
  ! mars_radiative_tend: Run the radiative process
  !=========================================================================
  use radiation,        only: radiation_tend
  use time_manager,     only: get_step_size
  use physics_types,    only: physics_ptend_init
  !
  ! Input Parameters
  !------------------
  type(physics_state),intent(in)   :: state
  type(physics_ptend),intent(out)  :: ptend
  type(physics_buffer_desc), pointer      :: pbuf(:)
  type(cam_in_t),     intent(inout)       :: cam_in
  type(cam_out_t),    intent(inout)       :: cam_out
  real(r8),           intent(out)         :: net_flx(:)
  !
  ! local
  !------------------
  real(r8)                                :: ztodt
  logical                                 :: lq(pcnst)               ! Calc tendencies?
  integer                                 :: k
  integer                                 :: lchnk                   ! chunk identifier
  integer                                 :: ncol                    ! number of atmospheric columns
  integer i


  lchnk       = state%lchnk
  ncol        = state%ncol

  ! initialize individual parameterization tendencies
  !---------------------------------------------------
  lq(:) = .false.
  call physics_ptend_init(ptend, state%psetcols, 'mars exo-rt radiative_tend',   &
       ls=.true., lu=.false., lv=.false., lq=lq)

  ! ExoRT radiation tendency
  call radiation_tend( state, ptend, pbuf, cam_out, cam_in, net_flx)

end subroutine mars_radiative_tend

!=======================================================================
subroutine mars_restart_init(File,hdimids,hdimcnt)
  !
  ! mars_restart_init:
  !==========================================================================
  use radiation,          only: radiation_define_restart
  use subcol,             only: subcol_init_restart
  !
  ! Passed variables
  !--------------------
  type(file_desc_t),intent(inout):: File
  integer          ,intent(in)   :: hdimcnt
  integer          ,intent(in)   :: hdimids(1:hdimcnt)
  !
  ! Local values
  !--------------
  integer            :: i, ncol, ierr
  character(len=4)   :: num

    ierr = pio_def_var(File, 'FLWDS', pio_double, hdimids, flwds_desc)
    ierr = pio_def_var(File, 'SOLS', pio_double, hdimids, sols_desc)
    ierr = pio_def_var(File, 'SOLL', pio_double, hdimids, soll_desc)
    ierr = pio_def_var(File, 'SOLSD', pio_double, hdimids, solsd_desc)
    ierr = pio_def_var(File, 'SOLLD', pio_double, hdimids, solld_desc)

    ierr = pio_def_var(File, 'BCPHIDRY', pio_double, hdimids, bcphidry_desc)
    ierr = pio_def_var(File, 'BCPHODRY', pio_double, hdimids, bcphodry_desc)
    ierr = pio_def_var(File, 'OCPHIDRY', pio_double, hdimids, ocphidry_desc)
    ierr = pio_def_var(File, 'OCPHODRY', pio_double, hdimids, ocphodry_desc)
    ierr = pio_def_var(File, 'DSTDRY1',  pio_double, hdimids, dstdry1_desc)
    ierr = pio_def_var(File, 'DSTDRY2',  pio_double, hdimids, dstdry2_desc)
    ierr = pio_def_var(File, 'DSTDRY3',  pio_double, hdimids, dstdry3_desc)
    ierr = pio_def_var(File, 'DSTDRY4',  pio_double, hdimids, dstdry4_desc)

    if(co2_transport()) then
      ierr = pio_def_var(File, 'CO2PROG', pio_double, hdimids, co2prog_desc)
      ierr = pio_def_var(File, 'CO2DIAG', pio_double, hdimids, co2diag_desc)
    end if

    ! cam_import variables -- write the constituent surface fluxes as individual 2D arrays
    ! rather than as a single variable with a pcnst dimension.  Note that the cflx components
    ! are only needed for those constituents that are not passed to the coupler.  The restart
    ! for constituents passed through the coupler are handled by the .rs. restart file.  But
    ! we don't currently have a mechanism to know whether the constituent is handled by the
    ! coupler or not, so we write all of cflx to the CAM restart file.
    do i = 1, pcnst
      write(num,'(i4.4)') i
      ierr = pio_def_var(File, 'CFLX'//num,  pio_double, hdimids, cflx_desc(i))
    end do

    ierr = pio_def_var(File, 'wsx',  pio_double, hdimids, wsx_desc)
    ierr = pio_def_var(File, 'wsy',  pio_double, hdimids, wsy_desc)
    ierr = pio_def_var(File, 'shf',  pio_double, hdimids, shf_desc)

    call radiation_define_restart(file)

    if (is_subcol_on()) then
      call subcol_init_restart(file, hdimids)
    end if

!    call carma_restart_init(file)
end subroutine mars_restart_init
!=======================================================================


!=======================================================================
subroutine mars_restart_write (File, cam_in, cam_out)
  !
  ! mars_restart_write:
  !==========================================================================
  !
  use pio,                 only: pio_write_darray
  !
  ! Arguments
  !
  type(file_desc_t),   intent(inout) :: File
  type(cam_in_t),    intent(in)    :: cam_in(begchunk:endchunk)
  type(cam_out_t),   intent(in)    :: cam_out(begchunk:endchunk)
  !
  ! Local values
  !--------------
  type(io_desc_t),pointer :: iodesc
  real(r8):: tmpfield(pcols, begchunk:endchunk)
  integer :: i, m          ! loop index
  integer :: ncol          ! number of vertical columns
  integer :: ierr
  integer :: physgrid
  integer :: dims(3), gdims(3)
  integer :: nhdims

  ! cam_in/out variables
  ! This is a group of surface variables so can reuse dims
  dims(1) = pcols
  dims(2) = endchunk - begchunk + 1

  physgrid = cam_grid_id('physgrid')

  call cam_grid_dimensions(physgrid, gdims(1:2))

  if (gdims(2) == 1) then
     nhdims = 1
  else
     nhdims = 2
  end if
  call cam_grid_get_decomp(physgrid, dims(1:2), gdims(1:nhdims), &
       pio_double, iodesc)

      ! cam_out components
      do i = begchunk, endchunk
        ncol = cam_out(i)%ncol
        tmpfield(:ncol, i) = cam_out(i)%flwds(:ncol)
        ! Only have to do this once (cam_in/out vars all same shape)
        if (ncol < pcols) then
          tmpfield(ncol+1:, i) = fillvalue
        end if
      end do
      call pio_write_darray(File, flwds_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
        ncol = cam_out(i)%ncol
        tmpfield(:ncol, i) = cam_out(i)%sols(:ncol)
      end do
      call pio_write_darray(File, sols_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
        ncol = cam_out(i)%ncol
        tmpfield(:ncol, i) = cam_out(i)%soll(:ncol)
      end do
      call pio_write_darray(File, soll_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
        ncol = cam_out(i)%ncol
        tmpfield(:ncol, i) = cam_out(i)%solsd(:ncol)
      end do
      call pio_write_darray(File, solsd_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
        ncol = cam_out(i)%ncol
        tmpfield(:ncol, i) = cam_out(i)%solld(:ncol)
      end do
      call pio_write_darray(File, solld_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
        ncol = cam_out(i)%ncol
        tmpfield(:ncol, i) = cam_out(i)%bcphidry(:ncol)
      end do
      call pio_write_darray(File, bcphidry_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
        ncol = cam_out(i)%ncol
        tmpfield(:ncol, i) = cam_out(i)%bcphodry(:ncol)
      end do
      call pio_write_darray(File, bcphodry_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
        ncol = cam_out(i)%ncol
        tmpfield(:ncol, i) = cam_out(i)%ocphidry(:ncol)
      end do
      call pio_write_darray(File, ocphidry_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
        ncol = cam_out(i)%ncol
        tmpfield(:ncol, i) = cam_out(i)%ocphodry(:ncol)
      end do
      call pio_write_darray(File, ocphodry_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
        ncol = cam_out(i)%ncol
        tmpfield(:ncol, i) = cam_out(i)%dstdry1(:ncol)
      end do
      call pio_write_darray(File, dstdry1_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
        ncol = cam_out(i)%ncol
        tmpfield(:ncol, i) = cam_out(i)%dstdry2(:ncol)
      end do
      call pio_write_darray(File, dstdry2_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
        ncol = cam_out(i)%ncol
        tmpfield(:ncol, i) = cam_out(i)%dstdry3(:ncol)
      end do
      call pio_write_darray(File, dstdry3_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
        ncol = cam_out(i)%ncol
        tmpfield(:ncol, i) = cam_out(i)%dstdry4(:ncol)
      end do
      call pio_write_darray(File, dstdry4_desc, iodesc, tmpfield, ierr)

      if (co2_transport()) then
        do i = begchunk, endchunk
          ncol = cam_out(i)%ncol
          tmpfield(:ncol, i) = cam_out(i)%co2prog(:ncol)
        end do
        call pio_write_darray(File, co2prog_desc, iodesc, tmpfield, ierr)

        do i = begchunk, endchunk
          ncol = cam_out(i)%ncol
          tmpfield(:ncol, i) = cam_out(i)%co2diag(:ncol)
        end do
        call pio_write_darray(File, co2diag_desc, iodesc, tmpfield, ierr)
      end if

      ! cam_in components
      do m = 1, pcnst
        do i = begchunk, endchunk
          ncol = cam_in(i)%ncol
          tmpfield(:ncol, i) = cam_in(i)%cflx(:ncol, m)
        end do
        call pio_write_darray(File, cflx_desc(m), iodesc, tmpfield, ierr)
      end do

      do i = begchunk, endchunk
        ncol = cam_in(i)%ncol
        tmpfield(:ncol,i) = cam_in(i)%wsx(:ncol)
      end do
      call pio_write_darray(File, wsx_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
        ncol = cam_in(i)%ncol
        tmpfield(:ncol,i) = cam_in(i)%wsy(:ncol)
      end do
      call pio_write_darray(File, wsy_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
        ncol = cam_in(i)%ncol
        tmpfield(:ncol,i) = cam_in(i)%shf(:ncol)
      end do
      call pio_write_darray(File, shf_desc, iodesc, tmpfield, ierr)

      call radiation_write_restart(file)
!      call carma_restart_write(file)

end subroutine mars_restart_write
!=======================================================================


!=======================================================================
subroutine mars_restart_read(File, cam_in, cam_out)
  !
  ! mars_restart_read:
  !==========================================================================
  use subcol,              only: subcol_read_restart
  use pio,                 only: pio_read_darray
  !
  ! Arguments
  !
  type(file_desc_t), intent(inout)   :: File
  type(cam_in_t),            pointer :: cam_in(:)
  type(cam_out_t),           pointer :: cam_out(:)

  !
  ! Local workspace
  !
  real(r8), allocatable :: tmpfield2(:,:)
  integer :: i, c, m           ! loop index
  integer :: ierr             ! I/O status
  type(io_desc_t), pointer :: iodesc
  type(var_desc_t)         :: vardesc
  integer                  :: csize, vsize
  character(len=4)         :: num
  integer                  :: dims(3), gdims(3), nhdims
  integer                  :: err_handling
  integer                  :: physgrid
  !-----------------------------------------------------------------------

  ! subcol_read_restart must be called before pbuf_read_restart
  if (is_subcol_on()) then
     call subcol_read_restart(File)
  end if

  csize=endchunk-begchunk+1
  dims(1) = pcols
  dims(2) = csize

  physgrid = cam_grid_id('physgrid')

  call cam_grid_dimensions(physgrid, gdims(1:2))

  if (gdims(2) == 1) then
     nhdims = 1
  else
     nhdims = 2
  end if
  call cam_grid_get_decomp(physgrid, dims(1:2), gdims(1:nhdims), pio_double, &
       iodesc)

  ! When we are ready: data for chemistry
  ! call chem_read_restart(File)

  ! call read_prescribed_ozone_restart(File)
  !  call read_prescribed_ghg_restart(File)
  !  call read_prescribed_aero_restart(File)
  !  call read_prescribed_volcaero_restart(File)

     allocate(tmpfield2(pcols, begchunk:endchunk))
     tmpfield2 = fillvalue

     ierr = pio_inq_varid(File, 'FLWDS', vardesc)
     call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
     do c=begchunk,endchunk
       do i=1,pcols
         cam_out(c)%flwds(i) = tmpfield2(i, c)
       end do
     end do

     ierr = pio_inq_varid(File, 'SOLS', vardesc)
     call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
     do c=begchunk,endchunk
       do i=1,pcols
         cam_out(c)%sols(i) = tmpfield2(i, c)
       end do
     end do

     ierr = pio_inq_varid(File, 'SOLL', vardesc)
     call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
     do c=begchunk,endchunk
       do i=1,pcols
         cam_out(c)%soll(i) = tmpfield2(i, c)
       end do
     end do

     ierr = pio_inq_varid(File, 'SOLSD', vardesc)
     call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
     do c=begchunk,endchunk
       do i=1,pcols
         cam_out(c)%solsd(i) = tmpfield2(i, c)
       end do
     end do

     ierr = pio_inq_varid(File, 'SOLLD', vardesc)
     call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
     do c=begchunk,endchunk
       do i=1,pcols
         cam_out(c)%solld(i) = tmpfield2(i, c)
       end do
     end do

     ierr = pio_inq_varid(File, 'BCPHIDRY', vardesc)
     call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
     do c=begchunk,endchunk
       do i=1,pcols
         cam_out(c)%bcphidry(i) = tmpfield2(i, c)
       end do
     end do

     ierr = pio_inq_varid(File, 'BCPHODRY', vardesc)
     call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
     do c=begchunk,endchunk
       do i=1,pcols
         cam_out(c)%bcphodry(i) = tmpfield2(i, c)
       end do
     end do

     ierr = pio_inq_varid(File, 'OCPHIDRY', vardesc)
     call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
     do c=begchunk,endchunk
       do i=1,pcols
         cam_out(c)%ocphidry(i) = tmpfield2(i, c)
       end do
     end do

     ierr = pio_inq_varid(File, 'OCPHODRY', vardesc)
     call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
     do c=begchunk,endchunk
       do i=1,pcols
         cam_out(c)%ocphodry(i) = tmpfield2(i, c)
       end do
     end do

     ierr = pio_inq_varid(File, 'DSTDRY1', vardesc)
     call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
     do c=begchunk,endchunk
       do i=1,pcols
         cam_out(c)%dstdry1(i) = tmpfield2(i, c)
       end do
     end do

     ierr = pio_inq_varid(File, 'DSTDRY2', vardesc)
     call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
     do c=begchunk,endchunk
       do i=1,pcols
         cam_out(c)%dstdry2(i) = tmpfield2(i, c)
       end do
     end do

     ierr = pio_inq_varid(File, 'DSTDRY3', vardesc)
     call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
     do c=begchunk,endchunk
       do i=1,pcols
         cam_out(c)%dstdry3(i) = tmpfield2(i, c)
       end do
     end do

     ierr = pio_inq_varid(File, 'DSTDRY4', vardesc)
     call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
     do c=begchunk,endchunk
       do i=1,pcols
         cam_out(c)%dstdry4(i) = tmpfield2(i, c)
       end do
     end do

     if (co2_transport()) then
       ierr = pio_inq_varid(File, 'CO2PROG', vardesc)
       call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
       do c=begchunk,endchunk
         do i=1,pcols
           cam_out(c)%co2prog(i) = tmpfield2(i, c)
         end do
       end do

       ierr = pio_inq_varid(File, 'CO2DIAG', vardesc)
       call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
       do c=begchunk,endchunk
         do i=1,pcols
           cam_out(c)%co2diag(i) = tmpfield2(i, c)
         end do
       end do
     end if

     ! Reading the CFLX* components from the restart is optional for
     ! backwards compatibility.
     ! These components are only needed if they are not handled by the
     ! coupling layer restart (the ".rs." file), and if the values are
     ! used in the tphysbc physics before the tphysac code has a chance
     ! to update the values that are coming from boundary datasets.
     do m = 1, pcnst

       write(num,'(i4.4)') m

       call pio_seterrorhandling(File, PIO_BCAST_ERROR, err_handling)
       ierr = pio_inq_varid(File, 'CFLX'//num, vardesc)
       call pio_seterrorhandling(File, err_handling)

       if (ierr == PIO_NOERR) then ! CFLX variable found on restart file
         call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
         do c= begchunk, endchunk
           do i = 1, pcols
             cam_in(c)%cflx(i,m) = tmpfield2(i, c)
           end do
         end do
       end if

     end do

     call pio_seterrorhandling(File, PIO_BCAST_ERROR, err_handling)
     ierr = pio_inq_varid(File, 'wsx', vardesc)
     if (ierr == PIO_NOERR) then ! variable found on restart file
       call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
       do c= begchunk, endchunk
         do i = 1, pcols
           cam_in(c)%wsx(i) = tmpfield2(i, c)
         end do
       end do
     end if
     ierr = pio_inq_varid(File, 'wsy', vardesc)
     if (ierr == PIO_NOERR) then ! variable found on restart file
       call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
       do c= begchunk, endchunk
         do i = 1, pcols
           cam_in(c)%wsy(i) = tmpfield2(i, c)
         end do
       end do
     end if
     ierr = pio_inq_varid(File, 'shf', vardesc)
     if (ierr == PIO_NOERR) then ! variable found on restart file
       call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
       do c= begchunk, endchunk
         do i = 1, pcols
           cam_in(c)%shf(i) = tmpfield2(i, c)
         end do
       end do
     endif
     call pio_seterrorhandling(File, err_handling)

     deallocate(tmpfield2)

     call radiation_read_restart(file)
!     call carma_restart_read(file)

end subroutine mars_restart_read
!============================================================================

!=======================================================================
  subroutine mars_grow_topo(state, pbuf, ztodt)
    !-----------------------------------------------------------------------
    ! Purpose: Scale topography fields (SGH, PHIS) from 0 to 1 over time
    !          to prevent initialization shock.
    !-----------------------------------------------------------------------
    use physics_types,  only: physics_state
    use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_get_index
    use time_manager,   only: get_curr_date, get_step_size, get_nstep
    use physconst,      only: gravit

    type(physics_state), intent(inout) :: state
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8), intent(in)               :: ztodt

    ! Local variables
    integer :: yr, mon, day, ncsec
    real(r8) :: time_days, scale_factor
    real(r8), pointer :: sgh(:,:), sgh30(:,:)
    integer :: sgh_idx, sgh30_idx
    integer :: i, k, dt, nstep

    ! 1. Calculate Scaling Factor
    if (mars_topo_ramp_days <= 0.0_r8) return ! Scaling disabled
    write(6,*)'In mars_grow_topo'

    call get_curr_date(yr, mon, day, ncsec)

    ! Calculate total simulation time in days (assuming start at day 0/1)
    ! Note: Adjust this logic if your run doesn't start at day 0.
    ! Simple approximation: Just rely on elapsed time steps if available,
    ! but here we approximate from date components or just use ztodt accumulation if stored.
    ! Ideally, pass 'nstep' or total time. For now, we assume simple start.

    ! Robust method: Use a global time tracker or assume start at time 0.
    ! Here we assume the user sets ramp_days relative to start of run.
    ! If this is a restart, we assume topography is already grown (scale=1).

    ! ** CRITICAL **: You likely need to pass 'nstep' to this routine or use time_manager
    ! 'get_curr_calday' if starting from day 0.
    ! For simplicity here, we assume standard startup.

    ! 1. Calculate Time from Steps (Robust across restarts/calendars)
    nstep = get_nstep()
    dt    = get_step_size()

    ! Convert seconds to days
    time_days = (real(nstep, r8) * dt) / 86400.0_r8

    if (time_days >= mars_topo_ramp_days) then
       scale_factor = 1.0_r8
    else
       scale_factor = max(0.0_r8, time_days / mars_topo_ramp_days)
    endif
    if (masterproc) write(6,*)'scale_factor=',scale_factor
    ! 2. Scale Physics State PHIS (Geopotential Height) Don't
    ! This affects surface pressure calc in physics and some drag routines
    ! state%phis(:) = state%phis(:) * scale_factor

    ! 3. Scale Physics Buffer Fields (SGH, SGH30)
    ! These drive the Gravity Wave Drag schemes directly.
    sgh_idx = pbuf_get_index('SGH')
    sgh30_idx = pbuf_get_index('SGH30')

    call pbuf_get_field(pbuf, sgh_idx, sgh)
    call pbuf_get_field(pbuf, sgh30_idx, sgh30)

    sgh(:,:)   = sgh(:,:)   * scale_factor
    sgh30(:,:) = sgh30(:,:) * scale_factor

    ! Optional: Print status periodically (e.g., column 1)
    if (masterproc .and. state%lchnk == begchunk) then
       write(iulog,*) 'MARS_TOPO_RAMP: Day=', time_days, ' Scale=', scale_factor
    endif

  end subroutine mars_grow_topo
end module mars_cam
