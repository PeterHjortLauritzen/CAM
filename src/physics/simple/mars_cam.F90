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
  use cam_history,    only: outfld
  use cam_logfile,    only: iulog
  use camsrfexch,     only: cam_in_t,cam_out_t
  use constituents,   only: pcnst
  use physconst,      only: cpair
  use physics_buffer, only: dtype_r8, pbuf_add_field, physics_buffer_desc
  use physics_types,  only: physics_ptend_init
  use physics_types,  only: physics_state, physics_ptend
  use pio             ,only: file_desc_t, var_desc_t, io_desc_t, pio_double, pio_def_var
  use pio             ,only: pio_write_darray, pio_read_darray, pio_inq_varid
  use ppgrid,         only: pcols, pver, pverp, begchunk, endchunk
  use shr_const_mod,   only: SHR_CONST_STEBOL, SHR_CONST_REARTH, SHR_CONST_KARMAN, SHR_CONST_TKTRIP
  use shr_const_mod,  only: pi => shr_const_pi
  use shr_kind_mod,   only: r8 => shr_kind_r8, cl=>shr_kind_cl
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
    use physconst,      only: cpair, mwh2o
    use constituents,   only: cnst_add
    use physics_buffer, only: pbuf_add_field, dtype_r8
    use radiation,      only: radiation_register

    integer                :: mm

    call radiation_register()

    call cnst_add('CO2', 44._r8, 800._r8, 1.e-12_r8, mm, fixed_ubc=.false., &
         longname='CO2', readiv=.true., is_convtran1=.true.)

    call cnst_add('N2', 28._r8, 800._r8, 1.e-12_r8, mm, fixed_ubc=.false., &
         longname='N2', readiv=.true., is_convtran1=.true.)

  end subroutine mars_register
  !==============================================================================


  !==============================================================================
  subroutine mars_readnl(nlfile)
    !
    ! mars_readnl: Read in parameters controlling Mars parameterizations.
    !=====================================================================
    use namelist_utils,only: find_group_name
    use units         ,only: getunit, freeunit
    use radiation,     only: radiation_readnl
    use rad_constituents,    only: rad_cnst_readnl

    ! Input Parameters
    !------------------
    character(len=*),intent(in):: nlfile
    !
    ! Local Values
    !--------------
    integer:: ierr,unitn

    character(len=*), parameter :: sub = 'mars_readnl'

!!$    namelist /mars_nl/ mars_Ts
!!$
!!$   ! Read in namelist values
!!$    !-------------------------
!!$    if(masterproc) then
!!$      unitn = getunit()
!!$      open(unitn,file=trim(nlfile),status='old')
!!$      call find_group_name(unitn,'mars_nl',status=ierr)
!!$      if(ierr == 0) then
!!$        read(unitn,mars_nl,iostat=ierr)
!!$        if(ierr /= 0) then
!!$          call endrun(sub//': ERROR reading namelist')
!!$        endif
!!$      endif
!!$      close(unitn)
!!$      call freeunit(unitn)
!!$    endif
!!$
!!$    ! Broadcast namelist values
!!$    !---------------------------
!!$    call mpi_bcast(mars_Ts  , 1, mpi_real8 , mstrid, mpicom, ierr)
!!$    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: mars_Wind_min")

    call radiation_readnl(nlfile)
    call rad_cnst_readnl(nlfile)

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
    use rad_constituents,    only: rad_cnst_init
    !
    ! Input Parameters
    !------------------
    type(physics_state)      ,pointer:: phys_state(:)
    type(physics_buffer_desc),pointer:: pbuf2d    (:,:)
    !
    ! Local Values
    !---------------

    ! For now just initialize radiation

    call radiation_init(pbuf2d)
    call rad_cnst_init()

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

subroutine mars_surface_init(cam_in)
  !
  !
  !==========================================================================
  !
  ! Passed variables
  !--------------------
  type(cam_in_t),     intent(inout)          :: cam_in
  !
  ! Local values
  !--------------
  cam_in%ts(:)=190.0_r8
  cam_in%asdir(:)=.5_r8
  cam_in%asdif(:)=.5_r8
  cam_in%aldir(:)=.5_r8
  cam_in%aldif(:)=.5_r8

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

  ! THIS IS A TEMPORARY TO SET SURFACE VALUES FOR TESTING - REMOVE WHEN HAVE SURFACE COMPONENT
  call mars_surface_init(cam_in)

  call radiation_tend( state, ptend, pbuf, cam_out, cam_in, net_flx)

end subroutine mars_radiative_tend

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


end module mars_cam
