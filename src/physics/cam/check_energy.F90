
module check_energy

!---------------------------------------------------------------------------------
! Purpose:
!
! Module to check
!   1. vertically integrated total energy and water conservation for each
!      column within the physical parameterizations
!
!   2. global mean total energy conservation between the physics output state
!      and the input state on the next time step.
!
!   3. add a globally uniform heating term to account for any change of total energy in 2.
  !
! Author: Byron Boville  Oct 31, 2002
!
! Modifications:
!   03.03.29  Boville  Add global energy check and fixer.
!
!---------------------------------------------------------------------------------

  use shr_kind_mod,    only: r8 => shr_kind_r8
  use ppgrid,          only: pcols, pver, begchunk, endchunk
  use spmd_utils,      only: masterproc

  use gmean_mod,       only: gmean
  use physconst,       only: gravit, rga, latvap, latice, cpair, rair
  use air_composition, only: cpairv, rairv, cp_or_cv_dycore
  use physics_types,   only: physics_state, physics_tend, physics_ptend, physics_ptend_init
  use constituents,    only: cnst_get_ind, pcnst, cnst_name, cnst_get_type_byind
  use time_manager,    only: is_first_step
  use cam_logfile,     only: iulog

  implicit none
  private

! Public types:
  public check_tracers_data

! Public methods
  public :: check_energy_readnl    ! read namelist values
  public :: check_energy_register  ! register fields in physics buffer
  public :: check_energy_get_integrals ! get energy integrals computed in check_energy_gmean
  public :: check_energy_init      ! initialization of module
  public :: check_energy_timestep_init  ! timestep initialization of energy integrals and cumulative boundary fluxes
  public :: check_energy_chng      ! check changes in integrals against cumulative boundary fluxes
  public :: check_energy_gmean     ! global means of physics input and output total energy
  public :: check_energy_fix       ! add global mean energy difference as a heating
  public :: check_tracers_init      ! initialize tracer integrals and cumulative boundary fluxes
  public :: check_tracers_chng      ! check changes in integrals against cumulative boundary fluxes

  public :: tot_energy_phys ! calculate and output total energy and axial angular momentum diagnostics
  public :: pressure_enthalpy_adjustment
! Private module data

  logical  :: print_energy_errors = .false.

  real(r8) :: teout_glob           ! global mean energy of output state
  real(r8) :: teinp_glob           ! global mean energy of input state
  real(r8) :: tedif_glob           ! global mean energy difference
  real(r8) :: psurf_glob           ! global mean surface pressure
  real(r8) :: ptopb_glob           ! global mean top boundary pressure
  real(r8) :: heat_glob            ! global mean heating rate

! Physics buffer indices

  integer  :: teout_idx  = 0       ! teout index in physics buffer
  integer  :: dtcore_idx = 0       ! dtcore index in physics buffer
  integer  :: dqcore_idx = 0       ! dqcore index in physics buffer
  integer  :: ducore_idx = 0       ! ducore index in physics buffer
  integer  :: dvcore_idx = 0       ! dvcore index in physics buffer

  type check_tracers_data
     real(r8) :: tracer(pcols,pcnst)       ! initial vertically integrated total (kinetic + static) energy
     real(r8) :: tracer_tnd(pcols,pcnst)   ! cumulative boundary flux of total energy
     integer :: count(pcnst)               ! count of values with significant imbalances
  end type check_tracers_data


!===============================================================================
contains
!===============================================================================

subroutine check_energy_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_logical
   use cam_abortutils,  only: endrun

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: sub = 'check_energy_readnl'

   namelist /check_energy_nl/ print_energy_errors
   !-----------------------------------------------------------------------------

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'check_energy_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, check_energy_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(sub//': FATAL: reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

   call mpi_bcast(print_energy_errors, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: print_energy_errors")

   if (masterproc) then
      write(iulog,*) 'check_energy options:'
      write(iulog,*) '  print_energy_errors =', print_energy_errors
   end if

end subroutine check_energy_readnl

!===============================================================================

  subroutine check_energy_register()
!
! Register fields in the physics buffer.
!
!-----------------------------------------------------------------------

    use physics_buffer, only : pbuf_add_field, dtype_r8, dyn_time_lvls
    use physics_buffer, only : pbuf_register_subcol
    use subcol_utils,   only : is_subcol_on

!-----------------------------------------------------------------------

! Request physics buffer space for fields that persist across timesteps.

    call pbuf_add_field('TEOUT', 'global',dtype_r8 , (/pcols,dyn_time_lvls/),      teout_idx)
    call pbuf_add_field('DTCORE','global',dtype_r8,  (/pcols,pver,dyn_time_lvls/),dtcore_idx)
    ! DQCORE refers to dycore tendency of water vapor
    call pbuf_add_field('DQCORE','global',dtype_r8,  (/pcols,pver,dyn_time_lvls/),dqcore_idx)
    call pbuf_add_field('DUCORE','global',dtype_r8,  (/pcols,pver,dyn_time_lvls/),ducore_idx)
    call pbuf_add_field('DVCORE','global',dtype_r8,  (/pcols,pver,dyn_time_lvls/),dvcore_idx)
    if(is_subcol_on()) then
      call pbuf_register_subcol('TEOUT', 'phys_register', teout_idx)
      call pbuf_register_subcol('DTCORE', 'phys_register', dtcore_idx)
      call pbuf_register_subcol('DQCORE', 'phys_register', dqcore_idx)
      call pbuf_register_subcol('DUCORE', 'phys_register', ducore_idx)
      call pbuf_register_subcol('DVCORE', 'phys_register', dvcore_idx)
    end if

  end subroutine check_energy_register

!===============================================================================

subroutine check_energy_get_integrals( tedif_glob_out, heat_glob_out )

!-----------------------------------------------------------------------
! Purpose: Return energy integrals
!-----------------------------------------------------------------------

     real(r8), intent(out), optional :: tedif_glob_out
     real(r8), intent(out), optional :: heat_glob_out

!-----------------------------------------------------------------------

   if ( present(tedif_glob_out) ) then
      tedif_glob_out = tedif_glob
   endif
   if ( present(heat_glob_out) ) then
      heat_glob_out = heat_glob
   endif

end subroutine check_energy_get_integrals

!================================================================================================

  subroutine check_energy_init()
!
! Initialize the energy conservation module
!
!-----------------------------------------------------------------------
    use cam_history,       only: addfld, add_default, horiz_only
    use phys_control,      only: phys_getopts

    implicit none

    logical          :: history_budget, history_waccm
    integer          :: history_budget_histfile_num ! output history file number for budget fields

!-----------------------------------------------------------------------

    call phys_getopts( history_budget_out = history_budget, &
                       history_budget_histfile_num_out = history_budget_histfile_num, &
                       history_waccm_out = history_waccm )

! register history variables
    call addfld('TEINP',  horiz_only,  'A', 'J/m2', 'Total energy of physics input')
    call addfld('TEOUT',  horiz_only,  'A', 'J/m2', 'Total energy of physics output')
    call addfld('TEFIX',  horiz_only,  'A', 'J/m2', 'Total energy after fixer')
    call addfld('EFIX',   horiz_only,  'A', 'W/m2', 'Effective sensible heat flux due to energy fixer')
    call addfld('DTCORE', (/ 'lev' /), 'A', 'K/s' , 'T tendency due to dynamical core')
    call addfld('DQCORE', (/ 'lev' /), 'A', 'kg/kg/s' , 'Water vapor tendency due to dynamical core')

    if ( history_budget ) then
       call add_default ('DTCORE', history_budget_histfile_num, ' ')
    end if
    if ( history_waccm ) then
       call add_default ('DTCORE', 1, ' ')
    end if

  end subroutine check_energy_init

!===============================================================================

  subroutine check_energy_timestep_init(state, tend, pbuf, col_type)
    use cam_thermo,      only: get_hydrostatic_energy
    use physics_buffer,  only: physics_buffer_desc, pbuf_set_field
    use cam_abortutils,  only: endrun
    use dyn_tests_utils, only: vc_physics, vc_dycore, vc_height, vc_dry_pressure
    use physics_types,   only: phys_te_idx, dyn_te_idx
!-----------------------------------------------------------------------
! Compute initial values of energy and water integrals,
! zero cumulative tendencies
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    type(physics_state),   intent(inout)    :: state
    type(physics_tend ),   intent(inout)    :: tend
    type(physics_buffer_desc), pointer      :: pbuf(:)
    integer, optional                       :: col_type  ! Flag inidicating whether using grid or subcolumns
!---------------------------Local storage-------------------------------
    real(r8)              :: cp_or_cv(state%psetcols,pver)
    integer lchnk                                  ! chunk identifier
    integer ncol                                   ! number of atmospheric columns
!-----------------------------------------------------------------------

    lchnk = state%lchnk
    ncol  = state%ncol

    ! cp_or_cv needs to be allocated to a size which matches state and ptend
    ! If psetcols == pcols, cpairv is the correct size and just copy into cp_or_cv
    ! If psetcols > pcols and all cpairv match cpair, then assign the constant cpair

    if (state%psetcols == pcols) then
       cp_or_cv(:,:) = cpairv(:,:,lchnk)
    else if (state%psetcols > pcols .and. all(cpairv(:,:,lchnk) == cpair)) then
       cp_or_cv(1:ncol,:) = cpair
    else
       call endrun('check_energy_timestep_init: cpairv is not allowed to vary when subcolumns are turned on')
    end if
    !
    ! CAM physics total energy
    !
    call get_hydrostatic_energy(state%q(1:ncol,1:pver,1:pcnst),.true.,               &
         state%pdel(1:ncol,1:pver), cp_or_cv(1:ncol,1:pver),                         &
         state%u(1:ncol,1:pver), state%v(1:ncol,1:pver), state%T(1:ncol,1:pver),     &
         vc_physics, ptop=state%pintdry(1:ncol,1), phis = state%phis(1:ncol),&
         te = state%te_ini(1:ncol,phys_te_idx), H2O = state%tw_ini(1:ncol,phys_te_idx))
    !
    ! Dynamical core total energy
    !
    state%temp_ini(:ncol,:) = state%T(:ncol,:)
    state%z_ini(:ncol,:)    = state%zm(:ncol,:)
    if (vc_dycore == vc_height) then
      !
      ! MPAS specific hydrostatic energy computation (internal energy)
      !
      if (state%psetcols == pcols) then
        cp_or_cv(:ncol,:) = cp_or_cv_dycore(:ncol,:,lchnk)
      else
        cp_or_cv(:ncol,:) = cpair-rair
      endif

      call get_hydrostatic_energy(state%q(1:ncol,1:pver,1:pcnst),.true.,               &
           state%pdel(1:ncol,1:pver), cp_or_cv(1:ncol,1:pver),                         &
           state%u(1:ncol,1:pver), state%v(1:ncol,1:pver), state%T(1:ncol,1:pver),     &
           vc_dycore, ptop=state%pintdry(1:ncol,1), phis = state%phis(1:ncol),         &
           z_mid = state%z_ini(1:ncol,:),                                              &
           te = state%te_ini(1:ncol,dyn_te_idx), H2O = state%tw_ini(1:ncol,dyn_te_idx))
    else if (vc_dycore == vc_dry_pressure) then
      !
      ! SE specific hydrostatic energy (enthalpy)
      !
      if (state%psetcols == pcols) then
        cp_or_cv(:ncol,:) = cp_or_cv_dycore(:ncol,:,lchnk)
      else
        cp_or_cv(:ncol,:) = cpair
      endif
      call get_hydrostatic_energy(state%q(1:ncol,1:pver,1:pcnst),.true.,               &
           state%pdel(1:ncol,1:pver), cp_or_cv(1:ncol,1:pver),                         &
           state%u(1:ncol,1:pver), state%v(1:ncol,1:pver), state%T(1:ncol,1:pver),     &
           vc_dry_pressure, ptop=state%pintdry(1:ncol,1), phis = state%phis(1:ncol),   &
           te = state%te_ini(1:ncol,dyn_te_idx), H2O = state%tw_ini(1:ncol,dyn_te_idx))
    else
      !
      ! dycore energy is the same as physics
      !
      state%te_ini(1:ncol,dyn_te_idx) = state%te_ini(1:ncol,phys_te_idx)
      state%tw_ini(1:ncol,dyn_te_idx) = state%tw_ini(1:ncol,phys_te_idx)
    end if
    state%te_cur(:ncol,:) = state%te_ini(:ncol,:)
    state%tw_cur(:ncol,:) = state%tw_ini(:ncol,:)

! zero cummulative boundary fluxes
    tend%te_tnd(:ncol) = 0._r8
    tend%te_sen(:ncol) = 0._r8
    tend%te_lat(:ncol) = 0._r8
    !xxx should have enthalpy flux here
    tend%tw_tnd(:ncol) = 0._r8

    state%count = 0

! initialize physics buffer
    if (is_first_step()) then
       call pbuf_set_field(pbuf, teout_idx, state%te_ini(:,dyn_te_idx), col_type=col_type)
    end if

  end subroutine check_energy_timestep_init

!===============================================================================

  subroutine check_energy_chng(state, tend, name, nstep, ztodt,        &
       flx_vap, flx_cnd, flx_ice, flx_sen)
    use cam_thermo,      only: get_hydrostatic_energy
    use dyn_tests_utils, only: vc_physics, vc_dycore, vc_height, vc_dry_pressure
    use cam_abortutils,  only: endrun
    use physics_types,   only: phys_te_idx, dyn_te_idx
!-----------------------------------------------------------------------
! Check that the energy and water change matches the boundary fluxes
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    type(physics_state)    , intent(inout) :: state
    type(physics_tend )    , intent(inout) :: tend
    character*(*),intent(in) :: name               ! parameterization name for fluxes
    integer , intent(in   ) :: nstep               ! current timestep number
    real(r8), intent(in   ) :: ztodt               ! 2 delta t (model time increment)
    real(r8), intent(in   ) :: flx_vap(:)          ! (pcols) - boundary flux of vapor         (kg/m2/s)
    real(r8), intent(in   ) :: flx_cnd(:)          ! (pcols) -boundary flux of liquid+ice    (m/s) (precip?)
    real(r8), intent(in   ) :: flx_ice(:)          ! (pcols) -boundary flux of ice           (m/s) (snow?)
    real(r8), intent(in   ) :: flx_sen(:)          ! (pcols) -boundary flux of sensible heat (w/m2)

!******************** BAB ******************************************************
!******* Note that the precip and ice fluxes are in precip units (m/s). ********
!******* I would prefer to have kg/m2/s.                                ********
!******* I would also prefer liquid (not total) and ice fluxes          ********
!*******************************************************************************

!---------------------------Local storage-------------------------------

    real(r8) :: te_xpd(state%ncol)                 ! expected value (f0 + dt*boundary_flux)
    real(r8) :: te_dif(state%ncol)                 ! energy of input state - original energy
    real(r8) :: te_tnd(state%ncol)                 ! tendency from last process
    real(r8) :: te_rer(state%ncol)                 ! relative error in energy column

    real(r8) :: tw_xpd(state%ncol)                 ! expected value (w0 + dt*boundary_flux)
    real(r8) :: tw_dif(state%ncol)                 ! tw_inp - original water
    real(r8) :: tw_tnd(state%ncol)                 ! tendency from last process
    real(r8) :: tw_rer(state%ncol)                 ! relative error in water column

    real(r8) :: te(state%ncol)                     ! vertical integral of total energy
    real(r8) :: tw(state%ncol)                     ! vertical integral of total water
    real(r8) :: cp_or_cv(state%psetcols,pver)      ! cp or cv depending on vcoord
    real(r8) :: scaling(state%psetcols,pver)       ! scaling for conversion of temperature increment
    real(r8) :: temp(state%ncol,pver)              ! temperature

    real(r8) :: se(state%ncol)                     ! enthalpy or internal energy (J/m2)
    real(r8) :: po(state%ncol)                     ! surface potential or potential energy (J/m2)
    real(r8) :: ke(state%ncol)                     ! kinetic energy    (J/m2)
    real(r8) :: wv(state%ncol)                     ! column integrated vapor       (kg/m2)
    real(r8) :: liq(state%ncol)                    ! column integrated liquid      (kg/m2)
    real(r8) :: ice(state%ncol)                    ! column integrated ice         (kg/m2)

    integer lchnk                                  ! chunk identifier
    integer ncol                                   ! number of atmospheric columns
    integer  i                                     ! column index
!-----------------------------------------------------------------------

    lchnk = state%lchnk
    ncol  = state%ncol

    ! If psetcols == pcols, cpairv is the correct size and just copy into cp_or_cv
    ! If psetcols > pcols and all cpairv match cpair, then assign the constant cpair

    if (state%psetcols == pcols) then
       cp_or_cv(:,:) = cpairv(:,:,lchnk)
    else if (state%psetcols > pcols .and. all(cpairv(:,:,:) == cpair)) then
       cp_or_cv(:,:) = cpair
    else
       call endrun('check_energy_chng: cpairv is not allowed to vary when subcolumns are turned on')
    end if

    call get_hydrostatic_energy(state%q(1:ncol,1:pver,1:pcnst),.true.,               &
         state%pdel(1:ncol,1:pver), cp_or_cv(1:ncol,1:pver),                         &
         state%u(1:ncol,1:pver), state%v(1:ncol,1:pver), state%T(1:ncol,1:pver),     &
         vc_physics, ptop=state%pintdry(1:ncol,1), phis = state%phis(1:ncol),        &
         te = te(1:ncol), H2O = tw(1:ncol), se=se(1:ncol),po=po(1:ncol),             &
         ke=ke(1:ncol),wv=wv(1:ncol),liq=liq(1:ncol),ice=ice(1:ncol))
    ! compute expected values and tendencies
    do i = 1, ncol
       ! change in static energy and total water
       te_dif(i) = te(i) - state%te_cur(i,phys_te_idx)
       tw_dif(i) = tw(i) - state%tw_cur(i,phys_te_idx)

       ! expected tendencies from boundary fluxes for last process
       te_tnd(i) = flx_vap(i)*(latvap+latice) - (flx_cnd(i) - flx_ice(i))*1000._r8*latice + flx_sen(i)
       tw_tnd(i) = flx_vap(i) - flx_cnd(i) *1000._r8

       ! cummulative tendencies from boundary fluxes
       tend%te_tnd(i) = tend%te_tnd(i) + te_tnd(i)
       tend%tw_tnd(i) = tend%tw_tnd(i) + tw_tnd(i)
       tend%te_sen(i) = tend%te_sen(i)+flx_sen(i)
       tend%te_lat(i) = tend%te_lat(i)+flx_vap(i)*(latvap+latice)  - (flx_cnd(i) - flx_ice(i))*1000._r8*latice
       ! expected new values from previous state plus boundary fluxes
       te_xpd(i) = state%te_cur(i,phys_te_idx) + te_tnd(i)*ztodt
       tw_xpd(i) = state%tw_cur(i,phys_te_idx) + tw_tnd(i)*ztodt

       ! relative error, expected value - input state / previous state
       te_rer(i) = (te_xpd(i) - te(i)) / state%te_cur(i,phys_te_idx)
    end do
    ! relative error for total water (allow for dry atmosphere)
    tw_rer = 0._r8
    where (state%tw_cur(:ncol,phys_te_idx) > 0._r8)
       tw_rer(:ncol) = (tw_xpd(:ncol) - tw(:ncol)) / state%tw_cur(:ncol,1)
    end where

    ! error checking
    if (print_energy_errors) then
       if (any(abs(te_rer(1:ncol)) > 1.E-14_r8 .or. abs(tw_rer(1:ncol)) > 1.E-10_r8)) then
          do i = 1, ncol
             ! the relative error threshold for the water budget has been reduced to 1.e-10
             ! to avoid messages generated by QNEG3 calls
             ! PJR- change to identify if error in energy or water
             if (abs(te_rer(i)) > 1.E-14_r8 ) then
                state%count = state%count + 1
                write(iulog,*) "significant energy conservation error after ", name,        &
                      " count", state%count, " nstep", nstep, "chunk", lchnk, "col", i
                write(iulog,*) te(i),te_xpd(i),te_dif(i),tend%te_tnd(i)*ztodt,  &
                      te_tnd(i)*ztodt,te_rer(i)
             endif
             if ( abs(tw_rer(i)) > 1.E-10_r8) then
                state%count = state%count + 1
                write(iulog,*) "significant water conservation error after ", name,        &
                      " count", state%count, " nstep", nstep, "chunk", lchnk, "col", i
                write(iulog,*) tw(i),tw_xpd(i),tw_dif(i),tend%tw_tnd(i)*ztodt,  &
                      tw_tnd(i)*ztodt,tw_rer(i)
             end if
          end do
       end if
    end if

    ! copy new value to state

    do i = 1, ncol
      state%te_cur(i,phys_te_idx) = te(i)
      state%tw_cur(i,phys_te_idx) = tw(i)
    end do

    !
    ! Dynamical core total energy
    !
    if (vc_dycore == vc_height) then
      !
      ! compute cv if vertical coordinate is height: cv = cp - R
      !
      ! Note: cp_or_cv set above for pressure coordinate
      if (state%psetcols == pcols) then
        cp_or_cv(:ncol,:) = cp_or_cv_dycore(:ncol,:,lchnk)
      else
        cp_or_cv(:ncol,:) = cpair-rair
      endif
      scaling(:,:)   = cpairv(:,:,lchnk)/cp_or_cv(:,:) !cp/cv scaling
      temp(1:ncol,:) = state%temp_ini(1:ncol,:)+scaling(1:ncol,:)*(state%T(1:ncol,:)-state%temp_ini(1:ncol,:))
      call get_hydrostatic_energy(state%q(1:ncol,1:pver,1:pcnst),.true.,               &
           state%pdel(1:ncol,1:pver), cp_or_cv(1:ncol,1:pver),                         &
           state%u(1:ncol,1:pver), state%v(1:ncol,1:pver), temp(1:ncol,1:pver),        &
           vc_dycore, ptop=state%pintdry(1:ncol,1), phis = state%phis(1:ncol),         &
           z_mid = state%z_ini(1:ncol,:),                                              &
           te = state%te_cur(1:ncol,dyn_te_idx), H2O = state%tw_cur(1:ncol,dyn_te_idx))
    else if (vc_dycore == vc_dry_pressure) then
      !
      ! SE specific hydrostatic energy
      !
      if (state%psetcols == pcols) then
        cp_or_cv(:ncol,:) = cp_or_cv_dycore(:ncol,:,lchnk)
      else
        cp_or_cv(:ncol,:) = cpair
      endif
      !
      ! enthalpy scaling for energy consistency
      !
      scaling(:ncol,:) = cpairv(:ncol,:,lchnk)/cp_or_cv_dycore(:ncol,:,lchnk)
      temp(1:ncol,:)   = state%temp_ini(1:ncol,:)+scaling(1:ncol,:)*(state%T(1:ncol,:)-state%temp_ini(1:ncol,:))
      call get_hydrostatic_energy(state%q(1:ncol,1:pver,1:pcnst),.true.,               &
           state%pdel(1:ncol,1:pver), cp_or_cv(1:ncol,1:pver),                         &
           state%u(1:ncol,1:pver), state%v(1:ncol,1:pver), temp(1:ncol,1:pver),        &
           vc_dry_pressure, ptop=state%pintdry(1:ncol,1), phis = state%phis(1:ncol),   &
           te = state%te_cur(1:ncol,dyn_te_idx), H2O = state%tw_cur(1:ncol,dyn_te_idx))
    else
      state%te_cur(1:ncol,dyn_te_idx) = te(1:ncol)
      state%tw_cur(1:ncol,dyn_te_idx) = tw(1:ncol)
    end if
  end subroutine check_energy_chng


  subroutine check_energy_gmean(state, pbuf2d, dtime, nstep)

    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_get_chunk
    use physics_types,   only: dyn_te_idx
!-----------------------------------------------------------------------
! Compute global mean total energy of physics input and output states
! computed consistently with dynamical core vertical coordinate
! (under hydrostatic assumption)
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    type(physics_state), intent(in   ), dimension(begchunk:endchunk) :: state
    type(physics_buffer_desc),    pointer    :: pbuf2d(:,:)

    real(r8), intent(in) :: dtime        ! physics time step
    integer , intent(in) :: nstep        ! current timestep number

!---------------------------Local storage-------------------------------
    integer :: ncol                      ! number of active columns
    integer :: lchnk                     ! chunk index

    real(r8) :: te(pcols,begchunk:endchunk,4)
                                         ! total energy of input/output states (copy)
    real(r8) :: te_glob(4)               ! global means of total energy
    real(r8), pointer :: teout(:)
!-----------------------------------------------------------------------

    ! Copy total energy out of input and output states
    do lchnk = begchunk, endchunk
       ncol = state(lchnk)%ncol
       ! input energy using dynamical core energy formula
       te(:ncol,lchnk,1) = state(lchnk)%te_ini(:ncol,dyn_te_idx)
       ! output energy
       call pbuf_get_field(pbuf_get_chunk(pbuf2d,lchnk),teout_idx, teout)

       te(:ncol,lchnk,2) = teout(1:ncol)
       ! surface pressure for heating rate
       te(:ncol,lchnk,3) = state(lchnk)%pint(:ncol,pver+1)
       ! model top pressure for heating rate (not constant for z-based vertical coordinate!)
       te(:ncol,lchnk,4) = state(lchnk)%pint(:ncol,1)
    end do

    ! Compute global means of input and output energies and of
    ! surface pressure for heating rate (assume uniform ptop)
    call gmean(te, te_glob, 4)

    if (begchunk .le. endchunk) then
       teinp_glob = te_glob(1)
       teout_glob = te_glob(2)
       psurf_glob = te_glob(3)
       ptopb_glob = te_glob(4)

       ! Global mean total energy difference
       tedif_glob =  teinp_glob - teout_glob
       heat_glob  = -tedif_glob/dtime * gravit / (psurf_glob - ptopb_glob)
       if (masterproc) then
          write(iulog,'(1x,a9,1x,i8,5(1x,e25.17))') "nstep, te", nstep, teinp_glob, teout_glob, &
               heat_glob, psurf_glob, ptopb_glob
       end if
    else
       heat_glob = 0._r8
    end if  !  (begchunk .le. endchunk)

  end subroutine check_energy_gmean

!===============================================================================
  subroutine check_energy_fix(state, ptend, nstep, eshflx)

!-----------------------------------------------------------------------
! Add heating rate required for global mean total energy conservation
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    type(physics_state), intent(in   ) :: state
    type(physics_ptend), intent(out)   :: ptend

    integer , intent(in   ) :: nstep          ! time step number
    real(r8), intent(out  ) :: eshflx(pcols)  ! effective sensible heat flux

!---------------------------Local storage-------------------------------
    integer  :: i                        ! column
    integer  :: ncol                     ! number of atmospheric columns in chunk
!-----------------------------------------------------------------------
    ncol = state%ncol

    call physics_ptend_init(ptend, state%psetcols, 'chkenergyfix', ls=.true.)

#if ( defined OFFLINE_DYN )
    ! disable the energy fix for offline driver
    heat_glob = 0._r8
#endif
! add (-) global mean total energy difference as heating
    ptend%s(:ncol,:pver) = heat_glob

! compute effective sensible heat flux
    do i = 1, ncol
       eshflx(i) = heat_glob * (state%pint(i,pver+1) - state%pint(i,1)) * rga
    end do

    return
  end subroutine check_energy_fix


!===============================================================================
  subroutine check_tracers_init(state, tracerint)

!-----------------------------------------------------------------------
! Compute initial values of tracers integrals,
! zero cumulative tendencies
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------

    type(physics_state),   intent(in)    :: state
    type(check_tracers_data), intent(out)   :: tracerint

!---------------------------Local storage-------------------------------

    real(r8) :: tr(pcols)                          ! vertical integral of tracer
    real(r8) :: trpdel(pcols, pver)                ! pdel for tracer

    integer ncol                                   ! number of atmospheric columns
    integer  i,k,m                                 ! column, level,constituent indices
    integer :: ixcldice, ixcldliq                  ! CLDICE and CLDLIQ indices
    integer :: ixrain, ixsnow                      ! RAINQM and SNOWQM indices
    integer :: ixgrau                              ! GRAUQM index
!-----------------------------------------------------------------------

    ncol  = state%ncol
    call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
    call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
    call cnst_get_ind('RAINQM', ixrain,   abort=.false.)
    call cnst_get_ind('SNOWQM', ixsnow,   abort=.false.)
    call cnst_get_ind('GRAUQM', ixgrau,   abort=.false.)


    do m = 1,pcnst

       if ( any(m == (/ 1, ixcldliq, ixcldice, &
                           ixrain,   ixsnow, ixgrau /)) ) exit   ! dont process water substances
                                                                 ! they are checked in check_energy

       if (cnst_get_type_byind(m).eq.'dry') then
          trpdel(:ncol,:) = state%pdeldry(:ncol,:)
       else
          trpdel(:ncol,:) = state%pdel(:ncol,:)
       endif

       ! Compute vertical integrals of tracer
       tr = 0._r8
       do k = 1, pver
          do i = 1, ncol
             tr(i) = tr(i) + state%q(i,k,m)*trpdel(i,k)*rga
          end do
       end do

       ! Compute vertical integrals of frozen static tracers and total water.
       do i = 1, ncol
          tracerint%tracer(i,m) = tr(i)
       end do

       ! zero cummulative boundary fluxes
       tracerint%tracer_tnd(:ncol,m) = 0._r8

       tracerint%count(m) = 0

    end do

    return
  end subroutine check_tracers_init

!===============================================================================
  subroutine check_tracers_chng(state, tracerint, name, nstep, ztodt, cflx)

!-----------------------------------------------------------------------
! Check that the tracers and water change matches the boundary fluxes
! these checks are not save when there are tracers transformations, as
! they only check to see whether a mass change in the column is
! associated with a flux
!-----------------------------------------------------------------------

    use cam_abortutils, only: endrun


    implicit none

!------------------------------Arguments--------------------------------

    type(physics_state)    , intent(in   ) :: state
    type(check_tracers_data), intent(inout) :: tracerint! tracers integrals and boundary fluxes
    character*(*),intent(in) :: name               ! parameterization name for fluxes
    integer , intent(in   ) :: nstep               ! current timestep number
    real(r8), intent(in   ) :: ztodt               ! 2 delta t (model time increment)
    real(r8), intent(in   ) :: cflx(pcols,pcnst)       ! boundary flux of tracers       (kg/m2/s)

!---------------------------Local storage-------------------------------

    real(r8) :: tracer_inp(pcols,pcnst)                   ! total tracer of new (input) state
    real(r8) :: tracer_xpd(pcols,pcnst)                   ! expected value (w0 + dt*boundary_flux)
    real(r8) :: tracer_dif(pcols,pcnst)                   ! tracer_inp - original tracer
    real(r8) :: tracer_tnd(pcols,pcnst)                   ! tendency from last process
    real(r8) :: tracer_rer(pcols,pcnst)                   ! relative error in tracer column

    real(r8) :: tr(pcols)                           ! vertical integral of tracer
    real(r8) :: trpdel(pcols, pver)                       ! pdel for tracer

    integer lchnk                                  ! chunk identifier
    integer ncol                                   ! number of atmospheric columns
    integer  i,k                                   ! column, level indices
    integer :: ixcldice, ixcldliq                  ! CLDICE and CLDLIQ indices
    integer :: ixrain, ixsnow                      ! RAINQM and SNOWQM indices
    integer :: ixgrau                              ! GRAUQM index
    integer :: m                            ! tracer index
    character(len=8) :: tracname   ! tracername
!-----------------------------------------------------------------------

    lchnk = state%lchnk
    ncol  = state%ncol
    call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
    call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
    call cnst_get_ind('RAINQM', ixrain,   abort=.false.)
    call cnst_get_ind('SNOWQM', ixsnow,   abort=.false.)
    call cnst_get_ind('GRAUQM', ixgrau,   abort=.false.)

    do m = 1,pcnst

       if ( any(m == (/ 1, ixcldliq, ixcldice, &
                           ixrain,   ixsnow, ixgrau /)) ) exit   ! dont process water substances
                                                                 ! they are checked in check_energy
       tracname = cnst_name(m)
       if (cnst_get_type_byind(m).eq.'dry') then
          trpdel(:ncol,:) = state%pdeldry(:ncol,:)
       else
          trpdel(:ncol,:) = state%pdel(:ncol,:)
       endif

       ! Compute vertical integrals tracers
       tr = 0._r8
       do k = 1, pver
          do i = 1, ncol
             tr(i) = tr(i) + state%q(i,k,m)*trpdel(i,k)*rga
          end do
       end do

       ! Compute vertical integrals of tracer
       do i = 1, ncol
          tracer_inp(i,m) = tr(i)
       end do

       ! compute expected values and tendencies
       do i = 1, ncol
          ! change in tracers
          tracer_dif(i,m) = tracer_inp(i,m) - tracerint%tracer(i,m)

          ! expected tendencies from boundary fluxes for last process
          tracer_tnd(i,m) = cflx(i,m)

          ! cummulative tendencies from boundary fluxes
          tracerint%tracer_tnd(i,m) = tracerint%tracer_tnd(i,m) + tracer_tnd(i,m)

          ! expected new values from original values plus boundary fluxes
          tracer_xpd(i,m) = tracerint%tracer(i,m) + tracerint%tracer_tnd(i,m)*ztodt

          ! relative error, expected value - input value / original
          tracer_rer(i,m) = (tracer_xpd(i,m) - tracer_inp(i,m)) / tracerint%tracer(i,m)
       end do

!! final loop for error checking
!    do i = 1, ncol

!! error messages
!       if (abs(enrgy_rer(i)) > 1.E-14 .or. abs(water_rer(i)) > 1.E-14) then
!          tracerint%count = tracerint%count + 1
!          write(iulog,*) "significant conservations error after ", name,        &
!               " count", tracerint%count, " nstep", nstep, "chunk", lchnk, "col", i
!          write(iulog,*) enrgy_inp(i),enrgy_xpd(i),enrgy_dif(i),tracerint%enrgy_tnd(i)*ztodt,  &
!               enrgy_tnd(i)*ztodt,enrgy_rer(i)
!          write(iulog,*) water_inp(i),water_xpd(i),water_dif(i),tracerint%water_tnd(i)*ztodt,  &
!               water_tnd(i)*ztodt,water_rer(i)
!       end if
!    end do


       ! final loop for error checking
       if ( maxval(tracer_rer) > 1.E-14_r8 ) then
          write(iulog,*) "CHECK_TRACERS TRACER large rel error"
          write(iulog,*) tracer_rer
       endif

       do i = 1, ncol
          ! error messages
          if (abs(tracer_rer(i,m)) > 1.E-14_r8 ) then
             tracerint%count = tracerint%count + 1
             write(iulog,*) "CHECK_TRACERS TRACER significant conservation error after ", name,        &
                  " count", tracerint%count, " nstep", nstep, "chunk", lchnk, "col",i
             write(iulog,*)' process name, tracname, index ',  name, tracname, m
             write(iulog,*)" input integral              ",tracer_inp(i,m)
             write(iulog,*)" expected integral           ", tracer_xpd(i,m)
             write(iulog,*)" input - inital integral     ",tracer_dif(i,m)
             write(iulog,*)" cumulative tend      ",tracerint%tracer_tnd(i,m)*ztodt
             write(iulog,*)" process tend         ",tracer_tnd(i,m)*ztodt
             write(iulog,*)" relative error       ",tracer_rer(i,m)
             call endrun()
          end if
       end do
    end do

    return
  end subroutine check_tracers_chng

!#######################################################################

  subroutine tot_energy_phys(state, outfld_name_suffix,vc,enthalpy_flux_incr)
    use physconst,       only: rga,rearth,omega
    use cam_thermo,      only: get_hydrostatic_energy,thermo_budget_num_vars,thermo_budget_vars, &
         wvidx,wlidx,wiidx,seidx,poidx,keidx,moidx,mridx,ttidx,teidx
    use air_composition, only: cp_or_cv_dycore_init
    use cam_history,     only: outfld
    use dyn_tests_utils, only: vc_physics, vc_height, vc_dry_pressure

    use cam_abortutils,  only: endrun
    use cam_history_support, only: max_fieldname_len
    use cam_budget,      only: thermo_budget_history
!------------------------------Arguments--------------------------------

    type(physics_state), intent(inout) :: state
    character(len=*),    intent(in)    :: outfld_name_suffix ! suffix for "outfld"
    integer, optional,   intent(in)    :: vc                 ! vertical coordinate
    real(r8), optional,  intent(in)    :: enthalpy_flux_incr(pcols)
!---------------------------Local storage-------------------------------
    real(r8) :: se(pcols)                          ! Dry Static energy (J/m2)
    real(r8) :: po(pcols)                          ! surface potential or potential energy (J/m2)
    real(r8) :: ke(pcols)                          ! kinetic energy    (J/m2)
    real(r8) :: te(pcols)                          ! total energy incl. latent heat terms  (J/m2)
    real(r8) :: wv(pcols)                          ! column integrated vapor       (kg/m2)
    real(r8) :: liq(pcols)                         ! column integrated liquid      (kg/m2)
    real(r8) :: ice(pcols)                         ! column integrated ice         (kg/m2)
    real(r8) :: tt(pcols)                          ! column integrated test tracer (kg/m2)
    real(r8) :: mr(pcols)                          ! column integrated wind axial angular momentum (kg*m2/s)
    real(r8) :: mo(pcols)                          ! column integrated mass axial angular momentum (kg*m2/s)
    real(r8) :: tt_tmp,mr_tmp,mo_tmp,cos_lat
    real(r8) :: mr_cnst, mo_cnst
    real(r8) :: cp_or_cv(pcols,pver)               ! cp for pressure-based vcoord and cv for height vcoord
    real(r8) :: temp(pcols,pver)                   ! temperature
    real(r8) :: scaling(pcols,pver)                ! scaling for conversion of temperature increment
    real(r8) :: heating
    
    integer :: lchnk                               ! chunk identifier
    integer :: ncol                                ! number of atmospheric columns
    integer :: i,k                                 ! column, level indices
    integer :: vc_loc                              ! local vertical coordinate variable
    integer :: ixtt                                ! test tracer index
    character(len=max_fieldname_len) :: name_out(thermo_budget_num_vars)

!-----------------------------------------------------------------------

    if (.not.thermo_budget_history) return

    do i=1,thermo_budget_num_vars
       name_out(i)=trim(thermo_budget_vars(i))//'_'//trim(outfld_name_suffix)
    end do

    lchnk = state%lchnk
    ncol  = state%ncol

    if (present(vc)) then
      vc_loc = vc
    else
      vc_loc = vc_physics
    end if

    if (state%psetcols == pcols) then
      if (vc_loc == vc_height .or. vc_loc == vc_dry_pressure) then
        cp_or_cv(:ncol,:) = cp_or_cv_dycore(:ncol,:,lchnk)
      else
        cp_or_cv(:ncol,:) = cpairv(:ncol,:,lchnk)
      end if
    else
      call endrun('tot_energy_phys: energy diagnostics not implemented/tested for subcolumns')
    end if

    if (vc_loc == vc_height .or. vc_loc == vc_dry_pressure) then
      scaling(:ncol,:) = cpairv(:ncol,:,lchnk)/cp_or_cv(:ncol,:)!scaling for energy consistency
    else
      scaling(:ncol,:) = 1.0_r8 !internal energy / enthalpy same as CAM physics
    end if
    if (.not.(vc_loc == vc_dry_pressure)) then!xxx
       ! scale accumulated temperature increment for internal energy / enthalpy consistency
       temp(1:ncol,:) = state%temp_ini(1:ncol,:)+scaling(1:ncol,:)*(state%T(1:ncol,:)- state%temp_ini(1:ncol,:))
    else
       do k = 1, pver
          do i = 1, ncol
             heating = cpairv(i,k,lchnk)*state%pdel(i,k)*(state%T(i,k)- state%temp_ini(i,k))!xxx
             temp(i,k) = (heating/state%pdeldry(i,k)+state%temp_ini(i,k)*cp_or_cv_dycore_init(i,k,lchnk))/cp_or_cv_dycore(i,k,lchnk)!xxx
          end do
       end do
    end if
    call get_hydrostatic_energy(state%q(1:ncol,1:pver,1:pcnst),.true.,               &
         state%pdel(1:ncol,1:pver), cp_or_cv(1:ncol,1:pver),                         &
         state%u(1:ncol,1:pver), state%v(1:ncol,1:pver), temp(1:ncol,1:pver),        &
         vc_loc, ptop=state%pintdry(1:ncol,1), phis = state%phis(1:ncol),            &
         te = te(1:ncol), z_mid = state%z_ini(1:ncol,:), se = se(1:ncol),            &
         po = po(1:ncol), ke = ke(1:ncol), wv = wv(1:ncol), liq = liq(1:ncol),       &
         ice = ice(1:ncol))

    if (present(enthalpy_flux_incr)) then
       se(:ncol) = se(:ncol)-enthalpy_flux_incr(:ncol)
    end if
    
    call cnst_get_ind('TT_LW' , ixtt    , abort=.false.)
    tt    = 0._r8
    if (ixtt > 1) then
      if (name_out(ttidx) == 'TT_pAM'.or.name_out(ttidx) == 'TT_zAM') then
        !
        ! after dme_adjust mixing ratios are all wet
        !
        do k = 1, pver
          do i = 1, ncol
            tt_tmp   = state%q(i,k,ixtt)*state%pdel(i,k)*rga
            tt   (i) = tt(i)    + tt_tmp
          end do
        end do
      else
        do k = 1, pver
          do i = 1, ncol
            tt_tmp   = state%q(i,k,ixtt)*state%pdeldry(i,k)*rga
            tt   (i) = tt(i)    + tt_tmp
          end do
        end do
      end if
    end if

    call outfld(name_out(seidx)  ,se      , pcols   ,lchnk   )
    call outfld(name_out(poidx)  ,po      , pcols   ,lchnk   )
    call outfld(name_out(keidx)  ,ke      , pcols   ,lchnk   )
    call outfld(name_out(wvidx)  ,wv      , pcols   ,lchnk   )
    call outfld(name_out(wlidx)  ,liq     , pcols   ,lchnk   )
    call outfld(name_out(wiidx)  ,ice     , pcols   ,lchnk   )
    call outfld(name_out(ttidx)  ,tt      , pcols   ,lchnk   )
    call outfld(name_out(teidx)  ,te      , pcols   ,lchnk   )
    !
    ! Axial angular momentum diagnostics
    !
    ! Code follows
    !
    ! Lauritzen et al., (2014): Held-Suarez simulations with the Community Atmosphere Model
    ! Spectral Element (CAM-SE) dynamical core: A global axial angularmomentum analysis using Eulerian
    ! and floating Lagrangian vertical coordinates. J. Adv. Model. Earth Syst. 6,129-140,
    ! doi:10.1002/2013MS000268
    !
    ! MR is equation (6) without \Delta A and sum over areas (areas are in units of radians**2)
    ! MO is equation (7) without \Delta A and sum over areas (areas are in units of radians**2)
    !
    
    mr_cnst = rga*rearth**3
    mo_cnst = rga*omega*rearth**4
    
    mr = 0.0_r8
    mo = 0.0_r8
    do k = 1, pver
       do i = 1, ncol
          cos_lat = cos(state%lat(i))
          mr_tmp = mr_cnst*state%u(i,k)*state%pdel(i,k)*cos_lat
          mo_tmp = mo_cnst*state%pdel(i,k)*cos_lat**2
          
          mr(i) = mr(i) + mr_tmp
          mo(i) = mo(i) + mo_tmp
       end do
    end do
    
    call outfld(name_out(mridx)  ,mr, pcols,lchnk   )
    call outfld(name_out(moidx)  ,mo, pcols,lchnk   )

  end subroutine tot_energy_phys

  subroutine pressure_enthalpy_adjustment(ncol, lchnk, state, cam_in, cam_out, pbuf, ztodt, itim_old,&
       qini,totliqini,toticeini,tend)
    use cam_thermo,     only: compute_enthalpy_flux, get_hydrostatic_energy
    use air_composition,only: thermodynamic_active_species_liq_num!xxx
    use air_composition,only: thermodynamic_active_species_liq_idx!xxx
    use air_composition,only: thermodynamic_active_species_ice_num!xxx
    use air_composition,only: thermodynamic_active_species_ice_idx!xxx
    use camsrfexch,      only: cam_out_t, cam_in_t, get_prec_vars, get_enthalpy_flux, get_falling_precip_T
    use physics_types,   only: ifrain, ifsnow, ihrain, ihsnow
    use physics_buffer,  only: physics_buffer_desc, pbuf_set_field
    use dyn_tests_utils, only: vc_dycore, vc_height, vc_dry_pressure
    use dycore,          only: dycore_is
    use cam_thermo,      only: cam_thermo_water_update
    use air_composition, only: cpairv, cp_or_cv_dycore, cp_or_cv_dycore_init
    use physics_types,   only: physics_dme_adjust, set_dry_to_wet, dyn_te_idx
    use cam_budget,      only: thermo_budget_history
    use cam_abortutils,  only: endrun
    use cam_history,      only: outfld!xxx
    use physconst,       only: latvap, latice
    integer,             intent(in)    :: ncol, lchnk
    type(physics_state), intent(inout) :: state
    type(cam_in_t),      intent(inout) :: cam_in
    type(cam_out_t),     intent(inout) :: cam_out
    type(physics_buffer_desc), pointer :: pbuf(:)        
    real(r8),            intent(in)    :: ztodt
    integer,             intent(in)    :: itim_old
    
    real(r8), dimension(pcols,pver), intent(in) :: qini, totliqini, toticeini
    type(physics_tend )    , intent(inout) :: tend!xxx    
    real(r8) :: enthalpy_flux_incr(pcols)

    logical  :: moist_mixing_ratio_dycore
    real(r8) :: frain(pcols), fsnow(pcols)
    real(r8) :: te_start(pcols)
    real(r8) :: tot_wv(pcols), tot_ice(pcols), tot_liq(pcols)!xxx
    real(r8) :: heating(pcols)!xxx
    real(r8) :: Tprec(pcols) ! temperature of falling precipitation
    real(r8) :: Q_cpdry      ! accumulated heating for CAM physics
    real(r8) :: imbalance(pcols) ! energy imbalance
    
    integer  :: k,i,m_cnst,m
    tot_wv=0.0_r8!xxx
    do k = 1, pver!xxx
       do i = 1, ncol!xxx
          tot_wv(i)   = tot_wv(i)+(state%q(i,k,1)-qini(i,k))*state%pdel(i,k)*rga!xxx
       end do!xxx
    end do!xxx
    call outfld ('dtot_wv',tot_wv, pcols, lchnk)!xxx
    
    tot_liq=0.0_r8!xxx
    do k = 1, pver!xxx
       do i = 1, ncol!xxx
          do m_cnst=1,thermodynamic_active_species_liq_num
             m = thermodynamic_active_species_liq_idx(m_cnst)
             tot_liq(i)   = tot_liq(i)+(state%q(i,k,m))*state%pdel(i,k)*rga!xxx
          end do
          tot_liq(i)   = tot_liq(i)-totliqini(i,k)*state%pdel(i,k)*rga!xxx
       end do!xxx
    end do!xxx
    call outfld ('dtot_liq',tot_liq, pcols, lchnk)!xxx
    
    tot_ice=0.0_r8!xxx
    do k = 1, pver!xxx
       do i = 1, ncol!xxx
          do m_cnst=1,thermodynamic_active_species_ice_num
             m = thermodynamic_active_species_ice_idx(m_cnst)
             tot_ice(i)   = tot_ice(i)+state%q(i,k,m)*state%pdel(i,k)*rga!xxx
          end do
          tot_ice(i)   = tot_ice(i)-toticeini(i,k)*state%pdel(i,k)*rga!xxx
       end do!xxx
    end do!xxx
    call outfld ('dtot_ice',tot_ice, pcols, lchnk)!xxx

    heating=0.0_r8!xxx
    do k = 1, pver!xxx
       do i = 1, ncol!xxx
          heating(i)   = heating(i)+cpair*(state%T(i,k)-state%temp_ini(i,k))*state%pdel(i,k)/gravit
       end do!xxx
    end do!xxx
    call outfld ('heating',heating, pcols, lchnk)!xxx
    call outfld ('radiation',tend%flx_net, pcols, lchnk)!xxx
    call outfld ('te_sen',tend%te_sen, pcols, lchnk)!xxx
    call outfld ('te_lat',tend%te_lat, pcols, lchnk)!xxx

!    if (compute_enthalpy_flux) then



    
    call pbuf_set_field(pbuf, teout_idx, state%te_cur(:,dyn_te_idx), (/1,itim_old/),(/pcols,1/))!xxx remove

    if (compute_enthalpy_flux) then
       if (.not.dycore_is('SE')) then
          call endrun('FATAL: enthalpy flux is only supported with SE dycore')
       end if
       if (state%psetcols .ne. pcols) then
          call endrun('FATAL: enthalpy flux not supported for sub-columns')
       end if

       call physics_dme_adjust(state, qini, totliqini, toticeini, ztodt)
       call tot_energy_phys(state, 'phAC')
       call tot_energy_phys(state, 'dyAC', vc=vc_dycore)       
       !
       ! get tphysac enthalpy flux
       !
       call get_prec_vars(ncol,pbuf,state%hflx_ac(:,ifrain),state%hflx_ac(:,ifsnow))
       state%hflx_ac(:ncol,ifrain) = state%hflx_ac(:ncol,ifrain) - state%hflx_bc(:ncol,ifrain)
       state%hflx_ac(:ncol,ifsnow) = state%hflx_ac(:ncol,ifsnow) - state%hflx_bc(:ncol,ifsnow)
       call outfld ('FRAIN_tot',state%hflx_ac(:,ifrain), pcols, lchnk)!xxx diags will remove
       call outfld ('FSNOW_tot',state%hflx_ac(:,ifsnow), pcols, lchnk)!xxx diags will remove
       call outfld ('FEVAP_tot',cam_in%cflx(:,1), pcols, lchnk)       !xxx diags will remove
       call outfld ('FRAIN_AC',state%hflx_ac(:,ifrain), pcols, lchnk) !xxx diags will remove
       call outfld ('FSNOW_AC',state%hflx_ac(:,ifsnow), pcols, lchnk) !xxx diags will remove
       call outfld ('FEVAP_AC',cam_in%cflx(:,1), pcols, lchnk)        !xxx diags will remove
       
       Tprec = get_falling_precip_T(ncol,state%T,state%pdel)
       call get_enthalpy_flux(ncol,state%hflx_ac(:,ifrain),state%hflx_ac(:,ifsnow),&
            cam_in%cflx(:,1),Tprec,Tprec,&
            cam_in%ts,state%hflx_ac(:,ihrain),state%hflx_ac(:,ihsnow),cam_out%hevap)
       !
       ! add enthalpy flux to cummulative boundary fluxes
       !
       call outfld ('te_tnd',tend%te_tnd, pcols, lchnk)
       tend%te_tnd(:ncol) = tend%te_tnd(:ncol)+&
                            cam_out%hevap(:ncol)+state%hflx_bc(:ncol,ihsnow)+&
                            state%hflx_ac(:ncol,ihrain)+state%hflx_ac(:ncol,ihsnow)+&
                            state%hflx_bc(:ncol,ihrain)+state%hflx_bc(:ncol,ihsnow)

       call cam_thermo_water_update(state%q(:ncol,:,:), lchnk, ncol, vc_dycore,&
            to_dry_factor=state%pdel(:ncol,:)/state%pdeldry(:ncol,:))
       !
       ! add heating under variable latent heat assumption and compute temperature
       !
       do k = 1, pver
          do i = 1, ncol
             Q_cpdry = cpairv(i,k,lchnk)*state%pdel(i,k)*(state%T(i,k)- state%temp_ini(i,k))
             state%T(i,k) = (Q_cpdry/state%pdeldry(i,k)+state%temp_ini(i,k)*cp_or_cv_dycore_init(i,k,lchnk))/cp_or_cv_dycore(i,k,lchnk)
          end do
       end do
       !
       ! update pdel
       !
       
       !
       !
       !


       !
       ! compute column integrated energy
       !
       call get_hydrostatic_energy(state%q(1:ncol,1:pver,1:pcnst),.true.,              &
           state%pdel(1:ncol,1:pver), cp_or_cv_dycore(1:ncol,1:pver,lchnk),            &
           state%u(1:ncol,1:pver), state%v(1:ncol,1:pver), state%T(1:ncol,1:pver),     &
           vc_dry_pressure, ptop=state%pintdry(1:ncol,1), phis = state%phis(1:ncol),   &
           te = state%te_cur(1:ncol,dyn_te_idx))


       imbalance(:ncol) = (state%te_cur(1:ncol,dyn_te_idx)-state%te_ini(1:ncol,dyn_te_idx))/ztodt-tend%te_tnd(:ncol)

       call outfld ('imbalance',imbalance(:), pcols, lchnk)

       
       call outfld ('HRAIN_AC',state%hflx_ac(:,ihrain), pcols, lchnk)
       call outfld ('HSNOW_AC',state%hflx_ac(:,ihsnow), pcols, lchnk)
       call outfld ('HEVAP_AC',cam_out%hevap, pcols, lchnk)



       
!       enthalpy_flux_incr(:ncol)=ztodt*(&
!            cam_out%hevap(:ncol))
!            state%hflx_bc(:ncol,ihrain)+state%hflx_bc(:ncol,ihsnow))
!                   state%hflx_ac(:ncol,ihrain)+state%hflx_ac(:ncol,ihsnow)+&
!            state%hflx_bc(:ncol,ihrain)+state%hflx_bc(:ncol,ihsnow)+cam_out%hevap(:ncol))


       te_start(:ncol) = state%te_cur(:ncol,dyn_te_idx)!save total energy in each column




       call outfld ('FRHS_FLX',(cam_in%cflx(:,1)*(latvap+latice)-state%hflx_ac(:,ifrain)*latice+cam_in%shf), pcols, lchnk)
!xxx       call outfld ('FRHS_FLX',cam_in%cflx(:,1)*(latvap+latice)+state%hflx_ac(:,ifrain)*latice+cam_in%shf, pcols, lchnk)
       call outfld ('FRHS_FLXA',cam_in%cflx(:,1)*(latvap+latice), pcols, lchnk)
       call outfld ('FRHS_FLXB',state%hflx_ac(:,ifrain)*latice, pcols, lchnk)
       call outfld ('FRHS_FLXC',cam_in%shf, pcols, lchnk)


       


       
!       call tot_energy_phys(state, 'phAC',enthalpy_flux_incr=enthalpy_flux_incr)
!       call tot_energy_phys(state, 'dyAC',enthalpy_flux_incr=enthalpy_flux_incr, vc=vc_dycore)
    else
       call tot_energy_phys(state, 'phAC')
       call tot_energy_phys(state, 'dyAC', vc=vc_dycore)
       call pbuf_set_field(pbuf, teout_idx, state%te_cur(:,dyn_te_idx), (/1,itim_old/),(/pcols,1/))
    end if
    !
    ! update cp/cv for energy computation based in updated water variables
    !
!    call cam_thermo_water_update(state%q(:ncol,:,:), lchnk, ncol, vc_dycore,&
!         to_dry_factor=state%pdel(:ncol,:)/state%pdeldry(:ncol,:))
!
!    call physics_dme_adjust(state, qini, totliqini, toticeini, ztodt)
    
    moist_mixing_ratio_dycore = dycore_is('LR').or. dycore_is('FV3')
    if (moist_mixing_ratio_dycore) then
       !
       ! change dry mixing ratios for moist for moist mixing ratio dycores
       !
       ! Note: this operation will NOT be reverted with set_wet_to_dry after set_dry_to_wet call
       !
       call set_dry_to_wet(state)
    endif
    call tot_energy_phys(state, 'phAF',enthalpy_flux_incr=enthalpy_flux_incr)!xxx , ztodt)!xxx should be AF
    call tot_energy_phys(state, 'dyAF', vc=vc_dycore,enthalpy_flux_incr=enthalpy_flux_incr)!xxx ztodt, vc=vc_dycore)!xxx should be AF
!    if (compute_enthalpy_flux) then
!       call column_fixer_using_variable_latent_heat_energy_formula(ncol, lchnk, state, te_start, enthalpy_flux_incr)

!       call pbuf_set_field(pbuf, teout_idx, state%te_cur(:,dyn_te_idx), (/1,itim_old/),(/pcols,1/))
!    end if
    call tot_energy_phys(state, 'phAM', enthalpy_flux_incr=enthalpy_flux_incr)
    call tot_energy_phys(state, 'dyAM', vc=vc_dycore,enthalpy_flux_incr=enthalpy_flux_incr)
  end subroutine pressure_enthalpy_adjustment

  subroutine column_fixer_using_variable_latent_heat_energy_formula(ncol, lchnk, state, enthalpy_flux_incr,te_start)
    use air_composition, only: cp_or_cv_dycore
    use dyn_tests_utils, only: vc_dry_pressure
    use cam_thermo,      only: get_hydrostatic_energy
    use physics_types,   only: dyn_te_idx
    use physconst,       only: rga
    integer,             intent(in)    :: ncol, lchnk
    type(physics_state), intent(inout) :: state
    real(r8),            intent(in)    :: enthalpy_flux_incr(pcols)
    real(r8),            intent(in)    :: te_start(pcols)

    real(r8) :: temp(state%ncol,pver) ! temperature
    real(r8) :: scaling(pcols,pver)   ! scaling for conversion of temperature increment
    real(r8) :: te_discr(pcols)       ! total energy discrepancy in each column
    real(r8) :: column_enthalpy(pcols)! column integrated enthalpy
    real(r8) :: cpave(pcols)
    real(r8) :: Tfix
    integer  :: i,k
    !
    ! compute current dycore energy (updated by dme_adjust)
    !
    !
    ! enthalpy scaling for energy consistency
    !
    scaling(:ncol,:) = cpairv(:ncol,:,lchnk)/cp_or_cv_dycore(:ncol,:,lchnk)
    temp(1:ncol,:)   = state%temp_ini(1:ncol,:)+scaling(1:ncol,:)*(state%T(1:ncol,:)-state%temp_ini(1:ncol,:))
    call get_hydrostatic_energy(state%q(1:ncol,1:pver,1:pcnst),.true.,               &
         state%pdel(1:ncol,1:pver), cp_or_cv_dycore(1:ncol,1:pver,lchnk),            &
         state%u(1:ncol,1:pver), state%v(1:ncol,1:pver), temp(1:ncol,1:pver),        &
         vc_dry_pressure, ptop=state%pintdry(1:ncol,1), phis = state%phis(1:ncol),   &
         te = state%te_cur(1:ncol,dyn_te_idx), se = column_enthalpy(1:ncol),         &
         H2O = state%tw_cur(1:ncol,dyn_te_idx))
    !
    ! total energy discrepancy
    !
    te_discr(:ncol) = state%te_cur(:ncol,dyn_te_idx)-te_start(:ncol)-enthalpy_flux_incr(:ncol)
    !
    ! compute temperature increment (fixer) so that te_discr(:ncol)=0
    !
    cpave = 0._r8
    do k=1,pver
       do i=1,ncol
          cpave(i) = cpave(i) + cp_or_cv_dycore(i,k,lchnk) * state%pdel(i,k) * rga
       end do
    end do
    do i=1,ncol
       Tfix = te_discr(i)/cpave(i)
       do k=1,pver
          state%T(i,k) = (state%T(i,k)-Tfix)/scaling(i,k)
       end do
    end do
  end subroutine column_fixer_using_variable_latent_heat_energy_formula
end module check_energy
