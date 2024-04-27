module held_suarez_cam

  !-----------------------------------------------------------------------
  !
  ! Purpose: Implement idealized Held-Suarez forcings
  !    Held, I. M., and M. J. Suarez, 1994: 'A proposal for the
  !    intercomparison of the dynamical cores of atmospheric general
  !    circulation models.'
  !    Bulletin of the Amer. Meteor. Soc., vol. 75, pp. 1825-1830.
  !
  !-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols, pver

  implicit none
  private
  save

  public :: held_suarez_init, held_suarez_tend

  real(r8), parameter :: efoldf  =  1._r8  ! efolding time for wind dissipation
  real(r8), parameter :: efolda  = 40._r8  ! efolding time for T dissipation
  real(r8), parameter :: efolds  =  4._r8  ! efolding time for T dissipation
  real(r8), parameter :: sigmab  =  0.7_r8 ! threshold sigma level
  real(r8), parameter :: t00     = 200._r8 ! minimum reference temperature
  real(r8), parameter :: kf      = 1._r8/(86400._r8*efoldf) ! 1./efolding_time for wind dissipation

  real(r8), parameter :: onemsig = 1._r8 - sigmab ! 1. - sigma_reference

  real(r8), parameter :: ka      = 1._r8/(86400._r8 * efolda) ! 1./efolding_time for temperature diss.
  real(r8), parameter :: ks      = 1._r8/(86400._r8 * efolds)

!=======================================================================
contains
!=======================================================================

  subroutine held_suarez_init()
    use cam_history,        only: addfld, add_default
    use ref_pres,           only: psurf_ref
#ifdef planet_mars
    use held_suarez_1994_mars,   only: held_suarez_1994_mars_init
#else
    use held_suarez_1994,   only: held_suarez_1994_init
#endif

    character(len=512)            :: errmsg
    integer                       :: errflg

    ! Set model constant values
#ifdef planet_mars
    call held_suarez_1994_mars_init(psurf_ref, errmsg, errflg)
#else
    call held_suarez_1994_init(psurf_ref, errmsg, errflg)
#endif
    ! This field is added by radiation when full physics is used
    call addfld('QRS', (/ 'lev' /), 'A', 'K/s', &
         'Temperature tendency associated with the relaxation toward the equilibrium temperature profile')
    call add_default('QRS', 1, ' ')
 end subroutine held_suarez_init

  subroutine held_suarez_tend(state, ptend, ztodt)
    use air_composition,    only: cappav, cpairv
    use ref_pres,           only: pref_mid_norm
    use phys_grid,          only: get_rlat_all_p
    use physics_types,      only: physics_state, physics_ptend
    use physics_types,      only: physics_ptend_init
    use cam_abortutils,     only: endrun
    use cam_history,        only: outfld
#ifdef planet_mars
    use held_suarez_1994_mars,   only: held_suarez_1994_mars_run
    use std_atm_profile,    only: std_atm_height, std_atm_temp
#else
    use held_suarez_1994,   only: held_suarez_1994_run
#endif
    !
    ! Input arguments
    !
    type(physics_state), intent(inout) :: state
    real(r8),            intent(in)    :: ztodt            ! Two times model timestep (2 delta-t)
                                                           !
                                                           ! Output argument
                                                           !
    type(physics_ptend), intent(out)   :: ptend            ! Package tendencies
                                                           !
    !---------------------------Local workspace-----------------------------

    integer                            :: lchnk            ! chunk identifier
    integer                            :: ncol             ! number of atmospheric columns

    real(r8)                           :: clat(pcols)      ! latitudes(radians) for columns
    real(r8)                           :: pmid(pcols,pver) ! mid-point pressure
#ifdef planet_mars
    real(r8)                           :: zmid(pcols,pver) ! mid-point pseudo-height
    real(r8)                           :: tref(pcols,pver) ! reference temp
#endif
    integer                            :: i, k             ! Longitude, level indices

    character(len=64)                  :: scheme_name      ! CCPP-required variables (not used in CAM)
    character(len=512)            :: errmsg
    integer                       :: errflg

    !
    !-----------------------------------------------------------------------
    !

    lchnk = state%lchnk
    ncol  = state%ncol

    call get_rlat_all_p(lchnk, ncol, clat)
    do k = 1, pver
      do i = 1, ncol
        pmid(i,k) = state%pmid(i,k)
      end do
    end do
#ifdef planet_mars
    do i = 1, ncol
      call std_atm_height(state%pmid(i,:), zmid(i,:))
      call std_atm_temp  (zmid(i,:), clat(i), tref(i,:))
    end do
#endif

    ! initialize individual parameterization tendencies
    call physics_ptend_init(ptend, state%psetcols, 'held_suarez', ls=.true., lu=.true., lv=.true.)
#ifdef planet_mars
    call held_suarez_1994_mars_run(pver, ncol, pref_mid_norm, clat, &
         cpairv(1:ncol,:,lchnk), state%pmid(1:ncol,:),tref(1:ncol,:), &
         state%u(1:ncol,:), state%v(1:ncol,:), state%t(1:ncol,:), &
         ptend%u(1:ncol,:), ptend%v(1:ncol,:), ptend%s(1:ncol,:), &
         scheme_name, errmsg, errflg)
#else
    call held_suarez_1994_run(pver, ncol, pref_mid_norm, clat, cappav(1:ncol,:,lchnk), &
                              cpairv(1:ncol,:,lchnk), state%pmid(1:ncol,:),            &
                              state%u(1:ncol,:), state%v(1:ncol,:), state%t(1:ncol,:), &
                              ptend%u(1:ncol,:), ptend%v(1:ncol,:), ptend%s(1:ncol,:), &
                              scheme_name, errmsg, errflg)
#endif
    ! Note, we assume that there are no subcolumns in simple physics
    pmid(:ncol,:) = ptend%s(:ncol, :)/cpairv(:ncol,:,lchnk)
    if (pcols > ncol) then
      pmid(ncol+1:,:) = 0.0_r8
    end if
    call outfld('QRS', pmid, pcols, lchnk)

  end subroutine held_suarez_tend

end module held_suarez_cam
