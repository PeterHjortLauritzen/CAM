module enthalpy_flux_mod
  use shr_kind_mod,     only: r8 => shr_kind_r8
  use cam_abortutils,   only: endrun
  implicit none
  private
  
  public enthalpy_evap_tend
  public enthalpy_clubb_pumas_tend
  public enthalpy_zm_tend
  public get_enthalpy_flux, get_prec_vars

contains
  !
  ! Compute enthalpy flux associated with evaporation
  ! section 2.1.4 in https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2022MS003117
  !
  subroutine enthalpy_evap_tend(pbuf,state,ptend,evap_enthalpy_flux)
    use physics_types,      only: physics_ptend, physics_ptend_init, physics_state
    use physics_buffer,     only: pbuf_get_index, pbuf_get_field, physics_buffer_desc,physics_buffer_desc
    use ppgrid,             only: pcols,pver
    use physconst,          only: gravit
    type(physics_buffer_desc), pointer :: pbuf(:)
    type(physics_state), intent(in)    :: state
    type(physics_ptend), intent(out)   :: ptend
    real(r8),            intent(out)   :: evap_enthalpy_flux(pcols)

    integer  :: ncol, i, lchnk
    integer  :: hevap_iceref_idx
    real(r8), pointer :: hevap_iceref(:)
    
    ncol  = state%ncol
    lchnk = state%lchnk
    call physics_ptend_init(ptend,state%psetcols, "enthalpy_flux_evap",ls=.true.)
    hevap_iceref_idx = pbuf_get_index('hevap_iceref')
    if (hevap_iceref_idx>0) then
       call pbuf_get_field(pbuf, hevap_iceref_idx, hevap_iceref)
       ptend%s(:ncol,pver) = gravit*hevap_iceref(:ncol)/state%pdel(:ncol,pver)
       evap_enthalpy_flux = hevap_iceref
    else
       evap_enthalpy_flux = 0.0_r8
    end if
  end subroutine enthalpy_evap_tend
  !
  ! Compute enthalpy flux associated with falling precipitation from ZM
  ! section 2.1.4 in https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2022MS003117
  !
  subroutine enthalpy_zm_tend(pbuf,state,ptend,prect_enthalpy_flux)
    use physics_types,      only: physics_ptend, physics_ptend_init, physics_state
    use physics_buffer,     only: pbuf_get_index, pbuf_get_field, physics_buffer_desc,physics_buffer_desc
    use ppgrid,             only: pcols,pver
    use physconst,          only: gravit
    use cam_history,        only: outfld
    use physics_types,      only: ifrain, ifsnow, ihrain, ihsnow
    use physconst,          only: cpliq, cpice
    use air_composition,    only: cpairv
    type(physics_buffer_desc), pointer :: pbuf(:)
    type(physics_state), intent(inout) :: state
    type(physics_ptend), intent(out)   :: ptend
    real(r8),            intent(out)   :: prect_enthalpy_flux(pcols)

    integer  :: ncol, i, lchnk, k
    real(r8) :: ttend
    ! Temperature change due to deep convection.
    real(r8), pointer :: ttend_dp(:,:)
    integer           :: ttend_dp_idx = -1
    real(r8) :: flx(pcols,pver)  !xxx debug
    integer  :: kcloud(pcols)
    !
    ! initialization
    !
    ncol  = state%ncol
    lchnk = state%lchnk
    call physics_ptend_init(ptend,state%psetcols, "enthalpy_flux_zm",ls=.true.)
    !
    ! get falling precipitation fluxes
    !
    call get_prec_vars(ncol,pbuf,state%hflx_bc(:,ifrain),state%hflx_bc(:,ifsnow))
    !
    ! precipitation variables have been assumulated since last call to coupler hence
    ! subtract precipitation from tphys_ac from last time-step
    !
    state%hflx_bc(:ncol,ifrain) = state%hflx_bc(:,ifrain) - state%hflx_ac(:ncol,ifrain)
    state%hflx_bc(:ncol,ifsnow) = state%hflx_bc(:,ifsnow) - state%hflx_ac(:ncol,ifsnow)
    call outfld ('FRAIN_BC',state%hflx_bc(:,ifrain), pcols, lchnk) !xxx diags will remove
    call outfld ('FSNOW_BC',state%hflx_bc(:,ifsnow), pcols, lchnk) !xxx diags will remove
    !
    ! compute enthalpy fluxes from falling precipitation
    !
    do i=1,ncol
       state%hflx_bc(:,ihrain) = -state%hflx_bc(:ncol,ifrain)*state%T(i,pver)*cpliq !sign following atmosphere model convention: -ve for flux out of atmosphere
       state%hflx_bc(:,ihsnow) = -state%hflx_bc(:ncol,ifsnow)*state%T(i,pver)*cpice
    end do
    call outfld('HRAIN_BC', state%hflx_bc(:,ihrain),  pcols, lchnk)!xxx diags will remove
    call outfld('HSNOW_BC', state%hflx_bc(:,ihsnow),  pcols, lchnk)!xxx diags will remove

    
    kcloud = 30
    do i=1, ncol
       prect_enthalpy_flux(i) = state%hflx_bc(i,ihrain)+state%hflx_bc(i,ihsnow)
       !
       ! uniform tendency in each layer from surface to k=kcloud
       !
       ttend               = prect_enthalpy_flux(i)*gravit/SUM(state%pdel(i,kcloud(i):pver)*cpairv(i,kcloud(i):pver,lchnk))
       do k=kcloud(i),pver
          ptend%s(i,k)     = ttend*cpairv(i,k,lchnk)
       end do
    end do
    
    !
    ! get heating from ZM to ditribute heating in column
    !
    ! ttend_dp(:state%ncol,:pver) = ptend%s(:state%ncol,:pver)/cpair
    !
!    ttend_dp_idx = pbuf_get_index('TTEND_DP')
!    if (ttend_dp_idx<1) then
!       call endrun('enthalpy_zm_tend: TTEND_DP does not exist')
!    end if
!    call pbuf_get_field(pbuf, ttend_dp_idx, ttend_dp)

  end subroutine enthalpy_zm_tend

    !
  ! Compute enthalpy flux associated with falling precipitation from ZM
  ! section 2.1.4 in https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2022MS003117
  !
  subroutine enthalpy_clubb_pumas_tend(pbuf,state,ptend,prect_enthalpy_flux)
    use physics_types,      only: physics_ptend, physics_ptend_init, physics_state
    use physics_buffer,     only: pbuf_get_index, pbuf_get_field, physics_buffer_desc,physics_buffer_desc
    use ppgrid,             only: pcols,pver
    use physconst,          only: gravit
    use cam_history,        only: outfld
    use physics_types,      only: ifrain, ifsnow, ihrain, ihsnow
    use physconst,          only: cpliq, cpice
    use air_composition,    only: cpairv
    type(physics_buffer_desc), pointer :: pbuf(:)
    type(physics_state), intent(inout) :: state
    type(physics_ptend), intent(out)   :: ptend
    real(r8),            intent(out)   :: prect_enthalpy_flux(pcols)

    integer  :: ncol, i, lchnk, k
    real(r8) :: ttend
    integer  :: kcloud(pcols)
    !
    ! initialization
    !
    ncol  = state%ncol
    lchnk = state%lchnk
    call physics_ptend_init(ptend,state%psetcols, "enthalpy_flux_clubb_pumas",ls=.true.)
    !
    ! get falling precipitation fluxes
    !
    call get_prec_vars(ncol,pbuf,state%hflx_ac(:,ifrain),state%hflx_ac(:,ifsnow))
    !
    ! precipitation variables have been assumulated since last call to coupler hence
    ! subtract precipitation from tphys_ac from last time-step
    !
    state%hflx_ac(:ncol,ifrain) = state%hflx_ac(:,ifrain) - state%hflx_bc(:ncol,ifrain)
    state%hflx_ac(:ncol,ifsnow) = state%hflx_ac(:,ifsnow) - state%hflx_bc(:ncol,ifsnow)
    call outfld ('FRAIN_AC',state%hflx_ac(:,ifrain), pcols, lchnk) !xxx diags will remove
    call outfld ('FSNOW_AC',state%hflx_ac(:,ifsnow), pcols, lchnk) !xxx diags will remove
    !
    ! compute enthalpy fluxes from falling precipitation
    !
    do i=1,ncol
       state%hflx_ac(:,ihrain) = -state%hflx_ac(:ncol,ifrain)*state%T(i,pver)*cpliq !sign following atmosphere model convention: -ve for flux out of atmosphere
       state%hflx_ac(:,ihsnow) = -state%hflx_ac(:ncol,ifsnow)*state%T(i,pver)*cpice
    end do

    kcloud = 30
    do i=1, ncol
       prect_enthalpy_flux(i) = state%hflx_ac(i,ihrain)+state%hflx_ac(i,ihsnow)
       !
       ! uniform tendency in each layer from surface to k=kcloud
       !
       ttend               = prect_enthalpy_flux(i)*gravit/SUM(state%pdel(i,kcloud(i):pver)*cpairv(i,kcloud(i):pver,lchnk))
       do k=kcloud(i),pver
          ptend%s(i,k)     = ttend*cpairv(i,k,lchnk)
       end do
    end do
  end subroutine enthalpy_clubb_pumas_tend

   subroutine get_enthalpy_flux(ncol,frain,fsnow,fevap,Train,Tsnow,Tevap,hrain_iceref,hsnow_iceref,hevap_iceref)
     use physconst,        only: cpwv, cpliq, cpice
     use ppgrid,           only: pcols,pver
     integer,                    intent(in) :: ncol
     real(r8), dimension(pcols), intent(in) :: frain
     real(r8), dimension(pcols), intent(in) :: fsnow
     real(r8), dimension(pcols), intent(in) :: fevap
     real(r8), dimension(pcols), intent(in) :: Train
     real(r8), dimension(pcols), intent(in) :: Tsnow
     real(r8), dimension(pcols), intent(in) :: Tevap
     !
     !enthalpy flux using ice reference (atmosphere)
     !
     real(r8), dimension(pcols), intent(out) :: hrain_iceref,hsnow_iceref,hevap_iceref
     integer  :: i

     do i=1,ncol
!xxx        hsnow_iceref(i) = -fsnow(i)*280.0*cpice !sign following atmosphere model convention: -ve for flux out of atmopshere
!xxx        hrain_iceref(i) = -frain(i)*280.0*cpliq !sign following atmosphere model convention: -ve for flux out of atmosphere
!xxx        hevap_iceref(i) =  fevap(i)*280.0*cpwv  !sign following atmosphere model convention: +ve for flux into atmosphere
        hsnow_iceref(i) = -fsnow(i)*Tsnow(i)*cpice !sign following atmosphere model convention: -ve for flux out of atmopshere
        hrain_iceref(i) = -frain(i)*Train(i)*cpliq !sign following atmosphere model convention: -ve for flux out of atmosphere
        hevap_iceref(i) =  fevap(i)*Tevap(i)*cpwv  !sign following atmosphere model convention: +ve for flux into atmosphere
     end do
   end subroutine get_enthalpy_flux

subroutine get_prec_vars(ncol,pbuf,frain,fsnow,&
     precc_out,precl_out,precsc_out,precsl_out,rliqbc_out,rice_out)
     use ppgrid, only: pcols
     use physics_buffer,   only: pbuf_get_index, pbuf_get_field, physics_buffer_desc

     integer, intent(in) :: ncol
     type(physics_buffer_desc), pointer         :: pbuf(:)
     real(r8), dimension(ncol), intent(out):: frain!snow flux
     real(r8), dimension(ncol), intent(out):: fsnow!rain flux

     real(r8), dimension(pcols), optional, intent(out):: precc_out !total precipitation from convection
     real(r8), dimension(pcols), optional, intent(out):: precl_out !total large scale precipitation
     real(r8), dimension(pcols), optional, intent(out):: precsc_out!frozen precipitation from convection
     real(r8), dimension(pcols), optional, intent(out):: precsl_out!frozen large scale precipitation
     real(r8), dimension(pcols), optional, intent(out):: rliqbc_out!reserved liquid
     real(r8), dimension(pcols), optional, intent(out):: rice_out  !reserved ice

     integer :: i

     real(r8), pointer :: prec_dp(:)                 !total precipitation from from deep convection
     real(r8), pointer :: snow_dp(:)                 !frozen precipitation from deep convection
     real(r8), pointer :: prec_sh(:)                 !total precipitation from shallow convection
     real(r8), pointer :: snow_sh(:)                 !frozen precipitation from from shallow convection
     real(r8), pointer :: prec_sed(:)                !total precipitation from cloud sedimentation
     real(r8), pointer :: snow_sed(:)                !frozen precipitation from sedimentation
     real(r8), pointer :: prec_pcw(:)                !total precipitation from from microphysics
     real(r8), pointer :: snow_pcw(:)                !frozen precipitation from from microphysics
     real(r8), pointer :: rliqbc(:)                  !reserved liquid
     real(r8), pointer :: rice(:)                    !reserved ice

     real(r8), dimension(pcols):: precc, precl, precsc, precsl
     integer :: prec_dp_idx, snow_dp_idx, prec_sh_idx, snow_sh_idx
     integer :: prec_sed_idx,snow_sed_idx,prec_pcw_idx,snow_pcw_idx
     integer :: rliqbc_idx,rice_idx
     !
     ! get fields from pbuf
     !
     prec_dp_idx = pbuf_get_index('PREC_DP', errcode=i)
     snow_dp_idx = pbuf_get_index('SNOW_DP', errcode=i)
     prec_sh_idx = pbuf_get_index('PREC_SH', errcode=i)
     snow_sh_idx = pbuf_get_index('SNOW_SH', errcode=i)
     prec_sed_idx = pbuf_get_index('PREC_SED', errcode=i)
     snow_sed_idx = pbuf_get_index('SNOW_SED', errcode=i)
     prec_pcw_idx = pbuf_get_index('PREC_PCW', errcode=i)
     snow_pcw_idx = pbuf_get_index('SNOW_PCW', errcode=i)
     rliqbc_idx   = pbuf_get_index('RLIQBC', errcode=i)

     if (prec_dp_idx > 0) then
        call pbuf_get_field(pbuf, prec_dp_idx, prec_dp)
     end if
     if (snow_dp_idx > 0) then
        call pbuf_get_field(pbuf, snow_dp_idx, snow_dp)
     end if
     if (prec_sh_idx > 0) then
        call pbuf_get_field(pbuf, prec_sh_idx, prec_sh)
     end if
     if (snow_sh_idx > 0) then
        call pbuf_get_field(pbuf, snow_sh_idx, snow_sh)
     end if
     if (prec_sed_idx > 0) then
        call pbuf_get_field(pbuf, prec_sed_idx, prec_sed)
     end if
     if (snow_sed_idx > 0) then
        call pbuf_get_field(pbuf, snow_sed_idx, snow_sed)
     end if
     if (prec_pcw_idx > 0) then
        call pbuf_get_field(pbuf, prec_pcw_idx, prec_pcw)
     end if
     if (snow_pcw_idx > 0) then
        call pbuf_get_field(pbuf, snow_pcw_idx, snow_pcw)
     end if
     if (rliqbc_idx>0) then
        call pbuf_get_field(pbuf, rliqbc_idx, rliqbc)
     end if
     
     precc  = 0._r8
     precl  = 0._r8
     precsc = 0._r8
     precsl = 0._r8
     if (prec_dp_idx > 0) then
        precc(:ncol) = precc(:ncol) + prec_dp(:ncol)
     end if
     if (prec_sh_idx > 0) then
        precc(:ncol)  = precc(:ncol)  + prec_sh(:ncol)
     end if
     if (prec_sed_idx > 0) then
        precl(:ncol) = precl(1:ncol) + prec_sed(:ncol)
     end if
     if (prec_pcw_idx > 0) then
        precl(:ncol)  = precl(1:ncol) + prec_pcw(:ncol)
     end if
     if (snow_dp_idx > 0) then
        precsc(:ncol) = precsc(:ncol) + snow_dp(:ncol)
     end if
     if (snow_sh_idx > 0) then
        precsc(:ncol) = precsc(:ncol) + snow_sh(:ncol)
     end if
     if (snow_sed_idx > 0) then
        precsl(:ncol) = precsl(:ncol) + snow_sed(:ncol)
     end if
     if (snow_pcw_idx > 0) then
        precsl(:ncol)= precsl(:ncol) + snow_pcw(:ncol)
     end if

     do i=1,ncol
        precc(i)  = MAX(precc(i), 0.0_r8)
        precl(i)  = MAX(precl(i), 0.0_r8)
        precsc(i) = MAX(precsc(i),0.0_r8)
        precsl(i) = MAX(precsl(i),0.0_r8)
        if (precsc(i).gt.precc(i)) precsc(i)=precc(i)
        if (precsl(i).gt.precl(i)) precsl(i)=precl(i)
     end do
     if (present(precc_out )) precc_out (:ncol) = precc (:ncol)
     if (present(precl_out )) precl_out (:ncol) = precl (:ncol)
     if (present(precsc_out)) precsc_out(:ncol) = precsc(:ncol)
     if (present(precsl_out)) precsl_out(:ncol) = precsl(:ncol)
     if (present(rliqbc_out)) rliqbc_out(:ncol) = rliqbc(:ncol)
     
     
     fsnow(:) = 1000.0_r8*(precsc(:ncol)+precsl(:ncol))                           !snow flux
     frain(:) = 1000.0_r8*(precc (:ncol)-precsc(:ncol)+precl(:ncol)-precsl(:ncol))!rain flux
   end subroutine get_prec_vars
   
end module enthalpy_flux_mod
