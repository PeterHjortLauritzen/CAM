module radiation

!---------------------------------------------------------------------------------
!
! CAM interface to EXORT radiation parameterization
!
!---------------------------------------------------------------------------------

use cam_abortutils,      only: endrun
use cam_history,         only: addfld, add_default, horiz_only, outfld
use cam_logfile,         only: iulog
use camsrfexch,          only: cam_out_t, cam_in_t
use constituents,        only: pcnst
use pio,                 only: file_desc_t, var_desc_t,               &
                               pio_int, pio_double, pio_noerr,        &
                               pio_seterrorhandling, pio_bcast_error, &
                               pio_inq_varid, pio_def_var,            &
                               pio_put_var, pio_get_var, pio_put_att, &
                               pio_closefile, pio_nowrite
use physconst,           only: cappa, cpair,pstd
use phys_grid,           only: get_rlat_all_p,get_rlon_all_p
use physics_buffer,      only: pbuf_get_index, pbuf_set_field, physics_buffer_desc, &
                               pbuf_add_field, dtype_r8, pbuf_old_tim_idx, pbuf_get_field
use physics_types,       only: physics_state, physics_ptend
use physics_types,       only: physics_ptend_init
use ppgrid,              only: pcols, pver, pverp, begchunk, endchunk
use rad_constituents,    only: N_DIAG, rad_cnst_get_call_list
use radgrid,             only: ntot_wavlnrng,camtop,k_h2o, k_co2, k_ch4, k_c2h6, k_o3, k_o2, &
                               kh2h2,kn2h2,kn2n2,ko2o2,ko2n2, ko2co2,kco2co2_lw,kco2co2_sw,kco2h2,kco2ch4, &
                               kh2oself_mtckd, kh2ofrgn_mtckd, ngauss_8gpt, kc_npress, kc_ntemp, kmtckd_ntemp, &
                               kn2n2_ntemp, kn2h2_ntemp, kh2h2_ntemp, kco2co2_sw_ntemp, kco2co2_lw_ntemp, kco2h2_ntemp, &
                               kco2ch4_ntemp, ko2o2_ntemp, ko2n2_ntemp, ko2co2_ntemp, s0, solarflux, qcldliq, wcldliq, &
                               gcldliq, qcldice, wcldice, gcldice, nrei, nrel
use scamMod,             only: scm_crm_mode
use shr_kind_mod,        only: r8=>shr_kind_r8, cl=>shr_kind_cl
use spmd_utils,          only: masterproc,mpicom, mstrid=>masterprocid, mpi_integer, &
                               mpi_logical, mpi_real8, mpi_character
use time_manager,        only: get_nstep, is_first_restart_step, &
                               get_curr_calday, get_step_size, get_curr_date
use exoplanet_mod,       only: do_exo_rt_clearsky, exo_rad_step, do_exo_rt_spectral, exo_diurnal, &
                               exo_n2mmr, exo_h2mmr,exo_c2h6mmr,exo_ndays,exo_ch4mmr
implicit none
private
save

public :: &
   radiation_readnl,         &! read namelist variables
   radiation_register,       &! registers radiation physics buffer fields
   rad_is_active,            &! active radiation pkg
   radiation_define_restart, &! define variables for restart
   radiation_write_restart,  &! write variables to restart
   radiation_read_restart,   &! read variables from restart
   radiation_init,           &! initialize radiation
   radiation_tend,           &! compute heating rates and fluxes
   rad_out_t,                &! type for diagnostic outputs
   radiation_do

   ! Control variables set via namelist
    logical, public, protected :: do_exo_rt_spectral
    logical, public, protected :: do_exo_rt_clearsky
    logical, public, protected :: do_exo_rt_optimize_bands
    logical, public, protected :: do_exo_rt
    logical, public, protected :: do_exo_atmconst
    logical, public, protected :: do_exo_synchronous
    logical, public, protected :: do_exo_gw
    logical, public, protected :: graupel_in_rad

    integer, public, protected :: iradsw = -1     ! freq. of shortwave radiation calc in time steps (positive)
                               ! or hours (negative).
    integer, public, protected :: iradlw = -1     ! frequency of longwave rad. calc. in time steps (positive)
                               ! or hours (negative).
    integer, public, protected :: irad_always = 0 ! Specifies length of time in timesteps (positive)
                               ! or hours (negative) SW/LW radiation will be
                               ! run continuously from the start of an
                               ! initial or restart run
    integer, public, protected :: exo_rad_step = 3 ! freq. of radiation calc in time steps (positive)
                            ! or hours (negative).

    character(len=cl), public, protected :: k_h2o_file, &
    k_co2_file,          &
    k_o2_file,           &
    k_o3_file,           &
    k_ch4_file,          &
    k_c2h6_file,         &
    kh2o_mtckd_file,     &
    kn2n2cia_file,       &
    kn2h2cia_file,       &
    kh2h2cia_file,       &
    kco2co2cia_lw_file,  &
    kco2co2cia_sw_file,  &
    kco2h2cia_file,      &
    kco2ch4cia_file,     &
    ko2o2cia_file,       &
    ko2n2cia_file,       &
    ko2co2cia_file,      &
    cldoptsL_file,       &
    cldoptsI_file,       &
    solar_file

real(r8), public, protected :: nextsw_cday = -1._r8 ! future radiation calday for surface models
   ! Physics buffer indices
   integer :: qrs_idx      = 0
   integer :: qrl_idx      = 0
   integer :: rel_idx      = 0
   integer :: rei_idx      = 0
   integer :: su_idx       = 0
   integer :: sd_idx       = 0
   integer :: lu_idx       = 0
   integer :: ld_idx       = 0
   integer :: fsds_idx     = 0
   integer :: fsns_idx     = 0
   integer :: fsnt_idx     = 0
   integer :: flns_idx     = 0
   integer :: flnt_idx     = 0
   integer :: cicewp_idx   = 0
   integer :: cliqwp_idx   = 0
   integer :: cld_idx      = 0
   integer :: cldfsnow_idx = 0
   integer :: co2_idx      = 0
   integer :: n2_idx      = 0

   character(len=4) :: diag(0:N_DIAG) =(/'    ','_d1 ','_d2 ','_d3 ','_d4 ','_d5 ','_d6 ','_d7 ','_d8 ','_d9 ','_d10'/)

type rad_out_t
   real(r8) :: solin(pcols)         ! Solar incident flux

   real(r8) :: qrsc(pcols,pver)

   real(r8) :: fsntc(pcols)         ! Clear sky total column abs solar flux
   real(r8) :: fsntoa(pcols)        ! Net solar flux at TOA
   real(r8) :: fsntoac(pcols)       ! Clear sky net solar flux at TOA
   real(r8) :: fsutoa(pcols)        ! upwelling solar flux at TOA

   real(r8) :: fsnirt(pcols)        ! Near-IR flux absorbed at toa
   real(r8) :: fsnrtc(pcols)        ! Clear sky near-IR flux absorbed at toa
   real(r8) :: fsnirtsq(pcols)      ! Near-IR flux absorbed at toa >= 0.7 microns

   real(r8) :: fsn200(pcols)        ! fns interpolated to 200 mb
   real(r8) :: fsn200c(pcols)       ! fcns interpolated to 200 mb
   real(r8) :: fsnr(pcols)          ! fns interpolated to tropopause

   real(r8) :: fsnsc(pcols)         ! Clear sky surface abs solar flux
   real(r8) :: fsdsc(pcols)         ! Clear sky surface downwelling solar flux

   real(r8) :: qrlc(pcols,pver)

   real(r8) :: flntc(pcols)         ! Clear sky lw flux at model top
   real(r8) :: flut(pcols)          ! Upward flux at top of model
   real(r8) :: flutc(pcols)         ! Upward Clear Sky flux at top of model
   real(r8) :: lwcf(pcols)          ! longwave cloud forcing

   real(r8) :: fln200(pcols)        ! net longwave flux interpolated to 200 mb
   real(r8) :: fln200c(pcols)       ! net clearsky longwave flux interpolated to 200 mb
   real(r8) :: flnr(pcols)          ! net longwave flux interpolated to tropopause

   real(r8) :: flnsc(pcols)         ! Clear sky lw flux at srf (up-down)
   real(r8) :: fldsc(pcols)         ! Clear sky lw flux at srf (down)

   real(r8) :: tot_cld_vistau(pcols,pver)   ! gbx water+ice cloud optical depth (only during day, night = fillvalue)
   real(r8) :: tot_icld_vistau(pcols,pver)  ! in-cld water+ice cloud optical depth (only during day, night = fillvalue)
   real(r8) :: liq_icld_vistau(pcols,pver)  ! in-cld liq cloud optical depth (only during day, night = fillvalue)
   real(r8) :: ice_icld_vistau(pcols,pver)  ! in-cld ice cloud optical depth (only during day, night = fillvalue)
   real(r8) :: snow_icld_vistau(pcols,pver) ! snow in-cloud visible sw optical depth for output on history files
   real(r8) :: grau_icld_vistau(pcols,pver) ! Graupel in-cloud visible sw optical depth for output on history files

   real(r8) :: cld_tau_cloudsim(pcols,pver)
   real(r8) :: aer_tau400(pcols,0:pver)
   real(r8) :: aer_tau550(pcols,0:pver)
   real(r8) :: aer_tau700(pcols,0:pver)

end type rad_out_t

type(var_desc_t) :: nextsw_cday_desc

!========================================================================================
contains
!========================================================================================

function rad_is_active()
  !-----------------------------------------------------------------------
  logical :: rad_is_active
  !-----------------------------------------------------------------------
  rad_is_active = .true.
end function rad_is_active

!================================================================================================

subroutine radiation_readnl(nlfile)

   ! Read radiation_nl namelist group.

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit

   character(len=*), intent(in) :: nlfile
   ! Local variables
   integer :: unitn, ierr
   integer :: dtime      ! timestep size
   character(len=*), parameter :: sub = 'radiation_readnl'

   namelist /radiation_nl/  k_h2o_file,  k_co2_file,  k_o2_file,  k_o3_file,  k_ch4_file,  k_c2h6_file,  &
        kh2o_mtckd_file,  kn2n2cia_file,  kn2h2cia_file,  kh2h2cia_file,  kco2co2cia_lw_file, &
        kco2co2cia_sw_file,  kco2h2cia_file,  kco2ch4cia_file,  ko2o2cia_file,  ko2n2cia_file, &
        ko2co2cia_file,  cldoptsL_file,  cldoptsI_file,  solar_file,  do_exo_rt_spectral, &
        do_exo_rt_clearsky,  do_exo_rt_optimize_bands,  do_exo_rt,  do_exo_atmconst, graupel_in_rad, &
        do_exo_synchronous,  do_exo_gw,  iradsw,  iradlw,  irad_always,  exo_rad_step

    !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'radiation_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, radiation_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(sub // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

   iradsw=exo_rad_step
   iradlw=exo_rad_step

   ! Broadcast namelist variables
   call mpi_bcast(k_h2o_file, cl, mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: k_h2o_file")
   call mpi_bcast(k_co2_file, cl, mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: k_co2_file")
   call mpi_bcast(k_o2_file, cl, mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: k_o2_file")
   call mpi_bcast(k_o3_file, cl, mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: k_o3_file")
   call mpi_bcast(k_ch4_file, cl, mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: k_ch4_file")
   call mpi_bcast(k_c2h6_file, cl, mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: k_c2h6_file")
   call mpi_bcast(kh2o_mtckd_file, cl, mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: kh2o_mtckd_file")
   call mpi_bcast(kn2n2cia_file, cl, mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: kn2n2cia_file")
   call mpi_bcast(kn2h2cia_file, cl, mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: kn2h2cia_file")
   call mpi_bcast(kh2h2cia_file, cl, mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: kh2h2cia_file")
   call mpi_bcast(kco2co2cia_lw_file, cl, mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: kco2co2cia_lw_file")
   call mpi_bcast(kco2co2cia_sw_file, cl, mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: kco2co2cia_sw_file")
   call mpi_bcast(kco2h2cia_file, cl, mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: kco2h2cia_file")
   call mpi_bcast(kco2ch4cia_file, cl, mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: kco2ch4cia_file")
   call mpi_bcast(ko2o2cia_file, cl, mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: ko2o2cia_file")
   call mpi_bcast(ko2n2cia_file, cl, mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: ko2n2cia_file")
   call mpi_bcast(ko2co2cia_file, cl, mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: ko2co2cia_file")
   call mpi_bcast(cldoptsL_file, cl, mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: cldoptsL_file")
   call mpi_bcast(cldoptsI_file, cl, mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: cldoptsI_file")
   call mpi_bcast(solar_file, cl, mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: solar_file")
   call mpi_bcast(iradsw, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: iradsw")
   call mpi_bcast(iradlw, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: iradlw")
   call mpi_bcast(irad_always, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: irad_always")
   call mpi_bcast(exo_rad_step, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: exo_rad_step")
   call mpi_bcast(do_exo_rt_spectral, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_exo_rt_spectral")
   call mpi_bcast(graupel_in_rad, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: graupel_in_rad")
   call mpi_bcast(do_exo_rt_clearsky, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_exo_rt_clearsky")
   call mpi_bcast(do_exo_rt_optimize_bands, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_exo_rt_optimize_bands")
   call mpi_bcast(do_exo_rt, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_exo_rt")
   call mpi_bcast(do_exo_atmconst, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_exo_atmconst")
   call mpi_bcast(do_exo_synchronous, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_exo_synchronous")
   call mpi_bcast(do_exo_gw, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_exo_gw")


   ! Convert iradsw, iradlw and irad_always from hours to timesteps if necessary
   dtime  = get_step_size()
   if (iradsw      < 0) iradsw      = nint((-iradsw     *3600._r8)/dtime)
   if (iradlw      < 0) iradlw      = nint((-iradlw     *3600._r8)/dtime)
   if (irad_always < 0) irad_always = nint((-irad_always*3600._r8)/dtime)

   !-----------------------------------------------------------------------
   ! Print runtime options to log.
   !-----------------------------------------------------------------------

   if (masterproc) then
      write(iulog,*) 'EXORT radiation scheme parameters:'
      write(iulog,10) trim(k_h2o_file), trim(k_co2_file), trim(k_o2_file), trim(k_o3_file), trim(k_ch4_file), &
        trim(k_c2h6_file), trim(kh2o_mtckd_file), &
        trim(kn2n2cia_file), trim(kn2h2cia_file), trim(kh2h2cia_file), trim(kco2co2cia_lw_file), trim(kco2co2cia_sw_file), &
        trim(kco2h2cia_file), trim(kco2ch4cia_file), trim(ko2o2cia_file), trim(ko2n2cia_file), trim(ko2co2cia_file), trim(cldoptsL_file), &
        trim(cldoptsI_file), trim(solar_file), iradsw,iradlw,irad_always,exo_rad_step,do_exo_rt_spectral, &
        do_exo_rt_clearsky,do_exo_rt_optimize_bands, do_exo_rt,do_exo_atmconst,do_exo_synchronous,do_exo_gw
   end if

10 format('  k_h2o_file        : ',                                a/, &
          '  k_co2_file        : ',                                a/, &
          '  k_o3_file         : ',                                a/, &
          '  k_o3_file         : ',                                a/, &
          '  k_ch4_file        : ',                                a/, &
          '  k_c2h6_file       : ',                                a/, &
          '  kh2o_mtckd_file   : ',                                a/, &
          '  kn2n2cia_file     : ',                                a/, &
          '  kn2h2cia_file     : ',                                a/, &
          '  kh2h2cia_file     : ',                                a/, &
          '  kco2co2cia_lw_file: ',                                a/, &
          '  kco2co2cia_sw_file: ',                                a/, &
          '  kco2h2cia_file    : ',                                a/, &
          '  kco2ch4cia_file   : ',                                a/, &
          '  ko2o2cia_file     : ',                                a/, &
          '  ko2n2cia_file     : ',                                a/, &
          '  ko2co2cia_file    : ',                                a/, &
          '  cldoptsL_file     : ',                                a/, &
          '  cldoptsI_file     : ',                                a/, &
          '  solar_file        : ',                                a/, &
          '  iradsw            :                                 ',i5/, &
          '  iradlw            :                                 ',i5/, &
          '  irad_always       :                                 ',i5/, &
          '  exo_rad_step      :                                 ',i5/, &
          '  do_exo_rt_spectral:                                 ',l5/, &
          '  do_exo_rt_clearsky:                                 ',l5/, &
          '  do_exo_rt_optimize_bands:                           ',l5/, &
          '  do_exo_rt               :                           ',l5/, &
          '  do_exo_atmconst         :                           ',l5/, &
          '  do_exo_synchronous      :                           ',l5/, &
          '  do_exo_gw               :                           ',l5/)

end subroutine radiation_readnl

!========================================================================================

subroutine radiation_register

   ! Register radiation fields in the physics buffer

   use radiation_data, only: rad_data_register
   use constituents,    only: cnst_add

   integer :: mm         ! dummy index for constituent add


   call pbuf_add_field('QRS' , 'global',dtype_r8,(/pcols,pver/), qrs_idx) ! shortwave radiative heating rate
   call pbuf_add_field('QRL' , 'global',dtype_r8,(/pcols,pver/), qrl_idx) ! longwave  radiative heating rate

   call pbuf_add_field('FSDS' , 'global',dtype_r8,(/pcols/), fsds_idx) ! Surface solar downward flux
   call pbuf_add_field('FSNS' , 'global',dtype_r8,(/pcols/), fsns_idx) ! Surface net shortwave flux
   call pbuf_add_field('FSNT' , 'global',dtype_r8,(/pcols/), fsnt_idx) ! Top-of-model net shortwave flux

   call pbuf_add_field('FLNS' , 'global',dtype_r8,(/pcols/), flns_idx) ! Surface net longwave flux
   call pbuf_add_field('FLNT' , 'global',dtype_r8,(/pcols/), flnt_idx) ! Top-of-model net longwave flux

   ! If the namelist has been configured for preserving the spectral fluxes, then create
   ! physics buffer variables to store the results.
   if (do_exo_rt_spectral) then
      call pbuf_add_field('SU'  , 'global',dtype_r8,(/pcols,pverp,ntot_wavlnrng/), su_idx) ! shortwave upward flux (per band)
      call pbuf_add_field('SD'  , 'global',dtype_r8,(/pcols,pverp,ntot_wavlnrng/), sd_idx) ! shortwave downward flux (per band)
      call pbuf_add_field('LU'  , 'global',dtype_r8,(/pcols,pverp,ntot_wavlnrng/), lu_idx) ! longwave upward flux (per band)
      call pbuf_add_field('LD'  , 'global',dtype_r8,(/pcols,pverp,ntot_wavlnrng/), ld_idx) ! longwave downward flux (per band)
   end if

   ! exort requires these inputs.

   call pbuf_add_field('REL' , 'global',dtype_r8,(/pcols,pver/), rel_idx) ! liquid water droplet size
   call pbuf_add_field('REI' , 'global',dtype_r8,(/pcols,pver/), rei_idx) ! liquid water droplet size
   call pbuf_add_field('CLD' , 'global',dtype_r8,(/pcols,pver/), cld_idx) ! cloud
   call pbuf_add_field('CICEWP', 'global',dtype_r8,(/pcols,pver/), cicewp_idx) !
   call pbuf_add_field('CLIQWP', 'global',dtype_r8,(/pcols,pver/), cliqwp_idx) !

   call rad_data_register()

   call cnst_add('CO2', 44._r8, 800._r8, 1.e-12_r8, mm, fixed_ubc=.false., &
        longname='CO2', readiv=.true., is_convtran1=.true.)

   call cnst_add('N2', 28._r8, 800._r8, 1.e-12_r8, mm, fixed_ubc=.false., &
        longname='N2', readiv=.true., is_convtran1=.true.)

   call pbuf_add_field('CO2' , 'global',dtype_r8,(/pcols,pver/), co2_idx) ! shortwave radiative heating rate
   call pbuf_add_field('N2' , 'global',dtype_r8,(/pcols,pver/), n2_idx) ! longwave  radiative heating rate

end subroutine radiation_register

!================================================================================================

subroutine radiation_init(pbuf2d)

   ! Initialize the radiation parameterization, add fields to the history buffer

   use phys_control,           only: phys_getopts
   use exo_init_ref,           only: init_ref
   use exo_model_specific,     only: init_model_specific
   use spectral_output,        only: addfld_spectral_intervals
   use time_manager,           only: is_first_step
   use radiation_data,         only: rad_data_init
   use exo_radiation_mod,      only: init_planck
!jt   use cloud_rad_props,  only: cloud_rad_props_init


   ! arguments
   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   ! local variables
   integer :: icall
   logical :: active_calls(0:N_DIAG)
   integer :: nstep                       ! current timestep number
   logical :: history_amwg                ! output the variables used by the AMWG diag package
   logical :: history_vdiag               ! output the variables used by the AMWG variability diag package
   logical :: history_budget              ! output tendencies and state variables for CAM4
                                          ! temperature, water vapor, cloud ice and cloud
                                          ! liquid budgets.
   integer :: history_budget_histfile_num ! output history file number for budget fields
   integer :: err

   integer :: dtime

   !-----------------------------------------------------------------------
   !jt not sure if setting camtop is needed
   camtop = 1

   ! initialize solar spectrum, clouds, and kcoeff
   call init_solar
   call init_kcoeff
   call init_cldopts
   call init_ref
   call init_model_specific
   call init_planck
   call addfld_spectral_intervals
   call rad_data_init(pbuf2d) ! initialize output fields for offline driver

! No Clouds yet
   cld_idx      = pbuf_get_index('CLD')
!!!   cldfsnow_idx = pbuf_get_index('CLDFSNOW',errcode=err)

   if (is_first_step()) then
      call pbuf_set_field(pbuf2d, qrl_idx, 0._r8)
   end if

   ! Surface components to get radiation computed today
   if (.not. is_first_restart_step()) then
      nextsw_cday = get_curr_calday()
   end if

   call phys_getopts(history_amwg_out   = history_amwg,    &
                     history_vdiag_out  = history_vdiag,   &
                     history_budget_out = history_budget,  &
                     history_budget_histfile_num_out = history_budget_histfile_num)

   ! "irad_always" is number of time steps to execute radiation continuously from start of
   ! initial OR restart run
   nstep = get_nstep()
   if (irad_always > 0) then
      nstep       = get_nstep()
      irad_always = irad_always + nstep
   end if

!     No clouds yet
!!$   call addfld('O3colAbove',    horiz_only,   'A', 'DU', 'Column O3 above model top', sampling_seq='rad_lwsw')
!!$
!!$   call addfld('TOT_CLD_VISTAU',  (/ 'lev' /), 'A',   '1', 'Total gbx cloud extinction visible sw optical depth', &
!!$                                                       sampling_seq='rad_lwsw', flag_xyfill=.true.)
!!$   call addfld('TOT_ICLD_VISTAU', (/ 'lev' /), 'A',  '1', 'Total in-cloud extinction visible sw optical depth', &
!!$                                                       sampling_seq='rad_lwsw', flag_xyfill=.true.)
!!$   call addfld('LIQ_ICLD_VISTAU', (/ 'lev' /), 'A',  '1', 'Liquid in-cloud extinction visible sw optical depth', &
!!$                                                       sampling_seq='rad_lwsw', flag_xyfill=.true.)
!!$   call addfld('ICE_ICLD_VISTAU', (/ 'lev' /), 'A',  '1', 'Ice in-cloud extinction visible sw optical depth', &
!!$                                                       sampling_seq='rad_lwsw', flag_xyfill=.true.)
!!$
!!$   if (cldfsnow_idx > 0) then
!!$      call addfld('SNOW_ICLD_VISTAU', (/ 'lev' /), 'A', '1', 'Snow in-cloud extinction visible sw optical depth', &
!!$                                                       sampling_seq='rad_lwsw', flag_xyfill=.true.)
!!$   endif

   ! get list of active radiation calls
   call rad_cnst_get_call_list(active_calls)

   ! Add shortwave radiation fields to history master field list.

   do icall = 0, N_DIAG

      if (active_calls(icall)) then

         call addfld('SOLIN'//diag(icall),    horiz_only,   'A', 'W/m2', 'Solar insolation', sampling_seq='rad_lwsw')

         call addfld('QRS'//diag(icall),      (/ 'lev' /),  'A', 'K/s',  'Solar heating rate', sampling_seq='rad_lwsw')
         call addfld('QRSC'//diag(icall),     (/ 'lev' /),  'A', 'K/s',  'Clearsky solar heating rate',                     &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSNT'//diag(icall),     horiz_only,   'A', 'W/m2', 'Net solar flux at top of model',                  &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSNTC'//diag(icall),    horiz_only,   'A', 'W/m2', 'Clearsky net solar flux at top of model',         &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSNTOA'//diag(icall),   horiz_only,   'A', 'W/m2', 'Net solar flux at top of atmosphere',             &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSNTOAC'//diag(icall),  horiz_only,   'A', 'W/m2', 'Clearsky net solar flux at top of atmosphere',    &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('SWCF'//diag(icall),     horiz_only,   'A', 'W/m2', 'Shortwave cloud forcing',                         &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSUTOA'//diag(icall),   horiz_only,   'A', 'W/m2', 'Upwelling solar flux at top of atmosphere',       &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSNIRTOA'//diag(icall), horiz_only,   'A', 'W/m2',                                                    &
                               'Net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', sampling_seq='rad_lwsw')
         call addfld('FSNRTOAC'//diag(icall), horiz_only,   'A', 'W/m2',                                                    &
                      'Clearsky net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', sampling_seq='rad_lwsw')
         call addfld('FSNRTOAS'//diag(icall), horiz_only,   'A', 'W/m2',                                                    &
                              'Net near-infrared flux (>= 0.7 microns) at top of atmosphere', sampling_seq='rad_lwsw')

         call addfld('FSN200'//diag(icall),   horiz_only,   'A', 'W/m2', 'Net shortwave flux at 200 mb',                    &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSN200C'//diag(icall),  horiz_only,   'A', 'W/m2', 'Clearsky net shortwave flux at 200 mb',           &
                                                                                 sampling_seq='rad_lwsw')

         call addfld('FSNR'//diag(icall),     horiz_only,   'A', 'W/m2', 'Net solar flux at tropopause',                    &
                                                                                 sampling_seq='rad_lwsw')

         call addfld('SOLL'//diag(icall),     horiz_only,   'A', 'W/m2', 'Solar downward near infrared direct  to surface', &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('SOLS'//diag(icall),     horiz_only,   'A', 'W/m2', 'Solar downward visible direct  to surface',       &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('SOLLD'//diag(icall),    horiz_only,   'A', 'W/m2', 'Solar downward near infrared diffuse to surface', &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('SOLSD'//diag(icall),    horiz_only,   'A', 'W/m2', 'Solar downward visible diffuse to surface',       &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSNS'//diag(icall),     horiz_only,   'A', 'W/m2', 'Net solar flux at surface',                       &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSNSC'//diag(icall),    horiz_only,   'A', 'W/m2', 'Clearsky net solar flux at surface',              &
                                                                                 sampling_seq='rad_lwsw')

         call addfld('FSDS'//diag(icall),     horiz_only,   'A', 'W/m2', 'Downwelling solar flux at surface',               &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSDSC'//diag(icall),    horiz_only,   'A', 'W/m2', 'Clearsky downwelling solar flux at surface',      &
                                                                                 sampling_seq='rad_lwsw')

         call addfld('FUS'//diag(icall),      (/ 'ilev' /), 'I', 'W/m2', 'Shortwave upward flux')
         call addfld('FDS'//diag(icall),      (/ 'ilev' /), 'I', 'W/m2', 'Shortwave downward flux')
         call addfld('FUSC'//diag(icall),     (/ 'ilev' /), 'I', 'W/m2', 'Shortwave clear-sky upward flux')
         call addfld('FDSC'//diag(icall),     (/ 'ilev' /), 'I', 'W/m2', 'Shortwave clear-sky downward flux')

         if (history_amwg) then
            call add_default('SOLIN'//diag(icall),   1, ' ')
            call add_default('QRS'//diag(icall),     1, ' ')
            call add_default('FSNT'//diag(icall),    1, ' ')
            call add_default('FSNTC'//diag(icall),   1, ' ')
            call add_default('FSNTOA'//diag(icall),  1, ' ')
            call add_default('FSNTOAC'//diag(icall), 1, ' ')
            call add_default('SWCF'//diag(icall),    1, ' ')
            call add_default('FSNS'//diag(icall),    1, ' ')
            call add_default('FSNSC'//diag(icall),   1, ' ')
            call add_default('FSUTOA'//diag(icall),  1, ' ')
            call add_default('FSDSC'//diag(icall),   1, ' ')
            call add_default('FSDS'//diag(icall),    1, ' ')
         endif

      end if
   end do

   if (scm_crm_mode) then
      call add_default('FUS     ', 1, ' ')
      call add_default('FUSC    ', 1, ' ')
      call add_default('FDS     ', 1, ' ')
      call add_default('FDSC    ', 1, ' ')
   endif

   ! Add longwave radiation fields to history master field list.

   do icall = 0, N_DIAG

      if (active_calls(icall)) then

         call addfld('QRL'//diag(icall),     (/ 'lev' /), 'A', 'K/s',  'Longwave heating rate', sampling_seq='rad_lwsw')
         call addfld('QRLC'//diag(icall),    (/ 'lev' /), 'A', 'K/s',  'Clearsky longwave heating rate',                   &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FLNT'//diag(icall),    horiz_only,  'A', 'W/m2', 'Net longwave flux at top of model',                &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FLNTC'//diag(icall),   horiz_only,  'A', 'W/m2', 'Clearsky net longwave flux at top of model',       &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FLNTCLR'//diag(icall), horiz_only,  'A', 'W/m2', 'Clearsky ONLY points net longwave flux at top of model',&
                                                                           sampling_seq='rad_lwsw')
         call addfld('FREQCLR'//diag(icall), horiz_only,  'A', 'Frac', 'Frequency of Occurrence of Clearsky',       &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FLUT'//diag(icall),    horiz_only,  'A', 'W/m2', 'Upwelling longwave flux at top of model',          &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FLUTC'//diag(icall),   horiz_only,  'A', 'W/m2', 'Clearsky upwelling longwave flux at top of model', &
                                                                           sampling_seq='rad_lwsw')
         call addfld('LWCF'//diag(icall),    horiz_only,  'A', 'W/m2', 'Longwave cloud forcing', sampling_seq='rad_lwsw')

         call addfld('FLN200'//diag(icall),  horiz_only,  'A', 'W/m2', 'Net longwave flux at 200 mb',                      &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FLN200C'//diag(icall), horiz_only,  'A', 'W/m2', 'Clearsky net longwave flux at 200 mb',             &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FLNR'//diag(icall),    horiz_only,  'A', 'W/m2', 'Net longwave flux at tropopause',                  &
                                                                           sampling_seq='rad_lwsw')

         call addfld('FLNS'//diag(icall),    horiz_only,  'A', 'W/m2', 'Net longwave flux at surface',                     &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FLNSC'//diag(icall),   horiz_only,  'A', 'W/m2', 'Clearsky net longwave flux at surface',            &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FLDS'//diag(icall),    horiz_only,  'A', 'W/m2', 'Downwelling longwave flux at surface',             &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FLDSC'//diag(icall),   horiz_only,  'A', 'W/m2', 'Clearsky Downwelling longwave flux at surface',    &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FUL'//diag(icall),     (/ 'ilev' /),'I', 'W/m2', 'Longwave upward flux')
         call addfld('FDL'//diag(icall),     (/ 'ilev' /),'I', 'W/m2', 'Longwave downward flux')
         call addfld('FULC'//diag(icall),    (/ 'ilev' /),'I', 'W/m2', 'Longwave clear-sky upward flux')
         call addfld('FDLC'//diag(icall),    (/ 'ilev' /),'I', 'W/m2', 'Longwave clear-sky downward flux')
         call addfld('CO2'//diag(icall),     (/ 'lev' /), 'A', 'Kg/kg',  'CO2 mixing ratio')
         call addfld('N2'//diag(icall),     (/ 'lev' /), 'A', 'Kg/Kg',  'N2 mixing ratio')

         if (history_amwg) then
            call add_default('QRL'//diag(icall),   1, ' ')
            call add_default('FLNT'//diag(icall),  1, ' ')
            call add_default('FLNTC'//diag(icall), 1, ' ')
            call add_default('FLNTCLR'//diag(icall), 1, ' ')
            call add_default('FREQCLR'//diag(icall), 1, ' ')
            call add_default('FLUT'//diag(icall),  1, ' ')
            call add_default('FLUTC'//diag(icall), 1, ' ')
            call add_default('LWCF'//diag(icall),  1, ' ')
            call add_default('FLNS'//diag(icall),  1, ' ')
            call add_default('FLNSC'//diag(icall), 1, ' ')
            call add_default('FLDS'//diag(icall),  1, ' ')
         endif

      end if
   end do

   call addfld('EMIS', (/ 'lev' /), 'A', '1', 'Cloud longwave emissivity')

   if (scm_crm_mode) then
      call add_default ('FUL     ', 1, ' ')
      call add_default ('FULC    ', 1, ' ')
      call add_default ('FDL     ', 1, ' ')
      call add_default ('FDLC    ', 1, ' ')
   endif

   ! Heating rate needed for d(theta)/dt computation
   call addfld ('HR',(/ 'lev' /), 'A','K/s','Heating rate needed for d(theta)/dt computation')

   if ( history_budget .and. history_budget_histfile_num > 1 ) then
      call add_default ('QRL     ', history_budget_histfile_num, ' ')
      call add_default ('QRS     ', history_budget_histfile_num, ' ')
   end if

   if (history_vdiag) then
      call add_default('FLUT', 2, ' ')
      call add_default('FLUT', 3, ' ')
   end if

end subroutine radiation_init

!===============================================================================

function radiation_do(op, timestep)

   ! Return true if the specified operation is done this timestep.

   character(len=*), intent(in) :: op             ! name of operation
   integer, intent(in), optional:: timestep
   logical                      :: radiation_do   ! return value

   integer :: nstep
   !---------------------------------------------------------------------------
   radiation_do = .false.
   !------------------------------------------------------------------------
   !
   ! Start Code
   !
   if (present(timestep)) then
      nstep = timestep
   else
      nstep = get_nstep()
   end if

   select case (op)

   case ('sw','lw','exort') ! do a shortwave heating calc this timestep?
      radiation_do = nstep == 0  .or.  exo_rad_step == 1                     &
           .or. (mod(nstep,exo_rad_step) == 0  .and.  nstep /= 1)

   case default
      call endrun('radiation_do: unknown operation:'//op//' use "lw","sw",or,"exort"')

   end select


 end function radiation_do

!========================================================================================

subroutine radiation_define_restart(file)

   ! define variables to be written to restart file

   ! arguments
   type(file_desc_t), intent(inout) :: file

   ! local variables
   integer :: ierr
   !----------------------------------------------------------------------------

   call pio_seterrorhandling(File, PIO_BCAST_ERROR)

   ierr = pio_def_var(File, 'nextsw_cday', pio_double, nextsw_cday_desc)
   ierr = pio_put_att(File, nextsw_cday_desc, 'long_name', 'future radiation calday for surface models')

end subroutine radiation_define_restart

!===============================================================================

subroutine radiation_write_restart(file)

   ! write variables to restart file

   ! arguments
   type(file_desc_t), intent(inout) :: file

   ! local variables
   integer :: ierr
   !----------------------------------------------------------------------------

   ierr = pio_put_var(File, nextsw_cday_desc, (/ nextsw_cday /))

end subroutine radiation_write_restart

!===============================================================================

subroutine radiation_read_restart(file)

   ! read variables from restart file

   ! arguments
   type(file_desc_t), intent(inout) :: file

   ! local variables

   integer :: err_handling
   integer :: ierr
   real(r8) :: temp_var

   type(var_desc_t) :: vardesc
   !----------------------------------------------------------------------------

   ierr = pio_inq_varid(File, 'nextsw_cday', vardesc)
   ierr = pio_get_var(File, vardesc, temp_var)
   nextsw_cday = temp_var

end subroutine radiation_read_restart

!===============================================================================

subroutine radiation_tend( &
   state, ptend, pbuf, cam_out, cam_in, net_flx, rd_out)

   !-----------------------------------------------------------------------
   !
   ! Driver for radiation computation.
   !
   ! Revision history:
   !-----------------------------------------------------------------------
    ! mars_radiative_tend: Run the radiative process
    !=========================================================================
    use cam_control_mod, only: lambm0, obliqr, eccen, mvelpp
    use rad_constituents,    only: N_DIAG, rad_cnst_get_call_list, rad_cnst_get_gas, rad_cnst_out
    use radheat,             only: radheat_tend
    use exo_radiation_mod,    only: aerad_driver
    use shr_orb_mod,  only: shr_orb_decl, shr_orb_cosz
    !
    ! Passed Variables
    !------------------
    type(physics_state),intent(in)   :: state
    type(physics_ptend),intent(out)  :: ptend
    type(physics_buffer_desc), pointer      :: pbuf(:)
    type(cam_in_t),     intent(in)          :: cam_in
    type(cam_out_t),    intent(inout)       :: cam_out
    real(r8),           intent(out)         :: net_flx(:)
    type(rad_out_t), target, optional, intent(out) :: rd_out

   ! Local variables
   type(rad_out_t), pointer :: rd  ! allow rd_out to be optional by allocating a local object
                                   ! if the argument is not present
   logical  :: write_output

    real(r8):: T           (state%ncol,pver) ! T temporary
    real(r8):: qv          (state%ncol,pver) ! Q temporary
    real(r8):: dtdt_heating(state%ncol,pver) ! Longwave heating tendency K/s
    real(r8):: dtdt_solar  (state%ncol,pver) ! Shortwave heating tendency K/s
    real(r8):: Tsfc        (state%ncol)      ! Surface T
    real(r8):: Qsfc        (state%ncol)      ! Surface Q (saturated)
    logical :: lq(pcnst)                     ! Calc tendencies?
    integer :: lchnk                         ! chunk identifier
    integer :: ncol                          ! number of atmospheric columns
    integer :: k                             ! loop index
    real(r8):: ftem(pcols,pver)              ! Temporary workspace for outfld variables
    real(r8):: ftem2(pcols,pver)             ! Temporary workspace for outfld variables
    real(r8), pointer :: cicewp(:,:)   ! in-cloud cloud ice water path (from param_cldoptics_calc)
    real(r8), pointer :: cliqwp(:,:)   ! in-cloud cloud liquid water path (from param_cldoptics_calc)
    real(r8), pointer :: cldfrc(:,:)   ! cloud fraction
    real(r8), pointer :: qrs(:,:)      ! shortwave radiative heating rate
    real(r8), pointer :: qrl(:,:)      ! longwave  radiative heating rate
    real(r8), pointer :: rel(:,:)      ! effective liquid drop radius
    real(r8), pointer :: rei(:,:)      ! effective ice particle radius
    real(r8), pointer :: fsds(:)  ! Surface solar down flux
    real(r8), pointer :: fsns(:)  ! Surface solar absorbed flux
    real(r8), pointer :: fsnt(:)  ! Net column abs solar flux at model top
    real(r8), pointer :: flns(:)  ! Srf longwave cooling (up-down) flux
    real(r8), pointer :: flnt(:)  ! Net outgoing lw flux at model top

    real(r8), pointer, dimension(:,:,:) :: su => NULL()  ! shortwave spectral flux up
    real(r8), pointer, dimension(:,:,:) :: sd => NULL()  ! shortwave spectral flux down
    real(r8), pointer, dimension(:,:,:) :: lu => NULL()  ! longwave  spectral flux up
    real(r8), pointer, dimension(:,:,:) :: ld => NULL()  ! longwave  spectral flux down
    real(r8) :: calday                        ! current calendar day
    real(r8), dimension(pcols) :: clat        ! current latitudes(radians)
    real(r8), dimension(pcols) :: clon        ! current longitudes(radians)
    real(r8), dimension(pcols) ::  coszrs     ! Cosine solar zenith angle
    integer  :: nstep                         ! current timestep number
    logical  :: conserve_energy = .true.      ! flag to carry (QRS,QRL)*dp across time steps
    real(r8) :: frac_day
    real(r8) :: day_in_year
    logical :: do_exo_rad
    real(r8) :: eccf                            ! Earth/sun distance factor
    real(r8) :: delta                           ! Solar declination angle
    !------------------------------------
    integer :: ext_tslas_tog
    integer :: ext_tshadow_tog

    real(r8) :: ext_solar_azm_ang
    real(r8) :: ext_tazm_ang
    real(r8) :: ext_tslope_ang

    real(r8) :: ext_msdist
    real(r8) :: ext_rtgt
    ! null cloud place holders, for clear sky calculation
    real(r8), dimension(pcols,pver) :: cicewp_zero
    real(r8), dimension(pcols,pver) :: cliqwp_zero
    real(r8), dimension(pcols,pver) :: cldfrc_zero


    !the following would have dimension nazm_tshadow if shadows were taken into effect
    integer, parameter :: ext_nazm_tshadow = 1                  ! take shadows into effect
    real(r8), dimension(ext_nazm_tshadow) :: ext_cosz_horizon   ! cos of zenith angle of horizon
    real(r8), dimension(ext_nazm_tshadow) :: ext_TCx_obstruct
    real(r8), dimension(ext_nazm_tshadow) :: ext_TCz_obstruct
    !------------------------------------
    real(r8), pointer, dimension(:,:) :: h2ommr   ! h2o   mass mixing ratio
    real(r8), pointer, dimension(:,:) :: co2mmr   ! co2   mass mixing ratio
!jt    real(r8), pointer, dimension(:,:) :: ch4mmr   ! ch4   mass mixing ratio
    real(r8), dimension(pcols,pver) :: ch4mmr     ! h2    mass mixing ratio
    real(r8), dimension(pcols,pver) :: h2mmr      ! h2    mass mixing ratio
    real(r8), dimension(pcols,pver) :: n2mmr      ! n2    mass mixing ratio
    real(r8), dimension(pcols,pver) :: c2h6mmr    ! c2h6   mass mixing ratio

    !o2/o3--not sure where these will come from so will need to check!
    !real(r8), pointer, dimension(:,:) :: o2mmr   ! o2   mass mixing ratio
    !real(r8), pointer, dimension(:,:) :: o3mmr   ! o3   mass mixing ratio
    real(r8), dimension(pcols,pver) :: o2mmr      ! h2    mass mixing ratio
    real(r8), dimension(pcols,pver) :: o3mmr      ! h2    mass mixing ratio

    real(r8) :: rdtime
    real(r8) :: vis_dir
    real(r8) :: vis_dif
    real(r8) :: nir_dir
    real(r8) :: nir_dif
    real(r8) :: sol_toa
    real(r8), dimension(pver) :: sw_dTdt
    real(r8), dimension(pver) :: lw_dTdt
    real(r8), dimension(pverp) :: sw_upflux
    real(r8), dimension(pverp) :: sw_dnflux
    real(r8), dimension(pverp) :: lw_upflux
    real(r8), dimension(pverp) :: lw_dnflux
    real(r8), dimension(pverp,ntot_wavlnrng) :: sw_upflux_spec
    real(r8), dimension(pverp,ntot_wavlnrng) :: sw_dnflux_spec
    real(r8), dimension(pverp,ntot_wavlnrng) :: lw_upflux_spec
    real(r8), dimension(pverp,ntot_wavlnrng) :: lw_dnflux_spec
    real(r8), dimension(pcols) :: fsdtoa        ! Incoming solar flux at TOA
    real(r8), dimension(pcols) :: fsntoa        ! Net solar flux at TOA
    real(r8), dimension(pcols) :: fsntoac       ! Clear sky net solar flux at TOA
    real(r8), dimension(pcols) :: fsnirt        ! Near-IR flux absorbed at toa
    real(r8), dimension(pcols) :: fsnrtc        ! Clear sky near-IR flux absorbed at toa
    real(r8), dimension(pcols) :: fsnirtsq      ! Near-IR flux absorbed at toa >= 0.7 microns
    real(r8), dimension(pcols) :: fsntc         ! Clear sky total column abs solar flux
    real(r8), dimension(pcols) :: fsnsc         ! Clear sky surface abs solar flux
    real(r8), dimension(pcols) :: fsdsc         ! Clear sky surface downwelling solar flux
    real(r8), dimension(pcols) :: flut          ! Upward flux at top of model
    real(r8), dimension(pcols) :: lwcf          ! longwave cloud forcing
    real(r8), dimension(pcols) :: swcf          ! shortwave cloud forcing
    real(r8), dimension(pcols) :: flutc         ! Upward Clear Sky flux at top of model
    real(r8), dimension(pcols) :: flntc         ! Clear sky lw flux at model top
    real(r8), dimension(pcols) :: flnsc         ! Clear sky lw flux at srf (up-down)
    real(r8), dimension(pcols,pverp) :: fns     ! net shortwave flux
    real(r8), dimension(pcols,pverp) :: fcns    ! net clear-sky shortwave flux
    real(r8), dimension(pcols,pverp) :: fsn     ! net shortwave flux
    real(r8), dimension(pcols,pverp) :: fln     ! net longwave flux
    real(r8), dimension(pcols,pverp) :: fnl     ! net longwave flux
    real(r8), dimension(pcols,pverp) :: fcnl    ! net clear-sky longwave flux
    integer :: itim_old
    integer :: i
    ! summed fluxes
    ! broadband fluxes
    real(r8), dimension(pcols,pverp) :: lwup_rad
    real(r8), dimension(pcols,pverp) :: lwdown_rad
    real(r8), dimension(pcols,pverp) :: swup_rad
    real(r8), dimension(pcols,pverp) :: swdown_rad
    !jt ts hack
    real(r8), dimension(pver) :: hack_ts

    !------------------------------------------------------------------------
    !
    ! Start Code
    !

    ! initialize individual parameterization tendencies
    !---------------------------------------------------
    lchnk       = state%lchnk
    ncol        = state%ncol


    if (present(rd_out)) then
       rd => rd_out
       write_output = .false.
    else
       allocate(rd)
       write_output=.true.
    end if

    lq(:) = .false.
    call physics_ptend_init(ptend, state%psetcols, 'Mars radiative_tend',   &
         ls=.true., lu=.false., lv=.false., lq=lq)
    calday = get_curr_calday()

    itim_old = pbuf_old_tim_idx()
    nullify(cldfrc)
    call pbuf_get_field(pbuf, cld_idx,    cldfrc,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
    call pbuf_get_field(pbuf, qrs_idx, qrs)
    call pbuf_get_field(pbuf, qrl_idx, qrl)
    call pbuf_get_field(pbuf, rel_idx, rel)
    call pbuf_get_field(pbuf, rei_idx, rei)
    call pbuf_get_field(pbuf, cicewp_idx, cicewp)
    call pbuf_get_field(pbuf, cliqwp_idx, cliqwp)
    call pbuf_get_field(pbuf, fsnt_idx, fsnt)
    call pbuf_get_field(pbuf, fsds_idx, fsds)
    call pbuf_get_field(pbuf, fsns_idx, fsns)
    call pbuf_get_field(pbuf, flns_idx, flns)
    call pbuf_get_field(pbuf, flnt_idx, flnt)

    if (do_exo_rt_spectral) then
       call pbuf_get_field(pbuf, su_idx, su)
       call pbuf_get_field(pbuf, sd_idx, sd)
       call pbuf_get_field(pbuf, lu_idx, lu)
       call pbuf_get_field(pbuf, ld_idx, ld)
    end if

!!$       call phys_getopts(microp_scheme_out=microp_pgk)
!!$       rk_clouds = microp_pgk == 'RK'
!!$       mg_clouds = microp_pgk == 'MG'
!!$
!!$       if (mg_clouds) then
!!$          ! convert cloud condensate units
!!$          ! mg clouds deal in kg/m2, whereas rk clouds are in g/m2
!!$          ! the radiation scheme takes g/m2
!!$          cicewp(:,:) = cicewp(:,:)*1000.
!!$          cliqwp(:,:) = cliqwp(:,:)*1000.
!!$       endif

    !
    ! Cosine solar zenith angle for current time step
    !
    call get_rlat_all_p(lchnk, ncol, clat)
    call get_rlon_all_p(lchnk, ncol, clon)

    ! Wolf, length of day scaling for zenith angle calculation
    call get_curr_calday_rotation(frac_day, day_in_year)
    call zenith_rotation (frac_day, calday, clat, clon, coszrs, ncol)
    !call zenith (calday, clat, clon, coszrs, ncol)

    ! calculate Earth-Sun distance factor, scaled by eccentricity factor
    ! in the next version the orbital details will be more heavily modulated.
    ! NOTE:  SHR_CONST_MSDIST = normalized planet-star distance squared
    ! do standard orbital calculation, for determining sflux_frac
    call shr_orb_decl(calday, eccen, mvelpp, lambm0, obliqr, delta, eccf)
    ext_msdist=1.0/eccf

    !! Main calculation starts here !!
    !do_exo_rad determines if RT called on a given timesteps
    do_exo_rad = radiation_do('exort')

    if (do_exo_rad) then
       !write(*,*) "inside RT if statement"

       !-------for now, set these topography angles to 0 (ie, no topography blocking sun)
       !-------should be moved to inside do loop if we want to use them
       !-------azimuth can be calculated from cos(azim) = (cos(hr_ang)*cos(dec)*sin(lat)-sin(dec)*cos(lat))/cos(elev)
       !-------the azimuth is the angle east of south when the hour angle, h, is negative (morning)
       !-------and the angle west of south when the hour angle, h, is positive (afternoon).
       ext_solar_azm_ang = 0.     !solar azimuthal angle? should not be zero
       ext_tazm_ang = 0.
       ext_tslope_ang = 0.
       ext_tslas_tog = 0
       ext_tshadow_tog = 1        ! toggle shadowing
       ext_cosz_horizon(:) = 0.
       ext_TCx_obstruct(:) = 0.
       ext_TCz_obstruct(:) = 0.
       ext_rtgt = 1.

       sw_dTdt(:) = 0.    ! Initialize heating rate arrays
       lw_dTdt(:) = 0.    !

       lw_dnflux(:) = 0.   !
       lw_upflux(:) = 0.   ! Initialize entire arrays for summing below
       sw_upflux(:) = 0.   !
       sw_dnflux(:) = 0.   !

!!$      lwup_rad_spec(:,:,:) = 0.
!!$      lwdown_rad_spec(:,:,:) = 0.
!!$      swup_rad_spec(:,:,:) = 0.
!!$      swdown_rad_spec(:,:,:) = 0.

       if (do_exo_rt_spectral) then
          lu(:,:,:) = 0.
          ld(:,:,:) = 0.
          su(:,:,:) = 0.
          sd(:,:,:) = 0.
       end if

       qrs(:,:) = 0.
       qrl(:,:) = 0.
       lwup_rad(:,:) = 0.
       lwdown_rad(:,:) = 0.
       swup_rad(:,:) = 0.
       swdown_rad(:,:) = 0.

       nstep = get_nstep()

       ! Native CAM functions; returns pointer to mass mixing ratio for the gas specified
       ! gases that exist in the CESM physics buffer
       call rad_cnst_get_gas(0,'CO2', state, pbuf,  co2mmr)
!!$       call rad_cnst_get_gas(0,'CH4', state, pbuf,  ch4mmr)
       call rad_cnst_get_gas(0,'H2O', state, pbuf,  h2ommr) !H2O specific humidity
!!$       call rad_cnst_get_gas(0,'O2',  state, pbuf,  o2mmr)  !does CESM have this one?
!!$       call rad_cnst_get_gas(0,'O3',  state, pbuf,  o3mmr)  !does CESM have this one?

       ! well mixed species from exoplanet_mod.F90
       ! not in CESM physics buffer
       ch4mmr(:,:)   = exo_ch4mmr
       n2mmr(:,:)   = exo_n2mmr
       h2mmr(:,:)   = exo_h2mmr
       c2h6mmr(:,:) = exo_c2h6mmr
       o2mmr(:,:)   = 0._r8   !exo_o2mmr
       o3mmr(:,:)   = 0._r8   !exo_o3mmr   !use if CESM doesnt have these guys

       ! Do a parallel clearsky radiative calculation so we can calculate cloud forcings
       ! Setting do_exo_rt_clearsky to true, slows the code dramatically, use wisely and sparingly
       if (do_exo_rt_clearsky) then

          ! set clouds to zero everywhere
          cicewp_zero(:,:) = 0.0_r8
          cliqwp_zero(:,:) = 0.0_r8
          cldfrc_zero(:,:) = 0.0_r8

          !jt temporary hack for surface temperature - we need land to build to provide this.
          hack_ts(:) = 190.0_r8
          do i = 1, ncol

             call aerad_driver(h2ommr(i,:), co2mmr(i,:), &
                  ch4mmr(i,:), c2h6mmr(i,:), &
                  h2mmr(i,:),  n2mmr(i,:), o3mmr(i,:), o2mmr(i,:), &
                  cicewp_zero(i,:), cliqwp_zero(i,:), cldfrc_zero(i,:), &
                  rei(i,:), rel(i,:), &
!jt                  cam_in%ts(i), state%ps(i), state%pmid(i,:), &
                  hack_ts(i), state%ps(i), state%pmid(i,:), &
                  state%pdel(i,:), state%pdeldry(i,:), state%t(i,:), state%pint(i,:), state%pintdry(i,:), &
                  coszrs(i), ext_msdist, &
                  cam_in%asdir(i), cam_in%aldir(i), &
                  cam_in%asdif(i), cam_in%aldif(i), &
                  ext_rtgt, ext_solar_azm_ang, ext_tazm_ang, ext_tslope_ang,  &
                  ext_tslas_tog, ext_tshadow_tog, ext_nazm_tshadow, ext_cosz_horizon , &
                  ext_TCx_obstruct, ext_TCz_obstruct, state%zi(i,:), &
                  sw_dTdt, lw_dTdt, lw_dnflux, lw_upflux, sw_upflux, sw_dnflux,  &
                  lw_dnflux_spec, lw_upflux_spec, sw_upflux_spec, sw_dnflux_spec, &
                  vis_dir, vis_dif, nir_dir, nir_dif, sol_toa )


             ftem(i,:) = sw_dTdt(:)
             ftem2(i,:) = lw_dTdt(:)

             lwup_rad(i,:) = lw_upflux(:)
             lwdown_rad(i,:) = lw_dnflux(:)
             swup_rad(i,:) = sw_upflux(:)
             swdown_rad(i,:) = sw_dnflux(:)

             if (do_exo_rt_spectral) then
!!$            lwup_rad_spec(i,:,:)   = lw_upflux_spec(:,:)
!!$            lwdown_rad_spec(i,:,:) = lw_dnflux_spec(:,:)
!!$            swup_rad_spec(i,:,:)   = sw_upflux_spec(:,:)
!!$            swdown_rad_spec(i,:,:) = sw_dnflux_spec(:,:)

                lu(i,:,:)   = lw_upflux_spec(:,:)
                ld(i,:,:)   = lw_dnflux_spec(:,:)
                su(i,:,:)   = sw_upflux_spec(:,:)
                sd(i,:,:)   = sw_dnflux_spec(:,:)
             endif

             ! Fluxes sent to land model
             ! Note these values are overwritten by the full-sky rt calc
             cam_out%sols(i) = vis_dir
             cam_out%soll(i) = nir_dir
             cam_out%solsd(i) = vis_dif
             cam_out%solld(i) = nir_dif
             cam_out%flwds(i) = lw_dnflux(pverp)

          enddo   ! ncol loop

          fsn(:,:) = swdown_rad(:,:) - swup_rad(:,:)
          fln(:,:) = lwup_rad(:,:) - lwdown_rad(:,:)
          fsns(:) = fsn(:,pverp)
          flns(:) = fln(:,pverp)
          fsnt(:) = fsn(:,1)
          flnt(:) = fln(:,1)
          fsds(:) = swdown_rad(:,pverp)
          qrs(:ncol,:pver) = ftem(:ncol,:pver)
          qrl(:ncol,:pver) = ftem2(:ncol,:pver)

          !jt             call outfld('QRSC     ',qrs*SHR_CONST_CDAY  , pcols,lchnk)    ! [K/day]
!!$          call outfld('QRSC     ',qrs*exo_diurnal  , pcols,lchnk)    ! [K/day]
!!$          call outfld('FSDSC    ',fsds  ,pcols,lchnk)
!!$          call outfld('FSNTC    ',fsnt  ,pcols,lchnk)
!!$          call outfld('FSNSC    ',fsns  ,pcols,lchnk)
!!$          !jt             call outfld('QRLC     ',qrl*SHR_CONST_CDAY   ,pcols,lchnk)    ! [K/day]
!!$          call outfld('QRLC     ',qrl*exo_diurnal   ,pcols,lchnk)    ! [K/day]
!!$          call outfld('FLNTC    ',flnt  ,pcols,lchnk)
!!$          call outfld('FLUTC    ',lwup_rad(:,2)  ,pcols,lchnk)
!!$          call outfld('FLNSC    ',flns  ,pcols,lchnk)
!!$          call outfld('FULC     ',lwup_rad, pcols, lchnk)
!!$          call outfld('FDLC     ',lwdown_rad, pcols, lchnk)
!!$          call outfld('FUSC     ',swup_rad, pcols, lchnk)
!!$          call outfld('FDSC     ',swdown_rad, pcols, lchnk)
!jt          call outfld('SOLSC    ',cam_out%sols  ,pcols,lchnk)
!jt          call outfld('SOLLC    ',cam_out%soll  ,pcols,lchnk)
!jt          call outfld('SOLSDC   ',cam_out%solsd ,pcols,lchnk)
!jt          call outfld('SOLLDC   ',cam_out%solld ,pcols,lchnk)

!!$             call set_diags_clearsky()
!!$             call radiation_output_clearsky(lchnk, ncol, icall, rd, pbuf, cam_out)
          if (do_exo_rt_spectral) call outfld_spectral_flux_clearsky(lchnk, ld, lu, su, sd)

       endif  ! (do_exo_rt_clearsky)

       ! Do Column Radiative transfer calculation WITH clouds.

       ! For now set rel/rei to a reasonable value and init clouds to 0
          rel=3.0_r8
          rei=3.0_r8
          cicewp=0._r8
          cliqwp=0._r8
          cldfrc=0._r8

       do i = 1, ncol

          call aerad_driver(h2ommr(i,:), co2mmr(i,:), &
               ch4mmr(i,:), c2h6mmr(i,:), &
               h2mmr(i,:),  n2mmr(i,:), o3mmr(i,:), o2mmr(i,:), &
               cicewp(i,:), cliqwp(i,:), cldfrc(i,:), &
               rei(i,:), rel(i,:), &
               hack_ts(i), state%ps(i), state%pmid(i,:), &
!jt               cam_in%ts(i), state%ps(i), state%pmid(i,:), &
               state%pdel(i,:), state%pdeldry(i,:), state%t(i,:), state%pint(i,:), state%pintdry(i,:), &
               coszrs(i), ext_msdist, &
               cam_in%asdir(i), cam_in%aldir(i), &
               cam_in%asdif(i), cam_in%aldif(i), &
               ext_rtgt, ext_solar_azm_ang, ext_tazm_ang, ext_tslope_ang,  &
               ext_tslas_tog, ext_tshadow_tog, ext_nazm_tshadow, ext_cosz_horizon,  &
               ext_TCx_obstruct, ext_TCz_obstruct, state%zi(i,:), &
               sw_dTdt, lw_dTdt, lw_dnflux, lw_upflux, sw_upflux, sw_dnflux,  &
               lw_dnflux_spec, lw_upflux_spec, sw_upflux_spec, sw_dnflux_spec, &
               vis_dir, vis_dif, nir_dir, nir_dif, sol_toa )


          ftem(i,:) = sw_dTdt(:)
          ftem2(i,:) = lw_dTdt(:)

          lwup_rad(i,:) = lw_upflux(:)
          lwdown_rad(i,:) = lw_dnflux(:)
          swup_rad(i,:) = sw_upflux(:)
          swdown_rad(i,:) = sw_dnflux(:)
          fsdtoa(i) = sol_toa

          if (do_exo_rt_spectral) then
!!$           lwup_rad_spec(i,:,:)   = lw_upflux_spec(:,:)
!!$           lwdown_rad_spec(i,:,:) = lw_dnflux_spec(:,:)
!!$           swup_rad_spec(i,:,:)   = sw_upflux_spec(:,:)
!!$           swdown_rad_spec(i,:,:) = sw_dnflux_spec(:,:)
             lu(i,:,:)   = lw_upflux_spec(:,:)
             ld(i,:,:)   = lw_dnflux_spec(:,:)
             su(i,:,:)   = sw_upflux_spec(:,:)
             sd(i,:,:)   = sw_dnflux_spec(:,:)
          endif

          ! Fluxes sent to land model
          ! Note these values overwrite those use in do_exo_rt_clearsky calc
          cam_out%sols(i) = vis_dir
          cam_out%soll(i) = nir_dir
          cam_out%solsd(i) = vis_dif
          cam_out%solld(i) = nir_dif
          cam_out%flwds(i) = lw_dnflux(pverp)

       enddo   ! ncol loop

       fsn(:,:) = swdown_rad(:,:) - swup_rad(:,:)
       fln(:,:) = lwup_rad(:,:) - lwdown_rad(:,:)
       fsns(:) = fsn(:,pverp)
       flns(:) = fln(:,pverp)
       fsnt(:) = fsn(:,1)
       flnt(:) = fln(:,1)
       fsds(:) = swdown_rad(:,pverp)
       qrs(:ncol,:pver) = ftem(:ncol,:pver)
       qrl(:ncol,:pver) = ftem2(:ncol,:pver)

       !jt          call outfld('QRS     ',qrs*SHR_CONST_CDAY  , pcols,lchnk)    ! [K/day]
       call outfld('QRS     ',qrs*exo_diurnal  , pcols,lchnk)    ! [K/day]
!!$       call outfld('FSDS    ',fsds  ,pcols,lchnk)
!!$       call outfld('FSDTOA  ',fsdtoa  ,pcols,lchnk)
!!$       call outfld('FSNT    ',fsnt  ,pcols,lchnk)
!!$       call outfld('FSNS    ',fsns  ,pcols,lchnk)
!!$       !jt          call outfld('QRL     ',qrl*SHR_CONST_CDAY   ,pcols,lchnk)    ! [K/day]
       call outfld('QRL     ',qrl*exo_diurnal   ,pcols,lchnk)    ! [K/day]
!!$       call outfld('FLNT    ',flnt  ,pcols,lchnk)
!!$       call outfld('FLUT    ',lwup_rad(:,1)  ,pcols,lchnk)   ! was 2
!!$       call outfld('FLNS    ',flns  ,pcols,lchnk)
!!$       call outfld('FUL     ',lwup_rad, pcols, lchnk)
!!$       call outfld('FDL     ',lwdown_rad, pcols, lchnk)
!!$       call outfld('FUS     ',swup_rad, pcols, lchnk)
!!$       call outfld('FDS     ',swdown_rad, pcols, lchnk)
!!$       call outfld('SOLS    ',cam_out%sols  ,pcols,lchnk)
!!$       call outfld('SOLL    ',cam_out%soll  ,pcols,lchnk)
!!$       call outfld('SOLSD   ',cam_out%solsd ,pcols,lchnk)
!!$       call outfld('SOLLD   ',cam_out%solld ,pcols,lchnk)

       !jt          call set_diags_fullsky()
       !jt          call radiation_output_fullsky(lchnk, ncol, icall, rd, pbuf, cam_out)
       if (do_exo_rt_spectral) call outfld_spectral_flux_fullsky(lchnk, ld, lu, su, sd)

    else ! if (do_exo_rad) then

       ! Radiative flux calculations not done.  The quantity Q*dp is carried by the
       ! physics buffer across timesteps.  It must be converted to Q (dry static energy
       ! tendency) before being passed to radheat_tend.
       qrs(:ncol,:) = qrs(:ncol,:) / state%pdel(:ncol,:)
       qrl(:ncol,:) = qrl(:ncol,:) / state%pdel(:ncol,:)

    end if   ! (do_exo_rad)

!!$    ! output rad inputs and resulting heating rates
!!$    call output_rad_data(pbuf, state, cam_in, landm, coszrs(i))

    ! Compute net radiative heating tendency
    call radheat_tend(state, pbuf,  ptend, qrl*cpair, qrs*cpair, fsns, &
         fsnt, flns, flnt, cam_in%asdir, net_flx)

    ! Compute heating rate for dtheta/dt
    do k=1,pver
       do i=1,ncol
          ftem(i,k) = (qrs(i,k) + qrl(i,k))/cpair * &
               (pstd/state%pmid(i,k))**cappa
       enddo
    enddo
    call outfld('HR      ',ftem    ,pcols   ,lchnk   )

    ! convert radiative heating rates to Q*dp for energy conservation
    if (conserve_energy) then
       !DIR$ CONCURRENT
       do k =1 , pver
          !DIR$ CONCURRENT
          do i = 1, ncol
             qrs(i,k) = qrs(i,k)*state%pdel(i,k)
             qrl(i,k) = qrl(i,k)*state%pdel(i,k)
          enddo
       enddo
    endif ! (conserve_energy)

    dtdt_solar(:ncol,:) = 0._r8
contains

  subroutine get_curr_calday_rotation(frac_day, day_in_year, offset)

    ! WOLF, subroutine get_curr_calday_rotation
    ! modified to calculate zenith angle by scaling factor, exo_ndays
    ! allows modification to diurnal cycles without mucking with the calendar
    ! and history file systems

    implicit none
    integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
    ! Positive for future times, negative
    ! for previous times.

    real(r8) ,intent(out) :: frac_day, day_in_year

    integer ::&
         yr,    &! year
         mon,   &! month
         day,   &! day of month
         tod     ! time of day (seconds past 0Z)
    integer, parameter :: y0 = 1950
    integer  :: leap_days, y
    real(r8) :: day_earth, f

    if (present(offset)) then
       call get_curr_date(yr, mon, day, tod, offset)
       day_earth = get_curr_calday(offset)
    else
       call get_curr_date(yr, mon, day, tod)
       day_earth = get_curr_calday()
    endif

    leap_days = 0
    !do y=y0, yr-1
    !   if ( mod(y,4) == 0 )    leap_days = leap_days+1
    !   if ( mod(y,100) == 0 )  leap_days = leap_days-1
    !   if ( mod(y,400) == 0 )  leap_days = leap_days+1
    !enddo

    day_earth = day_earth + 365.*(yr) + leap_days - 1
    frac_day = day_earth/exo_ndays - FLOOR(day_earth / exo_ndays)

    ! length of year scaling not currently available
    !day_in_year = day_earth - (365.)*FLOOR(day_earth/(365.))

    !write(*,*) "day_earth, frac_day", day_earth, frac_day

  end subroutine get_curr_calday_rotation


  subroutine zenith_rotation(frac_day,  calday  ,clat    , clon   ,coszrs  ,ncol    )
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Compute cosine of solar zenith angle for albedo and radiation
    !   computations. Includes scaling factor to alter length of day.
    !   Modified by Wolf, to scale dirunal cycle with exo_ndays
    !
    ! Method:
    ! <Describe the algorithm(s) used in the routine.>
    ! <Also include any applicable external references.>
    !
    ! Author: J. Kiehl
    ! Edited by Wolf, E.T. 2017.
    !-----------------------------------------------------------------------
    use shr_kind_mod, only: r8 => shr_kind_r8
    use shr_orb_mod
    use cam_control_mod, only: lambm0, obliqr, eccen, mvelpp
    implicit none

    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    integer, intent(in) :: ncol                 ! number of positions
    real(r8), intent(in) :: frac_day
    real(r8), intent(in) :: calday              ! Calendar day, including fraction
    real(r8), intent(in) :: clat(ncol)          ! Current centered latitude (radians)
    real(r8), intent(in) :: clon(ncol)          ! Centered longitude (radians)
    !
    ! Output arguments
    !
    real(r8), intent(out) :: coszrs(ncol)       ! Cosine solar zenith angle
    !
    !---------------------------Local variables-----------------------------
    !
    integer i         ! Position loop index
    real(r8) delta    ! Solar declination angle  in radians
    real(r8) eccf     ! Earth orbit eccentricity factor
    !
    !-----------------------------------------------------------------------
    !
    call shr_orb_decl (calday  ,eccen     ,mvelpp  ,lambm0  ,obliqr  , &
         delta   ,eccf      )
    !
    ! Compute local cosine solar zenith angle,
    !
    do i=1,ncol
       !coszrs(i) = shr_orb_cosz( calday, clat(i), clon(i), delta )
       coszrs(i) = shr_orb_cosz( frac_day, clat(i), clon(i), delta )
    end do

  end subroutine zenith_rotation

!!$  !-------------------------------------------------------------------------------
!!$
!!$  subroutine set_diags_fullsky()
!!$
!!$    ! Transform RRTMGP output for CAM and compute heating rates.
!!$    ! SW fluxes from RRTMGP are on daylight columns only, so expand to
!!$    ! full chunks for output to CAM history.
!!$
!!$    integer :: i
!!$    real(r8), dimension(size(fsw%bnd_flux_dn,1), &
!!$         size(fsw%bnd_flux_dn,2), &
!!$         size(fsw%bnd_flux_dn,3)) :: flux_dn_diffuse
!!$    !-------------------------------------------------------------------------
!!$
!!$    ! Initialize to provide 0.0 values for night columns.
!!$    fns               = 0._r8 ! net sw flux
!!$    fsds              = 0._r8 ! downward sw flux at surface
!!$    rd%fsutoa         = 0._r8 ! upward sw flux at TOA
!!$    rd%fsntoa         = 0._r8 ! net sw at TOA
!!$    rd%solin          = 0._r8 ! solar irradiance at TOA
!!$    rd%flux_sw_up     = 0._r8
!!$    rd%flux_sw_dn     = 0._r8
!!$
!!$    qrs      = 0._r8
!!$    fsns     = 0._r8
!!$    fsnt     = 0._r8
!!$
!!$    do i = 1, nday
!!$       fns(idxday(i),ktopcam:)  = fsw%flux_net(i, ktoprad:)
!!$       fsds(idxday(i))          = fsw%flux_dn(i, nlay+1)
!!$       rd%fsutoa(idxday(i))     = fsw%flux_up(i, 1)
!!$       rd%fsntoa(idxday(i))     = fsw%flux_net(i, 1)
!!$       rd%flux_sw_up(idxday(i),ktopcam:)     = fsw%flux_up(i,ktoprad:)
!!$       rd%flux_sw_dn(idxday(i),ktopcam:)     = fsw%flux_dn(i,ktoprad:)
!!$    end do
!!$
!!$!    ! Compute heating rate as a dry static energy tendency.
!!$!    call heating_rate('SW', ncol, fns, qrs)
!!$!    call heating_rate('SW', ncol, fcns, rd%qrsc)
!!$
!!$    fsns(:ncol)     = fns(:ncol,pverp)    ! net sw flux at surface
!!$    fsnt(:ncol)     = fns(:ncol,ktopcam)  ! net sw flux at top-of-model (w/o extra layer)
!!$
!!$    cam_out%netsw(:ncol) = fsns(:ncol)
!!$
!!$    if (spectralflux) then
!!$       su  = 0._r8
!!$       sd  = 0._r8
!!$       do i = 1, nday
!!$          su(idxday(i),ktopcam:,:) = fsw%bnd_flux_up(i,ktoprad:,:)
!!$          sd(idxday(i),ktopcam:,:) = fsw%bnd_flux_dn(i,ktoprad:,:)
!!$       end do
!!$    end if
!!$
!!$    ! Export surface fluxes
!!$    ! sols(pcols)      Direct solar rad on surface (< 0.7)
!!$    ! soll(pcols)      Direct solar rad on surface (>= 0.7)
!!$    ! RRTMGP: Near-IR bands (1-10), 820-16000 cm-1, 0.625-12.195 microns
!!$    ! Put half of band 10 in each of the UV/visible and near-IR values,
!!$    ! since this band straddles 0.7 microns:
!!$    ! UV/visible bands 10-13, 16000-50000 cm-1, 0.200-0.625 micron
!!$
!!$    ! reset fluxes
!!$    cam_out%sols  = 0.0_r8
!!$    cam_out%soll  = 0.0_r8
!!$    cam_out%solsd = 0.0_r8
!!$    cam_out%solld = 0.0_r8
!!$
!!$    ! Calculate diffuse flux from total and direct
!!$    flux_dn_diffuse = fsw%bnd_flux_dn - fsw%bnd_flux_dn_dir
!!$
!!$    do i = 1, nday
!!$       cam_out%soll(idxday(i)) = sum(fsw%bnd_flux_dn_dir(i,nlay+1,1:9))      &
!!$            + 0.5_r8 * fsw%bnd_flux_dn_dir(i,nlay+1,10)
!!$
!!$       cam_out%sols(idxday(i)) = 0.5_r8 * fsw%bnd_flux_dn_dir(i,nlay+1,10)   &
!!$            + sum(fsw%bnd_flux_dn_dir(i,nlay+1,11:14))
!!$
!!$       cam_out%solld(idxday(i)) = sum(flux_dn_diffuse(i,nlay+1,1:9))         &
!!$            + 0.5_r8 * flux_dn_diffuse(i,nlay+1,10)
!!$
!!$       cam_out%solsd(idxday(i)) = 0.5_r8 * flux_dn_diffuse(i, nlay+1, 10)    &
!!$            + sum(flux_dn_diffuse(i,nlay+1,11:14))
!!$    end do
!!$
!!$    ! Set CAM LW diagnostics
!!$    !----------------------------------------------------------------------------
!!$
!!$    fnl = 0._r8
!!$    fcnl = 0._r8
!!$
!!$    ! RTE-RRTMGP convention for net is (down - up) **CAM assumes (up - down) !!
!!$    fnl(:ncol,ktopcam:)  = -1._r8 * flw%flux_net(    :, ktoprad:)
!!$
!!$    rd%flux_lw_up(:ncol,ktopcam:)     = flw%flux_up( :, ktoprad:)
!!$    rd%flux_lw_dn(:ncol,ktopcam:)     = flw%flux_dn( :, ktoprad:)
!!$
!!$!    call heating_rate('LW', ncol, fnl, qrl)
!!$!    call heating_rate('LW', ncol, fcnl, rd%qrlc)
!!$
!!$    flns(:ncol) = fnl(:ncol, pverp)
!!$    flnt(:ncol) = fnl(:ncol, ktopcam)
!!$
!!$
!!$    cam_out%flwds(:ncol) = flw%flux_dn(:, nlay+1)
!!$
!!$    rd%flut(:ncol)  = flw%flux_up(:, ktoprad)
!!$
!!$    if (spectralflux) then
!!$       lu  = 0._r8
!!$       ld  = 0._r8
!!$       lu(:ncol, ktopcam:, :)  = flw%bnd_flux_up(:, ktoprad:, :)
!!$       ld(:ncol, ktopcam:, :)  = flw%bnd_flux_dn(:, ktoprad:, :)
!!$    end if
!!$
!!$  end subroutine set_diags_fullsky
!!$
!!$  subroutine set_diags_clearsky()
!!$
!!$    ! Transform RRTMGP output for CAM and compute heating rates.
!!$    ! SW fluxes from RRTMGP are on daylight columns only, so expand to
!!$    ! full chunks for output to CAM history.
!!$
!!$    integer :: i
!!$    real(r8), dimension(size(fsw%bnd_flux_dn,1), &
!!$         size(fsw%bnd_flux_dn,2), &
!!$         size(fsw%bnd_flux_dn,3)) :: flux_dn_diffuse
!!$    !-------------------------------------------------------------------------
!!$
!!$    ! Initialize to provide 0.0 values for night columns.
!!$    fcns              = 0._r8 ! net sw clearsky flux
!!$    rd%fsdsc          = 0._r8 ! downward sw clearsky flux at surface
!!$    rd%fsntoac        = 0._r8 ! net sw clearsky flux at TOA
!!$    rd%flux_sw_clr_up = 0._r8
!!$    rd%flux_sw_clr_dn = 0._r8
!!$
!!$    rd%qrsc  = 0._r8
!!$    rd%fsnsc = 0._r8
!!$    rd%fsntc = 0._r8
!!$
!!$    do i = 1, nday
!!$       fcns(idxday(i),ktopcam:) = fswc%flux_net(i,ktoprad:)
!!$       rd%fsdsc(idxday(i))      = fswc%flux_dn(i, nlay+1)
!!$       rd%fsntoac(idxday(i))    = fswc%flux_net(i, 1)
!!$       rd%solin(idxday(i))      = fswc%flux_dn(i, 1)
!!$       rd%flux_sw_clr_up(idxday(i),ktopcam:) = fswc%flux_up(i,ktoprad:)
!!$       rd%flux_sw_clr_dn(idxday(i),ktopcam:) = fswc%flux_dn(i,ktoprad:)
!!$    end do
!!$
!!$!    ! Compute heating rate as a dry static energy tendency.
!!$!    call heating_rate('SW', ncol, fns, qrs)
!!$!    call heating_rate('SW', ncol, fcns, rd%qrsc)
!!$
!!$    rd%fsnsc(:ncol) = fcns(:ncol,pverp)   ! net sw clearsky flux at surface
!!$    rd%fsntc(:ncol) = fcns(:ncol,ktopcam) ! net sw clearsky flux at top
!!$
!!$!jt    ! Output fluxes at 200 mb
!!$!jt    call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fns,  rd%fsn200)
!!$!jt    call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fcns, rd%fsn200c)
!!$!jt    if (hist_fld_active('FSNR')) then
!!$!jt       do i = 1,ncol
!!$!jt          call vertinterp(1, 1, pverp, state%pint(i,:), p_trop(i), fns(i,:), rd%fsnr(i))
!!$!jt       end do
!!$!jt    end if
!!$
!!$    if (spectralflux) then
!!$       su  = 0._r8
!!$       sd  = 0._r8
!!$       do i = 1, nday
!!$          su(idxday(i),ktopcam:,:) = fsw%bnd_flux_up(i,ktoprad:,:)
!!$          sd(idxday(i),ktopcam:,:) = fsw%bnd_flux_dn(i,ktoprad:,:)
!!$       end do
!!$    end if
!!$
!!$    ! Export surface fluxes
!!$    ! sols(pcols)      Direct solar rad on surface (< 0.7)
!!$    ! soll(pcols)      Direct solar rad on surface (>= 0.7)
!!$    ! RRTMGP: Near-IR bands (1-10), 820-16000 cm-1, 0.625-12.195 microns
!!$    ! Put half of band 10 in each of the UV/visible and near-IR values,
!!$    ! since this band straddles 0.7 microns:
!!$    ! UV/visible bands 10-13, 16000-50000 cm-1, 0.200-0.625 micron
!!$
!!$    ! reset fluxes
!!$    cam_out%sols  = 0.0_r8
!!$    cam_out%soll  = 0.0_r8
!!$    cam_out%solsd = 0.0_r8
!!$    cam_out%solld = 0.0_r8
!!$
!!$    ! Calculate diffuse flux from total and direct
!!$    flux_dn_diffuse = fsw%bnd_flux_dn - fsw%bnd_flux_dn_dir
!!$
!!$    do i = 1, nday
!!$       cam_out%soll(idxday(i)) = sum(fsw%bnd_flux_dn_dir(i,nlay+1,1:9))      &
!!$            + 0.5_r8 * fsw%bnd_flux_dn_dir(i,nlay+1,10)
!!$
!!$       cam_out%sols(idxday(i)) = 0.5_r8 * fsw%bnd_flux_dn_dir(i,nlay+1,10)   &
!!$            + sum(fsw%bnd_flux_dn_dir(i,nlay+1,11:14))
!!$
!!$       cam_out%solld(idxday(i)) = sum(flux_dn_diffuse(i,nlay+1,1:9))         &
!!$            + 0.5_r8 * flux_dn_diffuse(i,nlay+1,10)
!!$
!!$       cam_out%solsd(idxday(i)) = 0.5_r8 * flux_dn_diffuse(i, nlay+1, 10)    &
!!$            + sum(flux_dn_diffuse(i,nlay+1,11:14))
!!$    end do
!!$
!!$    ! Set CAM LW diagnostics
!!$    !----------------------------------------------------------------------------
!!$
!!$    fnl = 0._r8
!!$    fcnl = 0._r8
!!$
!!$    ! RTE-RRTMGP convention for net is (down - up) **CAM assumes (up - down) !!
!!$    fnl(:ncol,ktopcam:)  = -1._r8 * flw%flux_net(    :, ktoprad:)
!!$    fcnl(:ncol,ktopcam:) = -1._r8 * flwc%flux_net(   :, ktoprad:)
!!$
!!$    rd%flux_lw_up(:ncol,ktopcam:)     = flw%flux_up( :, ktoprad:)
!!$    rd%flux_lw_clr_up(:ncol,ktopcam:) = flwc%flux_up(:, ktoprad:)
!!$    rd%flux_lw_dn(:ncol,ktopcam:)     = flw%flux_dn( :, ktoprad:)
!!$    rd%flux_lw_clr_dn(:ncol,ktopcam:) = flwc%flux_dn(:, ktoprad:)
!!$
!!$    call heating_rate('LW', ncol, fnl, qrl)
!!$    call heating_rate('LW', ncol, fcnl, rd%qrlc)
!!$
!!$    flns(:ncol) = fnl(:ncol, pverp)
!!$    flnt(:ncol) = fnl(:ncol, ktopcam)
!!$
!!$    rd%flnsc(:ncol) = fcnl(:ncol, pverp)
!!$    rd%flntc(:ncol) = fcnl(:ncol, ktopcam)    ! net lw flux at top-of-model
!!$
!!$    cam_out%flwds(:ncol) = flw%flux_dn(:, nlay+1)
!!$    rd%fldsc(:ncol)      = flwc%flux_dn(:, nlay+1)
!!$
!!$    rd%flut(:ncol)  = flw%flux_up(:, ktoprad)
!!$    rd%flutc(:ncol) = flwc%flux_up(:, ktoprad)
!!$
!!$    ! Output fluxes at 200 mb
!!$!jt    call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fnl,  rd%fln200)
!!$!jt    call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fcnl, rd%fln200c)
!!$!jt    if (hist_fld_active('FLNR')) then
!!$!jt       do i = 1,ncol
!!$!jt          call vertinterp(1, 1, pverp, state%pint(i,:), p_trop(i), fnl(i,:), rd%flnr(i))
!!$!jt       end do
!!$!jt    end if
!!$
!!$    if (spectralflux) then
!!$       lu  = 0._r8
!!$       ld  = 0._r8
!!$       lu(:ncol, ktopcam:, :)  = flw%bnd_flux_up(:, ktoprad:, :)
!!$       ld(:ncol, ktopcam:, :)  = flw%bnd_flux_dn(:, ktoprad:, :)
!!$    end if
!!$
!!$  end subroutine set_diags_clearsky
!!$
!!$  !-------------------------------------------------------------------------------
!!$  subroutine radiation_output_fullsky(lchnk, ncol, icall, rd, pbuf, cam_out)
!!$
!!$    integer,                intent(in) :: lchnk
!!$    integer,                intent(in) :: ncol
!!$    integer,                intent(in) :: icall  ! icall=0 for climate diagnostics
!!$    type(rad_out_t),        intent(in) :: rd
!!$    type(physics_buffer_desc), pointer :: pbuf(:)
!!$    type(cam_out_t),        intent(in) :: cam_out
!!$
!!$    ! local variables
!!$    real(r8), pointer :: qrl(:,:)
!!$    real(r8), pointer :: flnt(:)
!!$    real(r8), pointer :: flns(:)
!!$    real(r8), pointer :: qrs(:,:)
!!$    real(r8), pointer :: fsnt(:)
!!$    real(r8), pointer :: fsns(:)
!!$    real(r8), pointer :: fsds(:)
!!$    real(r8) :: ftem(pcols)
!!$    !----------------------------------------------------------------------------
!!$
!!$    call pbuf_get_field(pbuf, qrs_idx,  qrs)
!!$    call pbuf_get_field(pbuf, fsnt_idx, fsnt)
!!$    call pbuf_get_field(pbuf, fsns_idx, fsns)
!!$    call pbuf_get_field(pbuf, fsds_idx, fsds)
!!$
!!$    call outfld('SOLIN'//diag(icall),    rd%solin,      pcols, lchnk)
!!$
!!$    ! QRS is output as temperature tendency.
!!$    call outfld('QRS'//diag(icall),      qrs(:ncol,:)/cpair,     ncol, lchnk)
!!$
!!$    call outfld('FSNT'//diag(icall),    fsnt,               pcols, lchnk)
!!$    call outfld('FSNTOA'//diag(icall),  rd%fsntoa,          pcols, lchnk)
!!$
!!$    ftem(:ncol) = rd%fsntoa(:ncol) - rd%fsntoac(:ncol)
!!$    call outfld('SWCF'//diag(icall),     ftem,          pcols, lchnk)
!!$
!!$    call outfld('FSUTOA'//diag(icall),   rd%fsutoa,     pcols, lchnk)
!!$
!!$    call outfld('FSNIRTOA'//diag(icall), rd%fsnirt,     pcols, lchnk)
!!$    call outfld('FSNRTOAS'//diag(icall), rd%fsnirtsq,   pcols, lchnk)
!!$
!!$    call outfld('FSN200'//diag(icall),   rd%fsn200,     pcols, lchnk)
!!$
!!$    call outfld('FSNR'//diag(icall),     rd%fsnr,       pcols, lchnk)
!!$
!!$    call outfld('SOLS'//diag(icall),     cam_out%sols,  pcols, lchnk)
!!$    call outfld('SOLL'//diag(icall),     cam_out%soll,  pcols, lchnk)
!!$    call outfld('SOLSD'//diag(icall),    cam_out%solsd, pcols, lchnk)
!!$    call outfld('SOLLD'//diag(icall),    cam_out%solld, pcols, lchnk)
!!$
!!$    call outfld('FSNS'//diag(icall),     fsns,          pcols, lchnk)
!!$
!!$    call outfld('FSDS'//diag(icall),     fsds,          pcols, lchnk)
!!$
!!$    call outfld('FUS'//diag(icall),  rd%flux_sw_up,     pcols, lchnk)
!!$    call outfld('FDS'//diag(icall),  rd%flux_sw_dn,     pcols, lchnk)
!!$
!!$    !============lw========
!!$    call pbuf_get_field(pbuf, qrl_idx,  qrl)
!!$    call pbuf_get_field(pbuf, flnt_idx, flnt)
!!$    call pbuf_get_field(pbuf, flns_idx, flns)
!!$
!!$    call outfld('QRL'//diag(icall),     qrl(:ncol,:)/cpair,     ncol, lchnk)
!!$
!!$    call outfld('FLNT'//diag(icall),    flnt,          pcols, lchnk)
!!$
!!$    call outfld('FLUT'//diag(icall),    rd%flut,       pcols, lchnk)
!!$
!!$    ftem(:ncol) = rd%flutc(:ncol) - rd%flut(:ncol)
!!$    call outfld('LWCF'//diag(icall),    ftem,          pcols, lchnk)
!!$
!!$    call outfld('FLN200'//diag(icall),  rd%fln200,     pcols, lchnk)
!!$
!!$    call outfld('FLNR'//diag(icall),    rd%flnr,       pcols, lchnk)
!!$
!!$    call outfld('FLNS'//diag(icall),    flns,          pcols, lchnk)
!!$
!!$    call outfld('FLDS'//diag(icall),    cam_out%flwds, pcols, lchnk)
!!$
!!$    call outfld('FDL'//diag(icall),  rd%flux_lw_dn,     pcols, lchnk)
!!$    call outfld('FUL'//diag(icall),  rd%flux_lw_up,     pcols, lchnk)
!!$
!!$  end subroutine radiation_output_fullsky
!!$
!!$  subroutine radiation_output_clearsky(lchnk, ncol, icall, rd, pbuf, cam_out)
!!$
!!$    integer,                intent(in) :: lchnk
!!$    integer,                intent(in) :: ncol
!!$    integer,                intent(in) :: icall  ! icall=0 for climate diagnostics
!!$    type(rad_out_t),        intent(in) :: rd
!!$    type(physics_buffer_desc), pointer :: pbuf(:)
!!$    type(cam_out_t),        intent(in) :: cam_out
!!$
!!$
!!$    ! local variables
!!$    real(r8), pointer :: qrs(:,:)
!!$    real(r8), pointer :: fsnt(:)
!!$    real(r8), pointer :: fsns(:)
!!$    real(r8), pointer :: fsds(:)
!!$    real(r8), pointer :: qrl(:,:)
!!$    real(r8), pointer :: flnt(:)
!!$    real(r8), pointer :: flns(:)
!!$    real(r8) :: ftem(pcols)
!!$    !----------------------------------------------------------------------------
!!$
!!$    ! QRS is output as temperature tendency.
!!$    call outfld('QRSC'//diag(icall),     rd%qrsc(:ncol,:)/cpair, ncol, lchnk)
!!$    call outfld('FSNTC'//diag(icall),   rd%fsntc,           pcols, lchnk)
!!$    call outfld('FSNTOAC'//diag(icall), rd%fsntoac,         pcols, lchnk)
!!$    call outfld('FSNRTOAC'//diag(icall), rd%fsnrtc,     pcols, lchnk)
!!$    call outfld('FSN200C'//diag(icall),  rd%fsn200c,    pcols, lchnk)
!!$    call outfld('FSNSC'//diag(icall),    rd%fsnsc,      pcols, lchnk)
!!$    call outfld('FSDSC'//diag(icall),    rd%fsdsc,      pcols, lchnk)
!!$    call outfld('FUSC'//diag(icall), rd%flux_sw_clr_up, pcols, lchnk)
!!$    call outfld('FDSC'//diag(icall), rd%flux_sw_clr_dn, pcols, lchnk)
!!$    !---LW-----------------------------------------------------------------------
!!$
!!$    call outfld('QRLC'//diag(icall),    rd%qrlc(:ncol,:)/cpair, ncol, lchnk)
!!$    call outfld('FLNTC'//diag(icall),   rd%flntc,      pcols, lchnk)
!!$    call outfld('FLUTC'//diag(icall),   rd%flutc,      pcols, lchnk)
!!$    !jt   call outfld('FLN200C'//diag(icall), rd%fln200c,    pcols, lchnk)
!!$    call outfld('FLNSC'//diag(icall),   rd%flnsc,      pcols, lchnk)
!!$    call outfld('FLDSC'//diag(icall),   rd%fldsc,      pcols, lchnk)
!!$    call outfld('FDLC'//diag(icall), rd%flux_lw_clr_dn, pcols, lchnk)
!!$    call outfld('FULC'//diag(icall), rd%flux_lw_clr_up, pcols, lchnk)
!!$
!!$  end subroutine radiation_output_clearsky

  subroutine outfld_spectral_flux_fullsky(lchnk, lwdown_rad_spec_in, lwup_rad_spec_in, swup_rad_spec_in, swdown_rad_spec_in )

!------------------------------------------------------------------------
!
! Purpose:  Add spectral output to master field list
!
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!
! Arguments
!
    integer,  intent(in) :: lchnk
    real(r8), intent(in), dimension(pcols,pverp,ntot_wavlnrng) ::  lwdown_rad_spec_in
    real(r8), intent(in), dimension(pcols,pverp,ntot_wavlnrng) ::  lwup_rad_spec_in
    real(r8), intent(in), dimension(pcols,pverp,ntot_wavlnrng) ::  swup_rad_spec_in
    real(r8), intent(in), dimension(pcols,pverp,ntot_wavlnrng) ::  swdown_rad_spec_in

!------------------------------------------------------------------------
!
! Start Code
!

    ! shortwave fields
    call outfld('FUS_int01     ',swup_rad_spec_in(:,:,1), pcols, lchnk)
    call outfld('FDS_int01     ',swdown_rad_spec_in(:,:,1), pcols, lchnk)
    call outfld('FUS_int02     ',swup_rad_spec_in(:,:,2), pcols, lchnk)
    call outfld('FDS_int02     ',swdown_rad_spec_in(:,:,2), pcols, lchnk)
    call outfld('FUS_int03     ',swup_rad_spec_in(:,:,3), pcols, lchnk)
    call outfld('FDS_int03     ',swdown_rad_spec_in(:,:,3), pcols, lchnk)
    call outfld('FUS_int04     ',swup_rad_spec_in(:,:,4), pcols, lchnk)
    call outfld('FDS_int04     ',swdown_rad_spec_in(:,:,4), pcols, lchnk)
    call outfld('FUS_int05     ',swup_rad_spec_in(:,:,5), pcols, lchnk)
    call outfld('FDS_int05     ',swdown_rad_spec_in(:,:,5), pcols, lchnk)
    call outfld('FUS_int06     ',swup_rad_spec_in(:,:,6), pcols, lchnk)
    call outfld('FDS_int06     ',swdown_rad_spec_in(:,:,6), pcols, lchnk)
    call outfld('FUS_int07     ',swup_rad_spec_in(:,:,7), pcols, lchnk)
    call outfld('FDS_int07     ',swdown_rad_spec_in(:,:,7), pcols, lchnk)
    call outfld('FUS_int08     ',swup_rad_spec_in(:,:,8), pcols, lchnk)
    call outfld('FDS_int08     ',swdown_rad_spec_in(:,:,8), pcols, lchnk)
    call outfld('FUS_int09     ',swup_rad_spec_in(:,:,9), pcols, lchnk)
    call outfld('FDS_int09     ',swdown_rad_spec_in(:,:,9), pcols, lchnk)
    call outfld('FUS_int10     ',swup_rad_spec_in(:,:,10), pcols, lchnk)
    call outfld('FDS_int10     ',swdown_rad_spec_in(:,:,10), pcols, lchnk)
    call outfld('FUS_int11     ',swup_rad_spec_in(:,:,11), pcols, lchnk)
    call outfld('FDS_int11     ',swdown_rad_spec_in(:,:,11), pcols, lchnk)
    call outfld('FUS_int12     ',swup_rad_spec_in(:,:,12), pcols, lchnk)
    call outfld('FDS_int12     ',swdown_rad_spec_in(:,:,12), pcols, lchnk)
    call outfld('FUS_int13     ',swup_rad_spec_in(:,:,13), pcols, lchnk)
    call outfld('FDS_int13     ',swdown_rad_spec_in(:,:,13), pcols, lchnk)
    call outfld('FUS_int14     ',swup_rad_spec_in(:,:,14), pcols, lchnk)
    call outfld('FDS_int14     ',swdown_rad_spec_in(:,:,14), pcols, lchnk)
    call outfld('FUS_int15     ',swup_rad_spec_in(:,:,15), pcols, lchnk)
    call outfld('FDS_int15     ',swdown_rad_spec_in(:,:,15), pcols, lchnk)
    call outfld('FUS_int16     ',swup_rad_spec_in(:,:,16), pcols, lchnk)
    call outfld('FDS_int16     ',swdown_rad_spec_in(:,:,16), pcols, lchnk)
    call outfld('FUS_int17     ',swup_rad_spec_in(:,:,17), pcols, lchnk)
    call outfld('FDS_int17     ',swdown_rad_spec_in(:,:,17), pcols, lchnk)
    call outfld('FUS_int18     ',swup_rad_spec_in(:,:,18), pcols, lchnk)
    call outfld('FDS_int18     ',swdown_rad_spec_in(:,:,18), pcols, lchnk)
    call outfld('FUS_int19     ',swup_rad_spec_in(:,:,19), pcols, lchnk)
    call outfld('FDS_int19     ',swdown_rad_spec_in(:,:,19), pcols, lchnk)
    call outfld('FUS_int20     ',swup_rad_spec_in(:,:,20), pcols, lchnk)
    call outfld('FDS_int20     ',swdown_rad_spec_in(:,:,20), pcols, lchnk)
    call outfld('FUS_int21     ',swup_rad_spec_in(:,:,21), pcols, lchnk)
    call outfld('FDS_int21     ',swdown_rad_spec_in(:,:,21), pcols, lchnk)
    call outfld('FUS_int22     ',swup_rad_spec_in(:,:,22), pcols, lchnk)
    call outfld('FDS_int22     ',swdown_rad_spec_in(:,:,22), pcols, lchnk)
    call outfld('FUS_int23     ',swup_rad_spec_in(:,:,23), pcols, lchnk)
    call outfld('FDS_int23     ',swdown_rad_spec_in(:,:,23), pcols, lchnk)
    call outfld('FUS_int24     ',swup_rad_spec_in(:,:,24), pcols, lchnk)
    call outfld('FDS_int24     ',swdown_rad_spec_in(:,:,24), pcols, lchnk)
    call outfld('FUS_int25     ',swup_rad_spec_in(:,:,25), pcols, lchnk)
    call outfld('FDS_int25     ',swdown_rad_spec_in(:,:,25), pcols, lchnk)
    call outfld('FUS_int26     ',swup_rad_spec_in(:,:,26), pcols, lchnk)
    call outfld('FDS_int26     ',swdown_rad_spec_in(:,:,26), pcols, lchnk)
    call outfld('FUS_int27     ',swup_rad_spec_in(:,:,27), pcols, lchnk)
    call outfld('FDS_int27     ',swdown_rad_spec_in(:,:,27), pcols, lchnk)
    call outfld('FUS_int28     ',swup_rad_spec_in(:,:,28), pcols, lchnk)
    call outfld('FDS_int28     ',swdown_rad_spec_in(:,:,28), pcols, lchnk)
    call outfld('FUS_int29     ',swup_rad_spec_in(:,:,29), pcols, lchnk)
    call outfld('FDS_int29     ',swdown_rad_spec_in(:,:,29), pcols, lchnk)
    call outfld('FUS_int30     ',swup_rad_spec_in(:,:,30), pcols, lchnk)
    call outfld('FDS_int30     ',swdown_rad_spec_in(:,:,30), pcols, lchnk)
    call outfld('FUS_int31     ',swup_rad_spec_in(:,:,31), pcols, lchnk)
    call outfld('FDS_int31     ',swdown_rad_spec_in(:,:,31), pcols, lchnk)
    call outfld('FUS_int32     ',swup_rad_spec_in(:,:,32), pcols, lchnk)
    call outfld('FDS_int32     ',swdown_rad_spec_in(:,:,32), pcols, lchnk)
    call outfld('FUS_int33     ',swup_rad_spec_in(:,:,33), pcols, lchnk)
    call outfld('FDS_int33     ',swdown_rad_spec_in(:,:,33), pcols, lchnk)
    call outfld('FUS_int34     ',swup_rad_spec_in(:,:,34), pcols, lchnk)
    call outfld('FDS_int34     ',swdown_rad_spec_in(:,:,34), pcols, lchnk)
    call outfld('FUS_int35     ',swup_rad_spec_in(:,:,35), pcols, lchnk)
    call outfld('FDS_int35     ',swdown_rad_spec_in(:,:,35), pcols, lchnk)
    call outfld('FUS_int36     ',swup_rad_spec_in(:,:,36), pcols, lchnk)
    call outfld('FDS_int36     ',swdown_rad_spec_in(:,:,36), pcols, lchnk)
    call outfld('FUS_int37     ',swup_rad_spec_in(:,:,37), pcols, lchnk)
    call outfld('FDS_int37     ',swdown_rad_spec_in(:,:,37), pcols, lchnk)
    call outfld('FUS_int38     ',swup_rad_spec_in(:,:,38), pcols, lchnk)
    call outfld('FDS_int38     ',swdown_rad_spec_in(:,:,38), pcols, lchnk)
    call outfld('FUS_int39     ',swup_rad_spec_in(:,:,39), pcols, lchnk)
    call outfld('FDS_int39     ',swdown_rad_spec_in(:,:,39), pcols, lchnk)
    call outfld('FUS_int40     ',swup_rad_spec_in(:,:,40), pcols, lchnk)
    call outfld('FDS_int40     ',swdown_rad_spec_in(:,:,40), pcols, lchnk)
    call outfld('FUS_int41     ',swup_rad_spec_in(:,:,41), pcols, lchnk)
    call outfld('FDS_int41     ',swdown_rad_spec_in(:,:,41), pcols, lchnk)
    call outfld('FUS_int42     ',swup_rad_spec_in(:,:,42), pcols, lchnk)
    call outfld('FDS_int42     ',swdown_rad_spec_in(:,:,42), pcols, lchnk)
    call outfld('FUS_int43     ',swup_rad_spec_in(:,:,43), pcols, lchnk)
    call outfld('FDS_int43     ',swdown_rad_spec_in(:,:,43), pcols, lchnk)
    call outfld('FUS_int44     ',swup_rad_spec_in(:,:,44), pcols, lchnk)
    call outfld('FDS_int44     ',swdown_rad_spec_in(:,:,44), pcols, lchnk)
    call outfld('FUS_int45     ',swup_rad_spec_in(:,:,45), pcols, lchnk)
    call outfld('FDS_int45     ',swdown_rad_spec_in(:,:,45), pcols, lchnk)
    call outfld('FUS_int46     ',swup_rad_spec_in(:,:,46), pcols, lchnk)
    call outfld('FDS_int46     ',swdown_rad_spec_in(:,:,46), pcols, lchnk)
    call outfld('FUS_int47     ',swup_rad_spec_in(:,:,47), pcols, lchnk)
    call outfld('FDS_int47     ',swdown_rad_spec_in(:,:,47), pcols, lchnk)
    call outfld('FUS_int48     ',swup_rad_spec_in(:,:,48), pcols, lchnk)
    call outfld('FDS_int48     ',swdown_rad_spec_in(:,:,48), pcols, lchnk)
    call outfld('FUS_int49     ',swup_rad_spec_in(:,:,49), pcols, lchnk)
    call outfld('FDS_int49     ',swdown_rad_spec_in(:,:,49), pcols, lchnk)
    call outfld('FUS_int50     ',swup_rad_spec_in(:,:,50), pcols, lchnk)
    call outfld('FDS_int50     ',swdown_rad_spec_in(:,:,50), pcols, lchnk)
    call outfld('FUS_int51     ',swup_rad_spec_in(:,:,51), pcols, lchnk)
    call outfld('FDS_int51     ',swdown_rad_spec_in(:,:,51), pcols, lchnk)
    call outfld('FUS_int52     ',swup_rad_spec_in(:,:,52), pcols, lchnk)
    call outfld('FDS_int52     ',swdown_rad_spec_in(:,:,52), pcols, lchnk)
    call outfld('FUS_int53     ',swup_rad_spec_in(:,:,53), pcols, lchnk)
    call outfld('FDS_int53     ',swdown_rad_spec_in(:,:,53), pcols, lchnk)
    call outfld('FUS_int54     ',swup_rad_spec_in(:,:,54), pcols, lchnk)
    call outfld('FDS_int54     ',swdown_rad_spec_in(:,:,54), pcols, lchnk)
    call outfld('FUS_int55     ',swup_rad_spec_in(:,:,55), pcols, lchnk)
    call outfld('FDS_int55     ',swdown_rad_spec_in(:,:,55), pcols, lchnk)
    call outfld('FUS_int56     ',swup_rad_spec_in(:,:,56), pcols, lchnk)
    call outfld('FDS_int56     ',swdown_rad_spec_in(:,:,56), pcols, lchnk)
    call outfld('FUS_int57     ',swup_rad_spec_in(:,:,57), pcols, lchnk)
    call outfld('FDS_int57     ',swdown_rad_spec_in(:,:,57), pcols, lchnk)
    call outfld('FUS_int58     ',swup_rad_spec_in(:,:,58), pcols, lchnk)
    call outfld('FDS_int58     ',swdown_rad_spec_in(:,:,58), pcols, lchnk)
    call outfld('FUS_int59     ',swup_rad_spec_in(:,:,59), pcols, lchnk)
    call outfld('FDS_int59     ',swdown_rad_spec_in(:,:,59), pcols, lchnk)
    call outfld('FUS_int60     ',swup_rad_spec_in(:,:,60), pcols, lchnk)
    call outfld('FDS_int60     ',swdown_rad_spec_in(:,:,60), pcols, lchnk)
    call outfld('FUS_int61     ',swup_rad_spec_in(:,:,61), pcols, lchnk)
    call outfld('FDS_int61     ',swdown_rad_spec_in(:,:,61), pcols, lchnk)
    call outfld('FUS_int62     ',swup_rad_spec_in(:,:,62), pcols, lchnk)
    call outfld('FDS_int62     ',swdown_rad_spec_in(:,:,62), pcols, lchnk)
    call outfld('FUS_int63     ',swup_rad_spec_in(:,:,63), pcols, lchnk)
    call outfld('FDS_int63     ',swdown_rad_spec_in(:,:,63), pcols, lchnk)
    call outfld('FUS_int64     ',swup_rad_spec_in(:,:,64), pcols, lchnk)
    call outfld('FDS_int64     ',swdown_rad_spec_in(:,:,64), pcols, lchnk)
    call outfld('FUS_int65     ',swup_rad_spec_in(:,:,65), pcols, lchnk)
    call outfld('FDS_int65     ',swdown_rad_spec_in(:,:,65), pcols, lchnk)
    call outfld('FUS_int66     ',swup_rad_spec_in(:,:,66), pcols, lchnk)
    call outfld('FDS_int66     ',swdown_rad_spec_in(:,:,66), pcols, lchnk)
    call outfld('FUS_int67     ',swup_rad_spec_in(:,:,67), pcols, lchnk)
    call outfld('FDS_int67     ',swdown_rad_spec_in(:,:,67), pcols, lchnk)
    call outfld('FUS_int68     ',swup_rad_spec_in(:,:,68), pcols, lchnk)
    call outfld('FDS_int68     ',swdown_rad_spec_in(:,:,68), pcols, lchnk)

    call outfld('FUS_toa_int01     ',swup_rad_spec_in(:,1,1), pcols, lchnk)
    call outfld('FDS_toa_int01     ',swdown_rad_spec_in(:,1,1), pcols, lchnk)
    call outfld('FUS_toa_int02     ',swup_rad_spec_in(:,1,2), pcols, lchnk)
    call outfld('FDS_toa_int02     ',swdown_rad_spec_in(:,1,2), pcols, lchnk)
    call outfld('FUS_toa_int03     ',swup_rad_spec_in(:,1,3), pcols, lchnk)
    call outfld('FDS_toa_int03     ',swdown_rad_spec_in(:,1,3), pcols, lchnk)
    call outfld('FUS_toa_int04     ',swup_rad_spec_in(:,1,4), pcols, lchnk)
    call outfld('FDS_toa_int04     ',swdown_rad_spec_in(:,1,4), pcols, lchnk)
    call outfld('FUS_toa_int05     ',swup_rad_spec_in(:,1,5), pcols, lchnk)
    call outfld('FDS_toa_int05     ',swdown_rad_spec_in(:,1,5), pcols, lchnk)
    call outfld('FUS_toa_int06     ',swup_rad_spec_in(:,1,6), pcols, lchnk)
    call outfld('FDS_toa_int06     ',swdown_rad_spec_in(:,1,6), pcols, lchnk)
    call outfld('FUS_toa_int07     ',swup_rad_spec_in(:,1,7), pcols, lchnk)
    call outfld('FDS_toa_int07     ',swdown_rad_spec_in(:,1,7), pcols, lchnk)
    call outfld('FUS_toa_int08     ',swup_rad_spec_in(:,1,8), pcols, lchnk)
    call outfld('FDS_toa_int08     ',swdown_rad_spec_in(:,1,8), pcols, lchnk)
    call outfld('FUS_toa_int09     ',swup_rad_spec_in(:,1,9), pcols, lchnk)
    call outfld('FDS_toa_int09     ',swdown_rad_spec_in(:,1,9), pcols, lchnk)
    call outfld('FUS_toa_int10     ',swup_rad_spec_in(:,1,10), pcols, lchnk)
    call outfld('FDS_toa_int10     ',swdown_rad_spec_in(:,1,10), pcols, lchnk)
    call outfld('FUS_toa_int11     ',swup_rad_spec_in(:,1,11), pcols, lchnk)
    call outfld('FDS_toa_int11     ',swdown_rad_spec_in(:,1,11), pcols, lchnk)
    call outfld('FUS_toa_int12     ',swup_rad_spec_in(:,1,12), pcols, lchnk)
    call outfld('FDS_toa_int12     ',swdown_rad_spec_in(:,1,12), pcols, lchnk)
    call outfld('FUS_toa_int13     ',swup_rad_spec_in(:,1,13), pcols, lchnk)
    call outfld('FDS_toa_int13     ',swdown_rad_spec_in(:,1,13), pcols, lchnk)
    call outfld('FUS_toa_int14     ',swup_rad_spec_in(:,1,14), pcols, lchnk)
    call outfld('FDS_toa_int14     ',swdown_rad_spec_in(:,1,14), pcols, lchnk)
    call outfld('FUS_toa_int15     ',swup_rad_spec_in(:,1,15), pcols, lchnk)
    call outfld('FDS_toa_int15     ',swdown_rad_spec_in(:,1,15), pcols, lchnk)
    call outfld('FUS_toa_int16     ',swup_rad_spec_in(:,1,16), pcols, lchnk)
    call outfld('FDS_toa_int16     ',swdown_rad_spec_in(:,1,16), pcols, lchnk)
    call outfld('FUS_toa_int17     ',swup_rad_spec_in(:,1,17), pcols, lchnk)
    call outfld('FDS_toa_int17     ',swdown_rad_spec_in(:,1,17), pcols, lchnk)
    call outfld('FUS_toa_int18     ',swup_rad_spec_in(:,1,18), pcols, lchnk)
    call outfld('FDS_toa_int18     ',swdown_rad_spec_in(:,1,18), pcols, lchnk)
    call outfld('FUS_toa_int19     ',swup_rad_spec_in(:,1,19), pcols, lchnk)
    call outfld('FDS_toa_int19     ',swdown_rad_spec_in(:,1,19), pcols, lchnk)
    call outfld('FUS_toa_int20     ',swup_rad_spec_in(:,1,20), pcols, lchnk)
    call outfld('FDS_toa_int20     ',swdown_rad_spec_in(:,1,20), pcols, lchnk)
    call outfld('FUS_toa_int21     ',swup_rad_spec_in(:,1,21), pcols, lchnk)
    call outfld('FDS_toa_int21     ',swdown_rad_spec_in(:,1,21), pcols, lchnk)
    call outfld('FUS_toa_int22     ',swup_rad_spec_in(:,1,22), pcols, lchnk)
    call outfld('FDS_toa_int22     ',swdown_rad_spec_in(:,1,22), pcols, lchnk)
    call outfld('FUS_toa_int23     ',swup_rad_spec_in(:,1,23), pcols, lchnk)
    call outfld('FDS_toa_int23     ',swdown_rad_spec_in(:,1,23), pcols, lchnk)
    call outfld('FUS_toa_int24     ',swup_rad_spec_in(:,1,24), pcols, lchnk)
    call outfld('FDS_toa_int24     ',swdown_rad_spec_in(:,1,24), pcols, lchnk)
    call outfld('FUS_toa_int25     ',swup_rad_spec_in(:,1,25), pcols, lchnk)
    call outfld('FDS_toa_int25     ',swdown_rad_spec_in(:,1,25), pcols, lchnk)
    call outfld('FUS_toa_int26     ',swup_rad_spec_in(:,1,26), pcols, lchnk)
    call outfld('FDS_toa_int26     ',swdown_rad_spec_in(:,1,26), pcols, lchnk)
    call outfld('FUS_toa_int27     ',swup_rad_spec_in(:,1,27), pcols, lchnk)
    call outfld('FDS_toa_int27     ',swdown_rad_spec_in(:,1,27), pcols, lchnk)
    call outfld('FUS_toa_int28     ',swup_rad_spec_in(:,1,28), pcols, lchnk)
    call outfld('FDS_toa_int28     ',swdown_rad_spec_in(:,1,28), pcols, lchnk)
    call outfld('FUS_toa_int29     ',swup_rad_spec_in(:,1,29), pcols, lchnk)
    call outfld('FDS_toa_int29     ',swdown_rad_spec_in(:,1,29), pcols, lchnk)
    call outfld('FUS_toa_int30     ',swup_rad_spec_in(:,1,30), pcols, lchnk)
    call outfld('FDS_toa_int30     ',swdown_rad_spec_in(:,1,30), pcols, lchnk)
    call outfld('FUS_toa_int31     ',swup_rad_spec_in(:,1,31), pcols, lchnk)
    call outfld('FDS_toa_int31     ',swdown_rad_spec_in(:,1,31), pcols, lchnk)
    call outfld('FUS_toa_int32     ',swup_rad_spec_in(:,1,32), pcols, lchnk)
    call outfld('FDS_toa_int32     ',swdown_rad_spec_in(:,1,32), pcols, lchnk)
    call outfld('FUS_toa_int33     ',swup_rad_spec_in(:,1,33), pcols, lchnk)
    call outfld('FDS_toa_int33     ',swdown_rad_spec_in(:,1,33), pcols, lchnk)
    call outfld('FUS_toa_int34     ',swup_rad_spec_in(:,1,34), pcols, lchnk)
    call outfld('FDS_toa_int34     ',swdown_rad_spec_in(:,1,34), pcols, lchnk)
    call outfld('FUS_toa_int35     ',swup_rad_spec_in(:,1,35), pcols, lchnk)
    call outfld('FDS_toa_int35     ',swdown_rad_spec_in(:,1,35), pcols, lchnk)
    call outfld('FUS_toa_int36     ',swup_rad_spec_in(:,1,36), pcols, lchnk)
    call outfld('FDS_toa_int36     ',swdown_rad_spec_in(:,1,36), pcols, lchnk)
    call outfld('FUS_toa_int37     ',swup_rad_spec_in(:,1,37), pcols, lchnk)
    call outfld('FDS_toa_int37     ',swdown_rad_spec_in(:,1,37), pcols, lchnk)
    call outfld('FUS_toa_int38     ',swup_rad_spec_in(:,1,38), pcols, lchnk)
    call outfld('FDS_toa_int38     ',swdown_rad_spec_in(:,1,38), pcols, lchnk)
    call outfld('FUS_toa_int39     ',swup_rad_spec_in(:,1,39), pcols, lchnk)
    call outfld('FDS_toa_int39     ',swdown_rad_spec_in(:,1,39), pcols, lchnk)
    call outfld('FUS_toa_int40     ',swup_rad_spec_in(:,1,40), pcols, lchnk)
    call outfld('FDS_toa_int40     ',swdown_rad_spec_in(:,1,40), pcols, lchnk)
    call outfld('FUS_toa_int41     ',swup_rad_spec_in(:,1,41), pcols, lchnk)
    call outfld('FDS_toa_int41     ',swdown_rad_spec_in(:,1,41), pcols, lchnk)
    call outfld('FUS_toa_int42     ',swup_rad_spec_in(:,1,42), pcols, lchnk)
    call outfld('FDS_toa_int42     ',swdown_rad_spec_in(:,1,42), pcols, lchnk)
    call outfld('FUS_toa_int43     ',swup_rad_spec_in(:,1,43), pcols, lchnk)
    call outfld('FDS_toa_int43     ',swdown_rad_spec_in(:,1,43), pcols, lchnk)
    call outfld('FUS_toa_int44     ',swup_rad_spec_in(:,1,44), pcols, lchnk)
    call outfld('FDS_toa_int44     ',swdown_rad_spec_in(:,1,44), pcols, lchnk)
    call outfld('FUS_toa_int45     ',swup_rad_spec_in(:,1,45), pcols, lchnk)
    call outfld('FDS_toa_int45     ',swdown_rad_spec_in(:,1,45), pcols, lchnk)
    call outfld('FUS_toa_int46     ',swup_rad_spec_in(:,1,46), pcols, lchnk)
    call outfld('FDS_toa_int46     ',swdown_rad_spec_in(:,1,46), pcols, lchnk)
    call outfld('FUS_toa_int47     ',swup_rad_spec_in(:,1,47), pcols, lchnk)
    call outfld('FDS_toa_int47     ',swdown_rad_spec_in(:,1,47), pcols, lchnk)
    call outfld('FUS_toa_int48     ',swup_rad_spec_in(:,1,48), pcols, lchnk)
    call outfld('FDS_toa_int48     ',swdown_rad_spec_in(:,1,48), pcols, lchnk)
    call outfld('FUS_toa_int49     ',swup_rad_spec_in(:,1,49), pcols, lchnk)
    call outfld('FDS_toa_int49     ',swdown_rad_spec_in(:,1,49), pcols, lchnk)
    call outfld('FUS_toa_int50     ',swup_rad_spec_in(:,1,50), pcols, lchnk)
    call outfld('FDS_toa_int50     ',swdown_rad_spec_in(:,1,50), pcols, lchnk)
    call outfld('FUS_toa_int51     ',swup_rad_spec_in(:,1,51), pcols, lchnk)
    call outfld('FDS_toa_int51     ',swdown_rad_spec_in(:,1,51), pcols, lchnk)
    call outfld('FUS_toa_int52     ',swup_rad_spec_in(:,1,52), pcols, lchnk)
    call outfld('FDS_toa_int52     ',swdown_rad_spec_in(:,1,52), pcols, lchnk)
    call outfld('FUS_toa_int53     ',swup_rad_spec_in(:,1,53), pcols, lchnk)
    call outfld('FDS_toa_int53     ',swdown_rad_spec_in(:,1,53), pcols, lchnk)
    call outfld('FUS_toa_int54     ',swup_rad_spec_in(:,1,54), pcols, lchnk)
    call outfld('FDS_toa_int54     ',swdown_rad_spec_in(:,1,54), pcols, lchnk)
    call outfld('FUS_toa_int55     ',swup_rad_spec_in(:,1,55), pcols, lchnk)
    call outfld('FDS_toa_int55     ',swdown_rad_spec_in(:,1,55), pcols, lchnk)
    call outfld('FUS_toa_int56     ',swup_rad_spec_in(:,1,56), pcols, lchnk)
    call outfld('FDS_toa_int56     ',swdown_rad_spec_in(:,1,56), pcols, lchnk)
    call outfld('FUS_toa_int57     ',swup_rad_spec_in(:,1,57), pcols, lchnk)
    call outfld('FDS_toa_int57     ',swdown_rad_spec_in(:,1,57), pcols, lchnk)
    call outfld('FUS_toa_int58     ',swup_rad_spec_in(:,1,58), pcols, lchnk)
    call outfld('FDS_toa_int58     ',swdown_rad_spec_in(:,1,58), pcols, lchnk)
    call outfld('FUS_toa_int59     ',swup_rad_spec_in(:,1,59), pcols, lchnk)
    call outfld('FDS_toa_int59     ',swdown_rad_spec_in(:,1,59), pcols, lchnk)
    call outfld('FUS_toa_int60     ',swup_rad_spec_in(:,1,60), pcols, lchnk)
    call outfld('FDS_toa_int60     ',swdown_rad_spec_in(:,1,60), pcols, lchnk)
    call outfld('FUS_toa_int61     ',swup_rad_spec_in(:,1,61), pcols, lchnk)
    call outfld('FDS_toa_int61     ',swdown_rad_spec_in(:,1,61), pcols, lchnk)
    call outfld('FUS_toa_int62     ',swup_rad_spec_in(:,1,62), pcols, lchnk)
    call outfld('FDS_toa_int62     ',swdown_rad_spec_in(:,1,62), pcols, lchnk)
    call outfld('FUS_toa_int63     ',swup_rad_spec_in(:,1,63), pcols, lchnk)
    call outfld('FDS_toa_int63     ',swdown_rad_spec_in(:,1,63), pcols, lchnk)
    call outfld('FUS_toa_int64     ',swup_rad_spec_in(:,1,64), pcols, lchnk)
    call outfld('FDS_toa_int64     ',swdown_rad_spec_in(:,1,64), pcols, lchnk)
    call outfld('FUS_toa_int65     ',swup_rad_spec_in(:,1,65), pcols, lchnk)
    call outfld('FDS_toa_int65     ',swdown_rad_spec_in(:,1,65), pcols, lchnk)
    call outfld('FUS_toa_int66     ',swup_rad_spec_in(:,1,66), pcols, lchnk)
    call outfld('FDS_toa_int66     ',swdown_rad_spec_in(:,1,66), pcols, lchnk)
    call outfld('FUS_toa_int67     ',swup_rad_spec_in(:,1,67), pcols, lchnk)
    call outfld('FDS_toa_int67     ',swdown_rad_spec_in(:,1,67), pcols, lchnk)
    call outfld('FUS_toa_int68     ',swup_rad_spec_in(:,1,68), pcols, lchnk)
    call outfld('FDS_toa_int68     ',swdown_rad_spec_in(:,1,68), pcols, lchnk)

    ! longwave
    call outfld('FUL_int01     ',lwup_rad_spec_in(:,:,1), pcols, lchnk)
    call outfld('FDL_int01     ',lwdown_rad_spec_in(:,:,1), pcols, lchnk)
    call outfld('FUL_int02     ',lwup_rad_spec_in(:,:,2), pcols, lchnk)
    call outfld('FDL_int02     ',lwdown_rad_spec_in(:,:,2), pcols, lchnk)
    call outfld('FUL_int03     ',lwup_rad_spec_in(:,:,3), pcols, lchnk)
    call outfld('FDL_int03     ',lwdown_rad_spec_in(:,:,3), pcols, lchnk)
    call outfld('FUL_int04     ',lwup_rad_spec_in(:,:,4), pcols, lchnk)
    call outfld('FDL_int04     ',lwdown_rad_spec_in(:,:,4), pcols, lchnk)
    call outfld('FUL_int05     ',lwup_rad_spec_in(:,:,5), pcols, lchnk)
    call outfld('FDL_int05     ',lwdown_rad_spec_in(:,:,5), pcols, lchnk)
    call outfld('FUL_int06     ',lwup_rad_spec_in(:,:,6), pcols, lchnk)
    call outfld('FDL_int06     ',lwdown_rad_spec_in(:,:,6), pcols, lchnk)
    call outfld('FUL_int07     ',lwup_rad_spec_in(:,:,7), pcols, lchnk)
    call outfld('FDL_int07     ',lwdown_rad_spec_in(:,:,7), pcols, lchnk)
    call outfld('FUL_int08     ',lwup_rad_spec_in(:,:,8), pcols, lchnk)
    call outfld('FDL_int08     ',lwdown_rad_spec_in(:,:,8), pcols, lchnk)
    call outfld('FUL_int09     ',lwup_rad_spec_in(:,:,9), pcols, lchnk)
    call outfld('FDL_int09     ',lwdown_rad_spec_in(:,:,9), pcols, lchnk)
    call outfld('FUL_int10     ',lwup_rad_spec_in(:,:,10), pcols, lchnk)
    call outfld('FDL_int10     ',lwdown_rad_spec_in(:,:,10), pcols, lchnk)
    call outfld('FUL_int11     ',lwup_rad_spec_in(:,:,11), pcols, lchnk)
    call outfld('FDL_int11     ',lwdown_rad_spec_in(:,:,11), pcols, lchnk)
    call outfld('FUL_int12     ',lwup_rad_spec_in(:,:,12), pcols, lchnk)
    call outfld('FDL_int12     ',lwdown_rad_spec_in(:,:,12), pcols, lchnk)
    call outfld('FUL_int13     ',lwup_rad_spec_in(:,:,13), pcols, lchnk)
    call outfld('FDL_int13     ',lwdown_rad_spec_in(:,:,13), pcols, lchnk)
    call outfld('FUL_int14     ',lwup_rad_spec_in(:,:,14), pcols, lchnk)
    call outfld('FDL_int14     ',lwdown_rad_spec_in(:,:,14), pcols, lchnk)
    call outfld('FUL_int15     ',lwup_rad_spec_in(:,:,15), pcols, lchnk)
    call outfld('FDL_int15     ',lwdown_rad_spec_in(:,:,15), pcols, lchnk)
    call outfld('FUL_int16     ',lwup_rad_spec_in(:,:,16), pcols, lchnk)
    call outfld('FDL_int16     ',lwdown_rad_spec_in(:,:,16), pcols, lchnk)
    call outfld('FUL_int17     ',lwup_rad_spec_in(:,:,17), pcols, lchnk)
    call outfld('FDL_int17     ',lwdown_rad_spec_in(:,:,17), pcols, lchnk)
    call outfld('FUL_int18     ',lwup_rad_spec_in(:,:,18), pcols, lchnk)
    call outfld('FDL_int18     ',lwdown_rad_spec_in(:,:,18), pcols, lchnk)
    call outfld('FUL_int19     ',lwup_rad_spec_in(:,:,19), pcols, lchnk)
    call outfld('FDL_int19     ',lwdown_rad_spec_in(:,:,19), pcols, lchnk)
    call outfld('FUL_int20     ',lwup_rad_spec_in(:,:,20), pcols, lchnk)
    call outfld('FDL_int20     ',lwdown_rad_spec_in(:,:,20), pcols, lchnk)
    call outfld('FUL_int21     ',lwup_rad_spec_in(:,:,21), pcols, lchnk)
    call outfld('FDL_int21     ',lwdown_rad_spec_in(:,:,21), pcols, lchnk)
    call outfld('FUL_int22     ',lwup_rad_spec_in(:,:,22), pcols, lchnk)
    call outfld('FDL_int22     ',lwdown_rad_spec_in(:,:,22), pcols, lchnk)
    call outfld('FUL_int23     ',lwup_rad_spec_in(:,:,23), pcols, lchnk)
    call outfld('FDL_int23     ',lwdown_rad_spec_in(:,:,23), pcols, lchnk)
    call outfld('FUL_int24     ',lwup_rad_spec_in(:,:,24), pcols, lchnk)
    call outfld('FDL_int24     ',lwdown_rad_spec_in(:,:,24), pcols, lchnk)
    call outfld('FUL_int25     ',lwup_rad_spec_in(:,:,25), pcols, lchnk)
    call outfld('FDL_int25     ',lwdown_rad_spec_in(:,:,25), pcols, lchnk)
    call outfld('FUL_int26     ',lwup_rad_spec_in(:,:,26), pcols, lchnk)
    call outfld('FDL_int26     ',lwdown_rad_spec_in(:,:,26), pcols, lchnk)
    call outfld('FUL_int27     ',lwup_rad_spec_in(:,:,27), pcols, lchnk)
    call outfld('FDL_int27     ',lwdown_rad_spec_in(:,:,27), pcols, lchnk)
    call outfld('FUL_int28     ',lwup_rad_spec_in(:,:,28), pcols, lchnk)
    call outfld('FDL_int28     ',lwdown_rad_spec_in(:,:,28), pcols, lchnk)
    call outfld('FUL_int29     ',lwup_rad_spec_in(:,:,29), pcols, lchnk)
    call outfld('FDL_int29     ',lwdown_rad_spec_in(:,:,29), pcols, lchnk)
    call outfld('FUL_int30     ',lwup_rad_spec_in(:,:,30), pcols, lchnk)
    call outfld('FDL_int30     ',lwdown_rad_spec_in(:,:,30), pcols, lchnk)
    call outfld('FUL_int31     ',lwup_rad_spec_in(:,:,31), pcols, lchnk)
    call outfld('FDL_int31     ',lwdown_rad_spec_in(:,:,31), pcols, lchnk)
    call outfld('FUL_int32     ',lwup_rad_spec_in(:,:,32), pcols, lchnk)
    call outfld('FDL_int32     ',lwdown_rad_spec_in(:,:,32), pcols, lchnk)
    call outfld('FUL_int33     ',lwup_rad_spec_in(:,:,33), pcols, lchnk)
    call outfld('FDL_int33     ',lwdown_rad_spec_in(:,:,33), pcols, lchnk)
    call outfld('FUL_int34     ',lwup_rad_spec_in(:,:,34), pcols, lchnk)
    call outfld('FDL_int34     ',lwdown_rad_spec_in(:,:,34), pcols, lchnk)
    call outfld('FUL_int35     ',lwup_rad_spec_in(:,:,35), pcols, lchnk)
    call outfld('FDL_int35     ',lwdown_rad_spec_in(:,:,35), pcols, lchnk)
    call outfld('FUL_int36     ',lwup_rad_spec_in(:,:,36), pcols, lchnk)
    call outfld('FDL_int36     ',lwdown_rad_spec_in(:,:,36), pcols, lchnk)
    call outfld('FUL_int37     ',lwup_rad_spec_in(:,:,37), pcols, lchnk)
    call outfld('FDL_int37     ',lwdown_rad_spec_in(:,:,37), pcols, lchnk)
    call outfld('FUL_int38     ',lwup_rad_spec_in(:,:,38), pcols, lchnk)
    call outfld('FDL_int38     ',lwdown_rad_spec_in(:,:,38), pcols, lchnk)
    call outfld('FUL_int39     ',lwup_rad_spec_in(:,:,39), pcols, lchnk)
    call outfld('FDL_int39     ',lwdown_rad_spec_in(:,:,39), pcols, lchnk)
    call outfld('FUL_int40     ',lwup_rad_spec_in(:,:,40), pcols, lchnk)
    call outfld('FDL_int40     ',lwdown_rad_spec_in(:,:,40), pcols, lchnk)
    call outfld('FUL_int41     ',lwup_rad_spec_in(:,:,41), pcols, lchnk)
    call outfld('FDL_int41     ',lwdown_rad_spec_in(:,:,41), pcols, lchnk)
    call outfld('FUL_int42     ',lwup_rad_spec_in(:,:,42), pcols, lchnk)
    call outfld('FDL_int42     ',lwdown_rad_spec_in(:,:,42), pcols, lchnk)
    call outfld('FUL_int43     ',lwup_rad_spec_in(:,:,43), pcols, lchnk)
    call outfld('FDL_int43     ',lwdown_rad_spec_in(:,:,43), pcols, lchnk)
    call outfld('FUL_int44     ',lwup_rad_spec_in(:,:,44), pcols, lchnk)
    call outfld('FDL_int44     ',lwdown_rad_spec_in(:,:,44), pcols, lchnk)
    call outfld('FUL_int45     ',lwup_rad_spec_in(:,:,45), pcols, lchnk)
    call outfld('FDL_int45     ',lwdown_rad_spec_in(:,:,45), pcols, lchnk)
    call outfld('FUL_int46     ',lwup_rad_spec_in(:,:,46), pcols, lchnk)
    call outfld('FDL_int46     ',lwdown_rad_spec_in(:,:,46), pcols, lchnk)
    call outfld('FUL_int47     ',lwup_rad_spec_in(:,:,47), pcols, lchnk)
    call outfld('FDL_int47     ',lwdown_rad_spec_in(:,:,47), pcols, lchnk)
    call outfld('FUL_int48     ',lwup_rad_spec_in(:,:,48), pcols, lchnk)
    call outfld('FDL_int48     ',lwdown_rad_spec_in(:,:,48), pcols, lchnk)
    call outfld('FUL_int49     ',lwup_rad_spec_in(:,:,49), pcols, lchnk)
    call outfld('FDL_int49     ',lwdown_rad_spec_in(:,:,49), pcols, lchnk)
    call outfld('FUL_int50     ',lwup_rad_spec_in(:,:,50), pcols, lchnk)
    call outfld('FDL_int50     ',lwdown_rad_spec_in(:,:,50), pcols, lchnk)
    call outfld('FUL_int51     ',lwup_rad_spec_in(:,:,51), pcols, lchnk)
    call outfld('FDL_int51     ',lwdown_rad_spec_in(:,:,51), pcols, lchnk)
    call outfld('FUL_int52     ',lwup_rad_spec_in(:,:,52), pcols, lchnk)
    call outfld('FDL_int52     ',lwdown_rad_spec_in(:,:,52), pcols, lchnk)
    call outfld('FUL_int53     ',lwup_rad_spec_in(:,:,53), pcols, lchnk)
    call outfld('FDL_int53     ',lwdown_rad_spec_in(:,:,53), pcols, lchnk)
    call outfld('FUL_int54     ',lwup_rad_spec_in(:,:,54), pcols, lchnk)
    call outfld('FDL_int54     ',lwdown_rad_spec_in(:,:,54), pcols, lchnk)
    call outfld('FUL_int55     ',lwup_rad_spec_in(:,:,55), pcols, lchnk)
    call outfld('FDL_int55     ',lwdown_rad_spec_in(:,:,55), pcols, lchnk)
    call outfld('FUL_int56     ',lwup_rad_spec_in(:,:,56), pcols, lchnk)
    call outfld('FDL_int56     ',lwdown_rad_spec_in(:,:,56), pcols, lchnk)
    call outfld('FUL_int57     ',lwup_rad_spec_in(:,:,57), pcols, lchnk)
    call outfld('FDL_int57     ',lwdown_rad_spec_in(:,:,57), pcols, lchnk)
    call outfld('FUL_int58     ',lwup_rad_spec_in(:,:,58), pcols, lchnk)
    call outfld('FDL_int58     ',lwdown_rad_spec_in(:,:,58), pcols, lchnk)
    call outfld('FUL_int59     ',lwup_rad_spec_in(:,:,59), pcols, lchnk)
    call outfld('FDL_int59     ',lwdown_rad_spec_in(:,:,59), pcols, lchnk)
    call outfld('FUL_int60     ',lwup_rad_spec_in(:,:,60), pcols, lchnk)
    call outfld('FDL_int60     ',lwdown_rad_spec_in(:,:,60), pcols, lchnk)
    call outfld('FUL_int61     ',lwup_rad_spec_in(:,:,61), pcols, lchnk)
    call outfld('FDL_int61     ',lwdown_rad_spec_in(:,:,61), pcols, lchnk)
    call outfld('FUL_int62     ',lwup_rad_spec_in(:,:,62), pcols, lchnk)
    call outfld('FDL_int62     ',lwdown_rad_spec_in(:,:,62), pcols, lchnk)
    call outfld('FUL_int63     ',lwup_rad_spec_in(:,:,63), pcols, lchnk)
    call outfld('FDL_int63     ',lwdown_rad_spec_in(:,:,63), pcols, lchnk)
    call outfld('FUL_int64     ',lwup_rad_spec_in(:,:,64), pcols, lchnk)
    call outfld('FDL_int64     ',lwdown_rad_spec_in(:,:,64), pcols, lchnk)
    call outfld('FUL_int65     ',lwup_rad_spec_in(:,:,65), pcols, lchnk)
    call outfld('FDL_int65     ',lwdown_rad_spec_in(:,:,65), pcols, lchnk)
    call outfld('FUL_int66     ',lwup_rad_spec_in(:,:,66), pcols, lchnk)
    call outfld('FDL_int66     ',lwdown_rad_spec_in(:,:,66), pcols, lchnk)
    call outfld('FUL_int67     ',lwup_rad_spec_in(:,:,67), pcols, lchnk)
    call outfld('FDL_int67     ',lwdown_rad_spec_in(:,:,67), pcols, lchnk)
    call outfld('FUL_int68     ',lwup_rad_spec_in(:,:,68), pcols, lchnk)
    call outfld('FDL_int68     ',lwdown_rad_spec_in(:,:,68), pcols, lchnk)

    call outfld('FUL_toa_int01     ',lwup_rad_spec_in(:,1,1), pcols, lchnk)
    call outfld('FDL_toa_int01     ',lwdown_rad_spec_in(:,1,1), pcols, lchnk)
    call outfld('FUL_toa_int02     ',lwup_rad_spec_in(:,1,2), pcols, lchnk)
    call outfld('FDL_toa_int02     ',lwdown_rad_spec_in(:,1,2), pcols, lchnk)
    call outfld('FUL_toa_int03     ',lwup_rad_spec_in(:,1,3), pcols, lchnk)
    call outfld('FDL_toa_int03     ',lwdown_rad_spec_in(:,1,3), pcols, lchnk)
    call outfld('FUL_toa_int04     ',lwup_rad_spec_in(:,1,4), pcols, lchnk)
    call outfld('FDL_toa_int04     ',lwdown_rad_spec_in(:,1,4), pcols, lchnk)
    call outfld('FUL_toa_int05     ',lwup_rad_spec_in(:,1,5), pcols, lchnk)
    call outfld('FDL_toa_int05     ',lwdown_rad_spec_in(:,1,5), pcols, lchnk)
    call outfld('FUL_toa_int06     ',lwup_rad_spec_in(:,1,6), pcols, lchnk)
    call outfld('FDL_toa_int06     ',lwdown_rad_spec_in(:,1,6), pcols, lchnk)
    call outfld('FUL_toa_int07     ',lwup_rad_spec_in(:,1,7), pcols, lchnk)
    call outfld('FDL_toa_int07     ',lwdown_rad_spec_in(:,1,7), pcols, lchnk)
    call outfld('FUL_toa_int08     ',lwup_rad_spec_in(:,1,8), pcols, lchnk)
    call outfld('FDL_toa_int08     ',lwdown_rad_spec_in(:,1,8), pcols, lchnk)
    call outfld('FUL_toa_int09     ',lwup_rad_spec_in(:,1,9), pcols, lchnk)
    call outfld('FDL_toa_int09     ',lwdown_rad_spec_in(:,1,9), pcols, lchnk)
    call outfld('FUL_toa_int10     ',lwup_rad_spec_in(:,1,10), pcols, lchnk)
    call outfld('FDL_toa_int10     ',lwdown_rad_spec_in(:,1,10), pcols, lchnk)
    call outfld('FUL_toa_int11     ',lwup_rad_spec_in(:,1,11), pcols, lchnk)
    call outfld('FDL_toa_int11     ',lwdown_rad_spec_in(:,1,11), pcols, lchnk)
    call outfld('FUL_toa_int12     ',lwup_rad_spec_in(:,1,12), pcols, lchnk)
    call outfld('FDL_toa_int12     ',lwdown_rad_spec_in(:,1,12), pcols, lchnk)
    call outfld('FUL_toa_int13     ',lwup_rad_spec_in(:,1,13), pcols, lchnk)
    call outfld('FDL_toa_int13     ',lwdown_rad_spec_in(:,1,13), pcols, lchnk)
    call outfld('FUL_toa_int14     ',lwup_rad_spec_in(:,1,14), pcols, lchnk)
    call outfld('FDL_toa_int14     ',lwdown_rad_spec_in(:,1,14), pcols, lchnk)
    call outfld('FUL_toa_int15     ',lwup_rad_spec_in(:,1,15), pcols, lchnk)
    call outfld('FDL_toa_int15     ',lwdown_rad_spec_in(:,1,15), pcols, lchnk)
    call outfld('FUL_toa_int16     ',lwup_rad_spec_in(:,1,16), pcols, lchnk)
    call outfld('FDL_toa_int16     ',lwdown_rad_spec_in(:,1,16), pcols, lchnk)
    call outfld('FUL_toa_int17     ',lwup_rad_spec_in(:,1,17), pcols, lchnk)
    call outfld('FDL_toa_int17     ',lwdown_rad_spec_in(:,1,17), pcols, lchnk)
    call outfld('FUL_toa_int18     ',lwup_rad_spec_in(:,1,18), pcols, lchnk)
    call outfld('FDL_toa_int18     ',lwdown_rad_spec_in(:,1,18), pcols, lchnk)
    call outfld('FUL_toa_int19     ',lwup_rad_spec_in(:,1,19), pcols, lchnk)
    call outfld('FDL_toa_int19     ',lwdown_rad_spec_in(:,1,19), pcols, lchnk)
    call outfld('FUL_toa_int20     ',lwup_rad_spec_in(:,1,20), pcols, lchnk)
    call outfld('FDL_toa_int20     ',lwdown_rad_spec_in(:,1,20), pcols, lchnk)
    call outfld('FUL_toa_int21     ',lwup_rad_spec_in(:,1,21), pcols, lchnk)
    call outfld('FDL_toa_int21     ',lwdown_rad_spec_in(:,1,21), pcols, lchnk)
    call outfld('FUL_toa_int22     ',lwup_rad_spec_in(:,1,22), pcols, lchnk)
    call outfld('FDL_toa_int22     ',lwdown_rad_spec_in(:,1,22), pcols, lchnk)
    call outfld('FUL_toa_int23     ',lwup_rad_spec_in(:,1,23), pcols, lchnk)
    call outfld('FDL_toa_int23     ',lwdown_rad_spec_in(:,1,23), pcols, lchnk)
    call outfld('FUL_toa_int24     ',lwup_rad_spec_in(:,1,24), pcols, lchnk)
    call outfld('FDL_toa_int24     ',lwdown_rad_spec_in(:,1,24), pcols, lchnk)
    call outfld('FUL_toa_int25     ',lwup_rad_spec_in(:,1,25), pcols, lchnk)
    call outfld('FDL_toa_int25     ',lwdown_rad_spec_in(:,1,25), pcols, lchnk)
    call outfld('FUL_toa_int26     ',lwup_rad_spec_in(:,1,26), pcols, lchnk)
    call outfld('FDL_toa_int26     ',lwdown_rad_spec_in(:,1,26), pcols, lchnk)
    call outfld('FUL_toa_int27     ',lwup_rad_spec_in(:,1,27), pcols, lchnk)
    call outfld('FDL_toa_int27     ',lwdown_rad_spec_in(:,1,27), pcols, lchnk)
    call outfld('FUL_toa_int28     ',lwup_rad_spec_in(:,1,28), pcols, lchnk)
    call outfld('FDL_toa_int28     ',lwdown_rad_spec_in(:,1,28), pcols, lchnk)
    call outfld('FUL_toa_int29     ',lwup_rad_spec_in(:,1,29), pcols, lchnk)
    call outfld('FDL_toa_int29     ',lwdown_rad_spec_in(:,1,29), pcols, lchnk)
    call outfld('FUL_toa_int30     ',lwup_rad_spec_in(:,1,30), pcols, lchnk)
    call outfld('FDL_toa_int30     ',lwdown_rad_spec_in(:,1,30), pcols, lchnk)
    call outfld('FUL_toa_int31     ',lwup_rad_spec_in(:,1,31), pcols, lchnk)
    call outfld('FDL_toa_int31     ',lwdown_rad_spec_in(:,1,31), pcols, lchnk)
    call outfld('FUL_toa_int32     ',lwup_rad_spec_in(:,1,32), pcols, lchnk)
    call outfld('FDL_toa_int32     ',lwdown_rad_spec_in(:,1,32), pcols, lchnk)
    call outfld('FUL_toa_int33     ',lwup_rad_spec_in(:,1,33), pcols, lchnk)
    call outfld('FDL_toa_int33     ',lwdown_rad_spec_in(:,1,33), pcols, lchnk)
    call outfld('FUL_toa_int34     ',lwup_rad_spec_in(:,1,34), pcols, lchnk)
    call outfld('FDL_toa_int34     ',lwdown_rad_spec_in(:,1,34), pcols, lchnk)
    call outfld('FUL_toa_int35     ',lwup_rad_spec_in(:,1,35), pcols, lchnk)
    call outfld('FDL_toa_int35     ',lwdown_rad_spec_in(:,1,35), pcols, lchnk)
    call outfld('FUL_toa_int36     ',lwup_rad_spec_in(:,1,36), pcols, lchnk)
    call outfld('FDL_toa_int36     ',lwdown_rad_spec_in(:,1,36), pcols, lchnk)
    call outfld('FUL_toa_int37     ',lwup_rad_spec_in(:,1,37), pcols, lchnk)
    call outfld('FDL_toa_int37     ',lwdown_rad_spec_in(:,1,37), pcols, lchnk)
    call outfld('FUL_toa_int38     ',lwup_rad_spec_in(:,1,38), pcols, lchnk)
    call outfld('FDL_toa_int38     ',lwdown_rad_spec_in(:,1,38), pcols, lchnk)
    call outfld('FUL_toa_int39     ',lwup_rad_spec_in(:,1,39), pcols, lchnk)
    call outfld('FDL_toa_int39     ',lwdown_rad_spec_in(:,1,39), pcols, lchnk)
    call outfld('FUL_toa_int40     ',lwup_rad_spec_in(:,1,40), pcols, lchnk)
    call outfld('FDL_toa_int40     ',lwdown_rad_spec_in(:,1,40), pcols, lchnk)
    call outfld('FUL_toa_int41     ',lwup_rad_spec_in(:,1,41), pcols, lchnk)
    call outfld('FDL_toa_int41     ',lwdown_rad_spec_in(:,1,41), pcols, lchnk)
    call outfld('FUL_toa_int42     ',lwup_rad_spec_in(:,1,42), pcols, lchnk)
    call outfld('FDL_toa_int42     ',lwdown_rad_spec_in(:,1,42), pcols, lchnk)
    call outfld('FUL_toa_int43     ',lwup_rad_spec_in(:,1,43), pcols, lchnk)
    call outfld('FDL_toa_int43     ',lwdown_rad_spec_in(:,1,43), pcols, lchnk)
    call outfld('FUL_toa_int44     ',lwup_rad_spec_in(:,1,44), pcols, lchnk)
    call outfld('FDL_toa_int44     ',lwdown_rad_spec_in(:,1,44), pcols, lchnk)
    call outfld('FUL_toa_int45     ',lwup_rad_spec_in(:,1,45), pcols, lchnk)
    call outfld('FDL_toa_int45     ',lwdown_rad_spec_in(:,1,45), pcols, lchnk)
    call outfld('FUL_toa_int46     ',lwup_rad_spec_in(:,1,46), pcols, lchnk)
    call outfld('FDL_toa_int46     ',lwdown_rad_spec_in(:,1,46), pcols, lchnk)
    call outfld('FUL_toa_int47     ',lwup_rad_spec_in(:,1,47), pcols, lchnk)
    call outfld('FDL_toa_int47     ',lwdown_rad_spec_in(:,1,47), pcols, lchnk)
    call outfld('FUL_toa_int48     ',lwup_rad_spec_in(:,1,48), pcols, lchnk)
    call outfld('FDL_toa_int48     ',lwdown_rad_spec_in(:,1,48), pcols, lchnk)
    call outfld('FUL_toa_int49     ',lwup_rad_spec_in(:,1,49), pcols, lchnk)
    call outfld('FDL_toa_int49     ',lwdown_rad_spec_in(:,1,49), pcols, lchnk)
    call outfld('FUL_toa_int50     ',lwup_rad_spec_in(:,1,50), pcols, lchnk)
    call outfld('FDL_toa_int50     ',lwdown_rad_spec_in(:,1,50), pcols, lchnk)
    call outfld('FUL_toa_int51     ',lwup_rad_spec_in(:,1,51), pcols, lchnk)
    call outfld('FDL_toa_int51     ',lwdown_rad_spec_in(:,1,51), pcols, lchnk)
    call outfld('FUL_toa_int52     ',lwup_rad_spec_in(:,1,52), pcols, lchnk)
    call outfld('FDL_toa_int52     ',lwdown_rad_spec_in(:,1,52), pcols, lchnk)
    call outfld('FUL_toa_int53     ',lwup_rad_spec_in(:,1,53), pcols, lchnk)
    call outfld('FDL_toa_int53     ',lwdown_rad_spec_in(:,1,53), pcols, lchnk)
    call outfld('FUL_toa_int54     ',lwup_rad_spec_in(:,1,54), pcols, lchnk)
    call outfld('FDL_toa_int54     ',lwdown_rad_spec_in(:,1,54), pcols, lchnk)
    call outfld('FUL_toa_int55     ',lwup_rad_spec_in(:,1,55), pcols, lchnk)
    call outfld('FDL_toa_int55     ',lwdown_rad_spec_in(:,1,55), pcols, lchnk)
    call outfld('FUL_toa_int56     ',lwup_rad_spec_in(:,1,56), pcols, lchnk)
    call outfld('FDL_toa_int56     ',lwdown_rad_spec_in(:,1,56), pcols, lchnk)
    call outfld('FUL_toa_int57     ',lwup_rad_spec_in(:,1,57), pcols, lchnk)
    call outfld('FDL_toa_int57     ',lwdown_rad_spec_in(:,1,57), pcols, lchnk)
    call outfld('FUL_toa_int58     ',lwup_rad_spec_in(:,1,58), pcols, lchnk)
    call outfld('FDL_toa_int58     ',lwdown_rad_spec_in(:,1,58), pcols, lchnk)
    call outfld('FUL_toa_int59     ',lwup_rad_spec_in(:,1,59), pcols, lchnk)
    call outfld('FDL_toa_int59     ',lwdown_rad_spec_in(:,1,59), pcols, lchnk)
    call outfld('FUL_toa_int60     ',lwup_rad_spec_in(:,1,60), pcols, lchnk)
    call outfld('FDL_toa_int60     ',lwdown_rad_spec_in(:,1,60), pcols, lchnk)
    call outfld('FUL_toa_int61     ',lwup_rad_spec_in(:,1,61), pcols, lchnk)
    call outfld('FDL_toa_int61     ',lwdown_rad_spec_in(:,1,61), pcols, lchnk)
    call outfld('FUL_toa_int62     ',lwup_rad_spec_in(:,1,62), pcols, lchnk)
    call outfld('FDL_toa_int62     ',lwdown_rad_spec_in(:,1,62), pcols, lchnk)
    call outfld('FUL_toa_int63     ',lwup_rad_spec_in(:,1,63), pcols, lchnk)
    call outfld('FDL_toa_int63     ',lwdown_rad_spec_in(:,1,63), pcols, lchnk)
    call outfld('FUL_toa_int64     ',lwup_rad_spec_in(:,1,64), pcols, lchnk)
    call outfld('FDL_toa_int64     ',lwdown_rad_spec_in(:,1,64), pcols, lchnk)
    call outfld('FUL_toa_int65     ',lwup_rad_spec_in(:,1,65), pcols, lchnk)
    call outfld('FDL_toa_int65     ',lwdown_rad_spec_in(:,1,65), pcols, lchnk)
    call outfld('FUL_toa_int66     ',lwup_rad_spec_in(:,1,66), pcols, lchnk)
    call outfld('FDL_toa_int66     ',lwdown_rad_spec_in(:,1,66), pcols, lchnk)
    call outfld('FUL_toa_int67     ',lwup_rad_spec_in(:,1,67), pcols, lchnk)
    call outfld('FDL_toa_int67     ',lwdown_rad_spec_in(:,1,67), pcols, lchnk)
    call outfld('FUL_toa_int68     ',lwup_rad_spec_in(:,1,68), pcols, lchnk)
    call outfld('FDL_toa_int68     ',lwdown_rad_spec_in(:,1,68), pcols, lchnk)

end subroutine outfld_spectral_flux_fullsky



!============================================================================

  subroutine outfld_spectral_flux_clearsky (lchnk, lwdown_rad_spec_in, lwup_rad_spec_in, swup_rad_spec_in, swdown_rad_spec_in )

!------------------------------------------------------------------------
!
! Purpose:  Add spectral output to master field list
!
!------------------------------------------------------------------------
   use cam_history,       only: outfld
!------------------------------------------------------------------------
!
! Arguments
!
    integer,  intent(in) :: lchnk
    real(r8), intent(in), dimension(pcols,pverp,ntot_wavlnrng) ::  lwdown_rad_spec_in
    real(r8), intent(in), dimension(pcols,pverp,ntot_wavlnrng) ::  lwup_rad_spec_in
    real(r8), intent(in), dimension(pcols,pverp,ntot_wavlnrng) ::  swup_rad_spec_in
    real(r8), intent(in), dimension(pcols,pverp,ntot_wavlnrng) ::  swdown_rad_spec_in

!------------------------------------------------------------------------
!
! Start Code
!

    ! shortwave fields
    call outfld('FUSC_int01     ',swup_rad_spec_in(:,:,1), pcols, lchnk)
    call outfld('FDSC_int01     ',swdown_rad_spec_in(:,:,1), pcols, lchnk)
    call outfld('FUSC_int02     ',swup_rad_spec_in(:,:,2), pcols, lchnk)
    call outfld('FDSC_int02     ',swdown_rad_spec_in(:,:,2), pcols, lchnk)
    call outfld('FUSC_int03     ',swup_rad_spec_in(:,:,3), pcols, lchnk)
    call outfld('FDSC_int03     ',swdown_rad_spec_in(:,:,3), pcols, lchnk)
    call outfld('FUSC_int04     ',swup_rad_spec_in(:,:,4), pcols, lchnk)
    call outfld('FDSC_int04     ',swdown_rad_spec_in(:,:,4), pcols, lchnk)
    call outfld('FUSC_int05     ',swup_rad_spec_in(:,:,5), pcols, lchnk)
    call outfld('FDSC_int05     ',swdown_rad_spec_in(:,:,5), pcols, lchnk)
    call outfld('FUSC_int06     ',swup_rad_spec_in(:,:,6), pcols, lchnk)
    call outfld('FDSC_int06     ',swdown_rad_spec_in(:,:,6), pcols, lchnk)
    call outfld('FUSC_int07     ',swup_rad_spec_in(:,:,7), pcols, lchnk)
    call outfld('FDSC_int07     ',swdown_rad_spec_in(:,:,7), pcols, lchnk)
    call outfld('FUSC_int08     ',swup_rad_spec_in(:,:,8), pcols, lchnk)
    call outfld('FDSC_int08     ',swdown_rad_spec_in(:,:,8), pcols, lchnk)
    call outfld('FUSC_int09     ',swup_rad_spec_in(:,:,9), pcols, lchnk)
    call outfld('FDSC_int09     ',swdown_rad_spec_in(:,:,9), pcols, lchnk)
    call outfld('FUSC_int10     ',swup_rad_spec_in(:,:,10), pcols, lchnk)
    call outfld('FDSC_int10     ',swdown_rad_spec_in(:,:,10), pcols, lchnk)
    call outfld('FUSC_int11     ',swup_rad_spec_in(:,:,11), pcols, lchnk)
    call outfld('FDSC_int11     ',swdown_rad_spec_in(:,:,11), pcols, lchnk)
    call outfld('FUSC_int12     ',swup_rad_spec_in(:,:,12), pcols, lchnk)
    call outfld('FDSC_int12     ',swdown_rad_spec_in(:,:,12), pcols, lchnk)
    call outfld('FUSC_int13     ',swup_rad_spec_in(:,:,13), pcols, lchnk)
    call outfld('FDSC_int13     ',swdown_rad_spec_in(:,:,13), pcols, lchnk)
    call outfld('FUSC_int14     ',swup_rad_spec_in(:,:,14), pcols, lchnk)
    call outfld('FDSC_int14     ',swdown_rad_spec_in(:,:,14), pcols, lchnk)
    call outfld('FUSC_int15     ',swup_rad_spec_in(:,:,15), pcols, lchnk)
    call outfld('FDSC_int15     ',swdown_rad_spec_in(:,:,15), pcols, lchnk)
    call outfld('FUSC_int16     ',swup_rad_spec_in(:,:,16), pcols, lchnk)
    call outfld('FDSC_int16     ',swdown_rad_spec_in(:,:,16), pcols, lchnk)
    call outfld('FUSC_int17     ',swup_rad_spec_in(:,:,17), pcols, lchnk)
    call outfld('FDSC_int17    ',swdown_rad_spec_in(:,:,17), pcols, lchnk)
    call outfld('FUSC_int18     ',swup_rad_spec_in(:,:,18), pcols, lchnk)
    call outfld('FDSC_int18     ',swdown_rad_spec_in(:,:,18), pcols, lchnk)
    call outfld('FUSC_int19     ',swup_rad_spec_in(:,:,19), pcols, lchnk)
    call outfld('FDSC_int19     ',swdown_rad_spec_in(:,:,19), pcols, lchnk)
    call outfld('FUSC_int20     ',swup_rad_spec_in(:,:,20), pcols, lchnk)
    call outfld('FDSC_int20     ',swdown_rad_spec_in(:,:,20), pcols, lchnk)
    call outfld('FUSC_int21     ',swup_rad_spec_in(:,:,21), pcols, lchnk)
    call outfld('FDSC_int21     ',swdown_rad_spec_in(:,:,21), pcols, lchnk)
    call outfld('FUSC_int22     ',swup_rad_spec_in(:,:,22), pcols, lchnk)
    call outfld('FDSC_int22     ',swdown_rad_spec_in(:,:,22), pcols, lchnk)
    call outfld('FUSC_int23     ',swup_rad_spec_in(:,:,23), pcols, lchnk)
    call outfld('FDSC_int23     ',swdown_rad_spec_in(:,:,23), pcols, lchnk)
    call outfld('FUSC_int24     ',swup_rad_spec_in(:,:,24), pcols, lchnk)
    call outfld('FDSC_int24     ',swdown_rad_spec_in(:,:,24), pcols, lchnk)
    call outfld('FUSC_int25     ',swup_rad_spec_in(:,:,25), pcols, lchnk)
    call outfld('FDSC_int25     ',swdown_rad_spec_in(:,:,25), pcols, lchnk)
    call outfld('FUSC_int26     ',swup_rad_spec_in(:,:,26), pcols, lchnk)
    call outfld('FDSC_int26     ',swdown_rad_spec_in(:,:,26), pcols, lchnk)
    call outfld('FUSC_int27     ',swup_rad_spec_in(:,:,27), pcols, lchnk)
    call outfld('FDSC_int27     ',swdown_rad_spec_in(:,:,27), pcols, lchnk)
    call outfld('FUSC_int28     ',swup_rad_spec_in(:,:,28), pcols, lchnk)
    call outfld('FDSC_int28     ',swdown_rad_spec_in(:,:,28), pcols, lchnk)
    call outfld('FUSC_int29     ',swup_rad_spec_in(:,:,29), pcols, lchnk)
    call outfld('FDSC_int29     ',swdown_rad_spec_in(:,:,29), pcols, lchnk)
    call outfld('FUSC_int30     ',swup_rad_spec_in(:,:,30), pcols, lchnk)
    call outfld('FDSC_int30     ',swdown_rad_spec_in(:,:,30), pcols, lchnk)
    call outfld('FUSC_int31     ',swup_rad_spec_in(:,:,31), pcols, lchnk)
    call outfld('FDSC_int31     ',swdown_rad_spec_in(:,:,31), pcols, lchnk)
    call outfld('FUSC_int32     ',swup_rad_spec_in(:,:,32), pcols, lchnk)
    call outfld('FDSC_int32     ',swdown_rad_spec_in(:,:,32), pcols, lchnk)
    call outfld('FUSC_int33     ',swup_rad_spec_in(:,:,33), pcols, lchnk)
    call outfld('FDSC_int33     ',swdown_rad_spec_in(:,:,33), pcols, lchnk)
    call outfld('FUSC_int34     ',swup_rad_spec_in(:,:,34), pcols, lchnk)
    call outfld('FDSC_int34     ',swdown_rad_spec_in(:,:,34), pcols, lchnk)
    call outfld('FUSC_int35     ',swup_rad_spec_in(:,:,35), pcols, lchnk)
    call outfld('FDSC_int35     ',swdown_rad_spec_in(:,:,35), pcols, lchnk)
    call outfld('FUSC_int36     ',swup_rad_spec_in(:,:,36), pcols, lchnk)
    call outfld('FDSC_int36     ',swdown_rad_spec_in(:,:,36), pcols, lchnk)
    call outfld('FUSC_int37     ',swup_rad_spec_in(:,:,37), pcols, lchnk)
    call outfld('FDSC_int37     ',swdown_rad_spec_in(:,:,37), pcols, lchnk)
    call outfld('FUSC_int38     ',swup_rad_spec_in(:,:,38), pcols, lchnk)
    call outfld('FDSC_int38     ',swdown_rad_spec_in(:,:,38), pcols, lchnk)
    call outfld('FUSC_int39     ',swup_rad_spec_in(:,:,39), pcols, lchnk)
    call outfld('FDSC_int39     ',swdown_rad_spec_in(:,:,39), pcols, lchnk)
    call outfld('FUSC_int40     ',swup_rad_spec_in(:,:,40), pcols, lchnk)
    call outfld('FDSC_int40     ',swdown_rad_spec_in(:,:,40), pcols, lchnk)
    call outfld('FUSC_int41     ',swup_rad_spec_in(:,:,41), pcols, lchnk)
    call outfld('FDSC_int41     ',swdown_rad_spec_in(:,:,41), pcols, lchnk)
    call outfld('FUSC_int42     ',swup_rad_spec_in(:,:,42), pcols, lchnk)
    call outfld('FDSC_int42     ',swdown_rad_spec_in(:,:,42), pcols, lchnk)
    call outfld('FUSC_int43     ',swup_rad_spec_in(:,:,43), pcols, lchnk)
    call outfld('FDSC_int43     ',swdown_rad_spec_in(:,:,43), pcols, lchnk)
    call outfld('FUSC_int44     ',swup_rad_spec_in(:,:,44), pcols, lchnk)
    call outfld('FDSC_int44     ',swdown_rad_spec_in(:,:,44), pcols, lchnk)
    call outfld('FUSC_int45     ',swup_rad_spec_in(:,:,45), pcols, lchnk)
    call outfld('FDSC_int45     ',swdown_rad_spec_in(:,:,45), pcols, lchnk)
    call outfld('FUSC_int46     ',swup_rad_spec_in(:,:,46), pcols, lchnk)
    call outfld('FDSC_int46     ',swdown_rad_spec_in(:,:,46), pcols, lchnk)
    call outfld('FUSC_int47     ',swup_rad_spec_in(:,:,47), pcols, lchnk)
    call outfld('FDSC_int47     ',swdown_rad_spec_in(:,:,47), pcols, lchnk)
    call outfld('FUSC_int48     ',swup_rad_spec_in(:,:,48), pcols, lchnk)
    call outfld('FDSC_int48     ',swdown_rad_spec_in(:,:,48), pcols, lchnk)
    call outfld('FUSC_int49     ',swup_rad_spec_in(:,:,49), pcols, lchnk)
    call outfld('FDSC_int49     ',swdown_rad_spec_in(:,:,49), pcols, lchnk)
    call outfld('FUSC_int50     ',swup_rad_spec_in(:,:,50), pcols, lchnk)
    call outfld('FDSC_int50     ',swdown_rad_spec_in(:,:,50), pcols, lchnk)
    call outfld('FUSC_int51     ',swup_rad_spec_in(:,:,51), pcols, lchnk)
    call outfld('FDSC_int51     ',swdown_rad_spec_in(:,:,51), pcols, lchnk)
    call outfld('FUSC_int52     ',swup_rad_spec_in(:,:,52), pcols, lchnk)
    call outfld('FDSC_int52     ',swdown_rad_spec_in(:,:,52), pcols, lchnk)
    call outfld('FUSC_int53     ',swup_rad_spec_in(:,:,53), pcols, lchnk)
    call outfld('FDSC_int53     ',swdown_rad_spec_in(:,:,53), pcols, lchnk)
    call outfld('FUSC_int54     ',swup_rad_spec_in(:,:,54), pcols, lchnk)
    call outfld('FDSC_int54     ',swdown_rad_spec_in(:,:,54), pcols, lchnk)
    call outfld('FUSC_int55     ',swup_rad_spec_in(:,:,55), pcols, lchnk)
    call outfld('FDSC_int55     ',swdown_rad_spec_in(:,:,55), pcols, lchnk)
    call outfld('FUSC_int56     ',swup_rad_spec_in(:,:,56), pcols, lchnk)
    call outfld('FDSC_int56     ',swdown_rad_spec_in(:,:,56), pcols, lchnk)
    call outfld('FUSC_int57     ',swup_rad_spec_in(:,:,57), pcols, lchnk)
    call outfld('FDSC_int57     ',swdown_rad_spec_in(:,:,57), pcols, lchnk)
    call outfld('FUSC_int58     ',swup_rad_spec_in(:,:,58), pcols, lchnk)
    call outfld('FDSC_int58     ',swdown_rad_spec_in(:,:,58), pcols, lchnk)
    call outfld('FUSC_int59     ',swup_rad_spec_in(:,:,59), pcols, lchnk)
    call outfld('FDSC_int59     ',swdown_rad_spec_in(:,:,59), pcols, lchnk)
    call outfld('FUSC_int60     ',swup_rad_spec_in(:,:,60), pcols, lchnk)
    call outfld('FDSC_int60     ',swdown_rad_spec_in(:,:,60), pcols, lchnk)
    call outfld('FUSC_int61     ',swup_rad_spec_in(:,:,61), pcols, lchnk)
    call outfld('FDSC_int61     ',swdown_rad_spec_in(:,:,61), pcols, lchnk)
    call outfld('FUSC_int62     ',swup_rad_spec_in(:,:,62), pcols, lchnk)
    call outfld('FDSC_int62     ',swdown_rad_spec_in(:,:,62), pcols, lchnk)
    call outfld('FUSC_int63     ',swup_rad_spec_in(:,:,63), pcols, lchnk)
    call outfld('FDSC_int63     ',swdown_rad_spec_in(:,:,63), pcols, lchnk)
    call outfld('FUSC_int64     ',swup_rad_spec_in(:,:,64), pcols, lchnk)
    call outfld('FDSC_int64     ',swdown_rad_spec_in(:,:,64), pcols, lchnk)
    call outfld('FUSC_int65     ',swup_rad_spec_in(:,:,65), pcols, lchnk)
    call outfld('FDSC_int65     ',swdown_rad_spec_in(:,:,65), pcols, lchnk)
    call outfld('FUSC_int66     ',swup_rad_spec_in(:,:,66), pcols, lchnk)
    call outfld('FDSC_int66     ',swdown_rad_spec_in(:,:,66), pcols, lchnk)
    call outfld('FUSC_int67     ',swup_rad_spec_in(:,:,67), pcols, lchnk)
    call outfld('FDSC_int67     ',swdown_rad_spec_in(:,:,67), pcols, lchnk)
    call outfld('FUSC_int68     ',swup_rad_spec_in(:,:,68), pcols, lchnk)
    call outfld('FDSC_int68     ',swdown_rad_spec_in(:,:,68), pcols, lchnk)

    call outfld('FUSC_toa_int01     ',swup_rad_spec_in(:,1,1), pcols, lchnk)
    call outfld('FDSC_toa_int01     ',swdown_rad_spec_in(:,1,1), pcols, lchnk)
    call outfld('FUSC_toa_int02     ',swup_rad_spec_in(:,1,2), pcols, lchnk)
    call outfld('FDSC_toa_int02     ',swdown_rad_spec_in(:,1,2), pcols, lchnk)
    call outfld('FUSC_toa_int03     ',swup_rad_spec_in(:,1,3), pcols, lchnk)
    call outfld('FDSC_toa_int03     ',swdown_rad_spec_in(:,1,3), pcols, lchnk)
    call outfld('FUSC_toa_int04     ',swup_rad_spec_in(:,1,4), pcols, lchnk)
    call outfld('FDSC_toa_int04     ',swdown_rad_spec_in(:,1,4), pcols, lchnk)
    call outfld('FUSC_toa_int05     ',swup_rad_spec_in(:,1,5), pcols, lchnk)
    call outfld('FDSC_toa_int05     ',swdown_rad_spec_in(:,1,5), pcols, lchnk)
    call outfld('FUSC_toa_int06     ',swup_rad_spec_in(:,1,6), pcols, lchnk)
    call outfld('FDSC_toa_int06     ',swdown_rad_spec_in(:,1,6), pcols, lchnk)
    call outfld('FUSC_toa_int07     ',swup_rad_spec_in(:,1,7), pcols, lchnk)
    call outfld('FDSC_toa_int07     ',swdown_rad_spec_in(:,1,7), pcols, lchnk)
    call outfld('FUSC_toa_int08     ',swup_rad_spec_in(:,1,8), pcols, lchnk)
    call outfld('FDSC_toa_int08     ',swdown_rad_spec_in(:,1,8), pcols, lchnk)
    call outfld('FUSC_toa_int09     ',swup_rad_spec_in(:,1,9), pcols, lchnk)
    call outfld('FDSC_toa_int09     ',swdown_rad_spec_in(:,1,9), pcols, lchnk)
    call outfld('FUSC_toa_int10     ',swup_rad_spec_in(:,1,10), pcols, lchnk)
    call outfld('FDSC_toa_int10     ',swdown_rad_spec_in(:,1,10), pcols, lchnk)
    call outfld('FUSC_toa_int11     ',swup_rad_spec_in(:,1,11), pcols, lchnk)
    call outfld('FDSC_toa_int11     ',swdown_rad_spec_in(:,1,11), pcols, lchnk)
    call outfld('FUSC_toa_int12     ',swup_rad_spec_in(:,1,12), pcols, lchnk)
    call outfld('FDSC_toa_int12     ',swdown_rad_spec_in(:,1,12), pcols, lchnk)
    call outfld('FUSC_toa_int13     ',swup_rad_spec_in(:,1,13), pcols, lchnk)
    call outfld('FDSC_toa_int13     ',swdown_rad_spec_in(:,1,13), pcols, lchnk)
    call outfld('FUSC_toa_int14     ',swup_rad_spec_in(:,1,14), pcols, lchnk)
    call outfld('FDSC_toa_int14     ',swdown_rad_spec_in(:,1,14), pcols, lchnk)
    call outfld('FUSC_toa_int15     ',swup_rad_spec_in(:,1,15), pcols, lchnk)
    call outfld('FDSC_toa_int15     ',swdown_rad_spec_in(:,1,15), pcols, lchnk)
    call outfld('FUSC_toa_int16     ',swup_rad_spec_in(:,1,16), pcols, lchnk)
    call outfld('FDSC_toa_int16     ',swdown_rad_spec_in(:,1,16), pcols, lchnk)
    call outfld('FUSC_toa_int17     ',swup_rad_spec_in(:,1,17), pcols, lchnk)
    call outfld('FDSC_toa_int17     ',swdown_rad_spec_in(:,1,17), pcols, lchnk)
    call outfld('FUSC_toa_int18     ',swup_rad_spec_in(:,1,18), pcols, lchnk)
    call outfld('FDSC_toa_int18     ',swdown_rad_spec_in(:,1,18), pcols, lchnk)
    call outfld('FUSC_toa_int19     ',swup_rad_spec_in(:,1,19), pcols, lchnk)
    call outfld('FDSC_toa_int19     ',swdown_rad_spec_in(:,1,19), pcols, lchnk)
    call outfld('FUSC_toa_int20     ',swup_rad_spec_in(:,1,20), pcols, lchnk)
    call outfld('FDSC_toa_int20     ',swdown_rad_spec_in(:,1,20), pcols, lchnk)
    call outfld('FUSC_toa_int21     ',swup_rad_spec_in(:,1,21), pcols, lchnk)
    call outfld('FDSC_toa_int21     ',swdown_rad_spec_in(:,1,21), pcols, lchnk)
    call outfld('FUSC_toa_int22     ',swup_rad_spec_in(:,1,22), pcols, lchnk)
    call outfld('FDSC_toa_int22     ',swdown_rad_spec_in(:,1,22), pcols, lchnk)
    call outfld('FUSC_toa_int23     ',swup_rad_spec_in(:,1,23), pcols, lchnk)
    call outfld('FDSC_toa_int23     ',swdown_rad_spec_in(:,1,23), pcols, lchnk)
    call outfld('FUSC_toa_int24     ',swup_rad_spec_in(:,1,24), pcols, lchnk)
    call outfld('FDSC_toa_int24     ',swdown_rad_spec_in(:,1,24), pcols, lchnk)
    call outfld('FUSC_toa_int25     ',swup_rad_spec_in(:,1,25), pcols, lchnk)
    call outfld('FDSC_toa_int25     ',swdown_rad_spec_in(:,1,25), pcols, lchnk)
    call outfld('FUSC_toa_int26     ',swup_rad_spec_in(:,1,26), pcols, lchnk)
    call outfld('FDSC_toa_int26     ',swdown_rad_spec_in(:,1,26), pcols, lchnk)
    call outfld('FUSC_toa_int27     ',swup_rad_spec_in(:,1,27), pcols, lchnk)
    call outfld('FDSC_toa_int27     ',swdown_rad_spec_in(:,1,27), pcols, lchnk)
    call outfld('FUSC_toa_int28     ',swup_rad_spec_in(:,1,28), pcols, lchnk)
    call outfld('FDSC_toa_int28     ',swdown_rad_spec_in(:,1,28), pcols, lchnk)
    call outfld('FUSC_toa_int29     ',swup_rad_spec_in(:,1,29), pcols, lchnk)
    call outfld('FDSC_toa_int29     ',swdown_rad_spec_in(:,1,29), pcols, lchnk)
    call outfld('FUSC_toa_int30     ',swup_rad_spec_in(:,1,30), pcols, lchnk)
    call outfld('FDSC_toa_int30     ',swdown_rad_spec_in(:,1,30), pcols, lchnk)
    call outfld('FUSC_toa_int31     ',swup_rad_spec_in(:,1,31), pcols, lchnk)
    call outfld('FDSC_toa_int31     ',swdown_rad_spec_in(:,1,31), pcols, lchnk)
    call outfld('FUSC_toa_int32     ',swup_rad_spec_in(:,1,32), pcols, lchnk)
    call outfld('FDSC_toa_int32     ',swdown_rad_spec_in(:,1,32), pcols, lchnk)
    call outfld('FUSC_toa_int33     ',swup_rad_spec_in(:,1,33), pcols, lchnk)
    call outfld('FDSC_toa_int33     ',swdown_rad_spec_in(:,1,33), pcols, lchnk)
    call outfld('FUSC_toa_int34     ',swup_rad_spec_in(:,1,34), pcols, lchnk)
    call outfld('FDSC_toa_int34     ',swdown_rad_spec_in(:,1,34), pcols, lchnk)
    call outfld('FUSC_toa_int35     ',swup_rad_spec_in(:,1,35), pcols, lchnk)
    call outfld('FDSC_toa_int35     ',swdown_rad_spec_in(:,1,35), pcols, lchnk)
    call outfld('FUSC_toa_int36     ',swup_rad_spec_in(:,1,36), pcols, lchnk)
    call outfld('FDSC_toa_int36     ',swdown_rad_spec_in(:,1,36), pcols, lchnk)
    call outfld('FUSC_toa_int37     ',swup_rad_spec_in(:,1,37), pcols, lchnk)
    call outfld('FDSC_toa_int37     ',swdown_rad_spec_in(:,1,37), pcols, lchnk)
    call outfld('FUSC_toa_int38     ',swup_rad_spec_in(:,1,38), pcols, lchnk)
    call outfld('FDSC_toa_int38     ',swdown_rad_spec_in(:,1,38), pcols, lchnk)
    call outfld('FUSC_toa_int39     ',swup_rad_spec_in(:,1,39), pcols, lchnk)
    call outfld('FDSC_toa_int39     ',swdown_rad_spec_in(:,1,39), pcols, lchnk)
    call outfld('FUSC_toa_int40     ',swup_rad_spec_in(:,1,40), pcols, lchnk)
    call outfld('FDSC_toa_int40     ',swdown_rad_spec_in(:,1,40), pcols, lchnk)
    call outfld('FUSC_toa_int41     ',swup_rad_spec_in(:,1,41), pcols, lchnk)
    call outfld('FDSC_toa_int41     ',swdown_rad_spec_in(:,1,41), pcols, lchnk)
    call outfld('FUSC_toa_int42     ',swup_rad_spec_in(:,1,42), pcols, lchnk)
    call outfld('FDSC_toa_int42     ',swdown_rad_spec_in(:,1,42), pcols, lchnk)
    call outfld('FUSC_toa_int43     ',swup_rad_spec_in(:,1,43), pcols, lchnk)
    call outfld('FDSC_toa_int43     ',swdown_rad_spec_in(:,1,43), pcols, lchnk)
    call outfld('FUSC_toa_int44     ',swup_rad_spec_in(:,1,44), pcols, lchnk)
    call outfld('FDSC_toa_int44     ',swdown_rad_spec_in(:,1,44), pcols, lchnk)
    call outfld('FUSC_toa_int45     ',swup_rad_spec_in(:,1,45), pcols, lchnk)
    call outfld('FDSC_toa_int45     ',swdown_rad_spec_in(:,1,45), pcols, lchnk)
    call outfld('FUSC_toa_int46     ',swup_rad_spec_in(:,1,46), pcols, lchnk)
    call outfld('FDSC_toa_int46     ',swdown_rad_spec_in(:,1,46), pcols, lchnk)
    call outfld('FUSC_toa_int47     ',swup_rad_spec_in(:,1,47), pcols, lchnk)
    call outfld('FDSC_toa_int47     ',swdown_rad_spec_in(:,1,47), pcols, lchnk)
    call outfld('FUSC_toa_int48     ',swup_rad_spec_in(:,1,48), pcols, lchnk)
    call outfld('FDSC_toa_int48     ',swdown_rad_spec_in(:,1,48), pcols, lchnk)
    call outfld('FUSC_toa_int49     ',swup_rad_spec_in(:,1,49), pcols, lchnk)
    call outfld('FDSC_toa_int49     ',swdown_rad_spec_in(:,1,49), pcols, lchnk)
    call outfld('FUSC_toa_int50     ',swup_rad_spec_in(:,1,50), pcols, lchnk)
    call outfld('FDSC_toa_int50     ',swdown_rad_spec_in(:,1,50), pcols, lchnk)
    call outfld('FUSC_toa_int51     ',swup_rad_spec_in(:,1,51), pcols, lchnk)
    call outfld('FDSC_toa_int51     ',swdown_rad_spec_in(:,1,51), pcols, lchnk)
    call outfld('FUSC_toa_int52     ',swup_rad_spec_in(:,1,52), pcols, lchnk)
    call outfld('FDSC_toa_int52     ',swdown_rad_spec_in(:,1,52), pcols, lchnk)
    call outfld('FUSC_toa_int53     ',swup_rad_spec_in(:,1,53), pcols, lchnk)
    call outfld('FDSC_toa_int53     ',swdown_rad_spec_in(:,1,53), pcols, lchnk)
    call outfld('FUSC_toa_int54     ',swup_rad_spec_in(:,1,54), pcols, lchnk)
    call outfld('FDSC_toa_int54     ',swdown_rad_spec_in(:,1,54), pcols, lchnk)
    call outfld('FUSC_toa_int55     ',swup_rad_spec_in(:,1,55), pcols, lchnk)
    call outfld('FDSC_toa_int55     ',swdown_rad_spec_in(:,1,55), pcols, lchnk)
    call outfld('FUSC_toa_int56     ',swup_rad_spec_in(:,1,56), pcols, lchnk)
    call outfld('FDSC_toa_int56     ',swdown_rad_spec_in(:,1,56), pcols, lchnk)
    call outfld('FUSC_toa_int57     ',swup_rad_spec_in(:,1,57), pcols, lchnk)
    call outfld('FDSC_toa_int57     ',swdown_rad_spec_in(:,1,57), pcols, lchnk)
    call outfld('FUSC_toa_int58     ',swup_rad_spec_in(:,1,58), pcols, lchnk)
    call outfld('FDSC_toa_int58     ',swdown_rad_spec_in(:,1,58), pcols, lchnk)
    call outfld('FUSC_toa_int59     ',swup_rad_spec_in(:,1,59), pcols, lchnk)
    call outfld('FDSC_toa_int59     ',swdown_rad_spec_in(:,1,59), pcols, lchnk)
    call outfld('FUSC_toa_int60     ',swup_rad_spec_in(:,1,60), pcols, lchnk)
    call outfld('FDSC_toa_int60     ',swdown_rad_spec_in(:,1,60), pcols, lchnk)
    call outfld('FUSC_toa_int61     ',swup_rad_spec_in(:,1,61), pcols, lchnk)
    call outfld('FDSC_toa_int61     ',swdown_rad_spec_in(:,1,61), pcols, lchnk)
    call outfld('FUSC_toa_int62     ',swup_rad_spec_in(:,1,62), pcols, lchnk)
    call outfld('FDSC_toa_int62     ',swdown_rad_spec_in(:,1,62), pcols, lchnk)
    call outfld('FUSC_toa_int63     ',swup_rad_spec_in(:,1,63), pcols, lchnk)
    call outfld('FDSC_toa_int63     ',swdown_rad_spec_in(:,1,63), pcols, lchnk)
    call outfld('FUSC_toa_int64     ',swup_rad_spec_in(:,1,64), pcols, lchnk)
    call outfld('FDSC_toa_int64     ',swdown_rad_spec_in(:,1,64), pcols, lchnk)
    call outfld('FUSC_toa_int65     ',swup_rad_spec_in(:,1,65), pcols, lchnk)
    call outfld('FDSC_toa_int65     ',swdown_rad_spec_in(:,1,65), pcols, lchnk)
    call outfld('FUSC_toa_int66     ',swup_rad_spec_in(:,1,66), pcols, lchnk)
    call outfld('FDSC_toa_int66     ',swdown_rad_spec_in(:,1,66), pcols, lchnk)
    call outfld('FUSC_toa_int67     ',swup_rad_spec_in(:,1,67), pcols, lchnk)
    call outfld('FDSC_toa_int67     ',swdown_rad_spec_in(:,1,67), pcols, lchnk)
    call outfld('FUSC_toa_int68     ',swup_rad_spec_in(:,1,68), pcols, lchnk)
    call outfld('FDSC_toa_int68     ',swdown_rad_spec_in(:,1,68), pcols, lchnk)


    ! longwave
    call outfld('FULC_int01     ',lwup_rad_spec_in(:,:,1), pcols, lchnk)
    call outfld('FDLC_int01     ',lwdown_rad_spec_in(:,:,1), pcols, lchnk)
    call outfld('FULC_int02     ',lwup_rad_spec_in(:,:,2), pcols, lchnk)
    call outfld('FDLC_int02     ',lwdown_rad_spec_in(:,:,2), pcols, lchnk)
    call outfld('FULC_int03     ',lwup_rad_spec_in(:,:,3), pcols, lchnk)
    call outfld('FDLC_int03     ',lwdown_rad_spec_in(:,:,3), pcols, lchnk)
    call outfld('FULC_int04     ',lwup_rad_spec_in(:,:,4), pcols, lchnk)
    call outfld('FDLC_int04     ',lwdown_rad_spec_in(:,:,4), pcols, lchnk)
    call outfld('FULC_int05     ',lwup_rad_spec_in(:,:,5), pcols, lchnk)
    call outfld('FDLC_int05     ',lwdown_rad_spec_in(:,:,5), pcols, lchnk)
    call outfld('FULC_int06     ',lwup_rad_spec_in(:,:,6), pcols, lchnk)
    call outfld('FDLC_int06     ',lwdown_rad_spec_in(:,:,6), pcols, lchnk)
    call outfld('FULC_int07     ',lwup_rad_spec_in(:,:,7), pcols, lchnk)
    call outfld('FDLC_int07     ',lwdown_rad_spec_in(:,:,7), pcols, lchnk)
    call outfld('FULC_int08     ',lwup_rad_spec_in(:,:,8), pcols, lchnk)
    call outfld('FDLC_int08     ',lwdown_rad_spec_in(:,:,8), pcols, lchnk)
    call outfld('FULC_int09     ',lwup_rad_spec_in(:,:,9), pcols, lchnk)
    call outfld('FDLC_int09     ',lwdown_rad_spec_in(:,:,9), pcols, lchnk)
    call outfld('FULC_int10     ',lwup_rad_spec_in(:,:,10), pcols, lchnk)
    call outfld('FDLC_int10     ',lwdown_rad_spec_in(:,:,10), pcols, lchnk)
    call outfld('FULC_int11     ',lwup_rad_spec_in(:,:,11), pcols, lchnk)
    call outfld('FDLC_int11     ',lwdown_rad_spec_in(:,:,11), pcols, lchnk)
    call outfld('FULC_int12     ',lwup_rad_spec_in(:,:,12), pcols, lchnk)
    call outfld('FDLC_int12     ',lwdown_rad_spec_in(:,:,12), pcols, lchnk)
    call outfld('FULC_int13     ',lwup_rad_spec_in(:,:,13), pcols, lchnk)
    call outfld('FDLC_int13     ',lwdown_rad_spec_in(:,:,13), pcols, lchnk)
    call outfld('FULC_int14     ',lwup_rad_spec_in(:,:,14), pcols, lchnk)
    call outfld('FDLC_int14     ',lwdown_rad_spec_in(:,:,14), pcols, lchnk)
    call outfld('FULC_int15     ',lwup_rad_spec_in(:,:,15), pcols, lchnk)
    call outfld('FDLC_int15     ',lwdown_rad_spec_in(:,:,15), pcols, lchnk)
    call outfld('FULC_int16     ',lwup_rad_spec_in(:,:,16), pcols, lchnk)
    call outfld('FDLC_int16     ',lwdown_rad_spec_in(:,:,16), pcols, lchnk)
    call outfld('FULC_int17     ',lwup_rad_spec_in(:,:,17), pcols, lchnk)
    call outfld('FDLC_int17     ',lwdown_rad_spec_in(:,:,17), pcols, lchnk)
    call outfld('FULC_int18     ',lwup_rad_spec_in(:,:,18), pcols, lchnk)
    call outfld('FDLC_int18     ',lwdown_rad_spec_in(:,:,18), pcols, lchnk)
    call outfld('FULC_int19     ',lwup_rad_spec_in(:,:,19), pcols, lchnk)
    call outfld('FDLC_int19     ',lwdown_rad_spec_in(:,:,19), pcols, lchnk)
    call outfld('FULC_int20     ',lwup_rad_spec_in(:,:,20), pcols, lchnk)
    call outfld('FDLC_int20     ',lwdown_rad_spec_in(:,:,20), pcols, lchnk)
    call outfld('FULC_int21     ',lwup_rad_spec_in(:,:,21), pcols, lchnk)
    call outfld('FDLC_int21     ',lwdown_rad_spec_in(:,:,21), pcols, lchnk)
    call outfld('FULC_int22     ',lwup_rad_spec_in(:,:,22), pcols, lchnk)
    call outfld('FDLC_int22     ',lwdown_rad_spec_in(:,:,22), pcols, lchnk)
    call outfld('FULC_int23     ',lwup_rad_spec_in(:,:,23), pcols, lchnk)
    call outfld('FDLC_int23     ',lwdown_rad_spec_in(:,:,23), pcols, lchnk)
    call outfld('FULC_int24     ',lwup_rad_spec_in(:,:,24), pcols, lchnk)
    call outfld('FDLC_int24     ',lwdown_rad_spec_in(:,:,24), pcols, lchnk)
    call outfld('FULC_int25     ',lwup_rad_spec_in(:,:,25), pcols, lchnk)
    call outfld('FDLC_int25     ',lwdown_rad_spec_in(:,:,25), pcols, lchnk)
    call outfld('FULC_int26     ',lwup_rad_spec_in(:,:,26), pcols, lchnk)
    call outfld('FDLC_int26     ',lwdown_rad_spec_in(:,:,26), pcols, lchnk)
    call outfld('FULC_int27     ',lwup_rad_spec_in(:,:,27), pcols, lchnk)
    call outfld('FDLC_int27     ',lwdown_rad_spec_in(:,:,27), pcols, lchnk)
    call outfld('FULC_int28     ',lwup_rad_spec_in(:,:,28), pcols, lchnk)
    call outfld('FDLC_int28     ',lwdown_rad_spec_in(:,:,28), pcols, lchnk)
    call outfld('FULC_int29     ',lwup_rad_spec_in(:,:,29), pcols, lchnk)
    call outfld('FDLC_int29     ',lwdown_rad_spec_in(:,:,29), pcols, lchnk)
    call outfld('FULC_int30     ',lwup_rad_spec_in(:,:,30), pcols, lchnk)
    call outfld('FDLC_int30     ',lwdown_rad_spec_in(:,:,30), pcols, lchnk)
    call outfld('FULC_int31     ',lwup_rad_spec_in(:,:,31), pcols, lchnk)
    call outfld('FDLC_int31     ',lwdown_rad_spec_in(:,:,31), pcols, lchnk)
    call outfld('FULC_int32     ',lwup_rad_spec_in(:,:,32), pcols, lchnk)
    call outfld('FDLC_int32     ',lwdown_rad_spec_in(:,:,32), pcols, lchnk)
    call outfld('FULC_int33     ',lwup_rad_spec_in(:,:,33), pcols, lchnk)
    call outfld('FDLC_int33     ',lwdown_rad_spec_in(:,:,33), pcols, lchnk)
    call outfld('FULC_int34     ',lwup_rad_spec_in(:,:,34), pcols, lchnk)
    call outfld('FDLC_int34     ',lwdown_rad_spec_in(:,:,34), pcols, lchnk)
    call outfld('FULC_int35     ',lwup_rad_spec_in(:,:,35), pcols, lchnk)
    call outfld('FDLC_int35     ',lwdown_rad_spec_in(:,:,35), pcols, lchnk)
    call outfld('FULC_int36     ',lwup_rad_spec_in(:,:,36), pcols, lchnk)
    call outfld('FDLC_int36     ',lwdown_rad_spec_in(:,:,36), pcols, lchnk)
    call outfld('FULC_int37     ',lwup_rad_spec_in(:,:,37), pcols, lchnk)
    call outfld('FDLC_int37     ',lwdown_rad_spec_in(:,:,37), pcols, lchnk)
    call outfld('FULC_int38     ',lwup_rad_spec_in(:,:,38), pcols, lchnk)
    call outfld('FDLC_int38     ',lwdown_rad_spec_in(:,:,38), pcols, lchnk)
    call outfld('FULC_int39     ',lwup_rad_spec_in(:,:,39), pcols, lchnk)
    call outfld('FDLC_int39     ',lwdown_rad_spec_in(:,:,39), pcols, lchnk)
    call outfld('FULC_int40     ',lwup_rad_spec_in(:,:,40), pcols, lchnk)
    call outfld('FDLC_int40     ',lwdown_rad_spec_in(:,:,40), pcols, lchnk)
    call outfld('FULC_int41     ',lwup_rad_spec_in(:,:,41), pcols, lchnk)
    call outfld('FDLC_int41     ',lwdown_rad_spec_in(:,:,41), pcols, lchnk)
    call outfld('FULC_int42     ',lwup_rad_spec_in(:,:,42), pcols, lchnk)
    call outfld('FDLC_int42     ',lwdown_rad_spec_in(:,:,42), pcols, lchnk)
    call outfld('FULC_int43     ',lwup_rad_spec_in(:,:,43), pcols, lchnk)
    call outfld('FDLC_int43     ',lwdown_rad_spec_in(:,:,43), pcols, lchnk)
    call outfld('FULC_int44     ',lwup_rad_spec_in(:,:,44), pcols, lchnk)
    call outfld('FDLC_int44     ',lwdown_rad_spec_in(:,:,44), pcols, lchnk)
    call outfld('FULC_int45     ',lwup_rad_spec_in(:,:,45), pcols, lchnk)
    call outfld('FDLC_int45     ',lwdown_rad_spec_in(:,:,45), pcols, lchnk)
    call outfld('FULC_int46     ',lwup_rad_spec_in(:,:,46), pcols, lchnk)
    call outfld('FDLC_int46     ',lwdown_rad_spec_in(:,:,46), pcols, lchnk)
    call outfld('FULC_int47     ',lwup_rad_spec_in(:,:,47), pcols, lchnk)
    call outfld('FDLC_int47     ',lwdown_rad_spec_in(:,:,47), pcols, lchnk)
    call outfld('FULC_int48     ',lwup_rad_spec_in(:,:,48), pcols, lchnk)
    call outfld('FDLC_int48     ',lwdown_rad_spec_in(:,:,48), pcols, lchnk)
    call outfld('FULC_int49     ',lwup_rad_spec_in(:,:,49), pcols, lchnk)
    call outfld('FDLC_int49     ',lwdown_rad_spec_in(:,:,49), pcols, lchnk)
    call outfld('FULC_int50     ',lwup_rad_spec_in(:,:,50), pcols, lchnk)
    call outfld('FDLC_int50     ',lwdown_rad_spec_in(:,:,50), pcols, lchnk)
    call outfld('FULC_int51     ',lwup_rad_spec_in(:,:,51), pcols, lchnk)
    call outfld('FDLC_int51     ',lwdown_rad_spec_in(:,:,51), pcols, lchnk)
    call outfld('FULC_int52     ',lwup_rad_spec_in(:,:,52), pcols, lchnk)
    call outfld('FDLC_int52     ',lwdown_rad_spec_in(:,:,52), pcols, lchnk)
    call outfld('FULC_int53     ',lwup_rad_spec_in(:,:,53), pcols, lchnk)
    call outfld('FDLC_int53     ',lwdown_rad_spec_in(:,:,53), pcols, lchnk)
    call outfld('FULC_int54     ',lwup_rad_spec_in(:,:,54), pcols, lchnk)
    call outfld('FDLC_int54     ',lwdown_rad_spec_in(:,:,54), pcols, lchnk)
    call outfld('FULC_int55     ',lwup_rad_spec_in(:,:,55), pcols, lchnk)
    call outfld('FDLC_int55     ',lwdown_rad_spec_in(:,:,55), pcols, lchnk)
    call outfld('FULC_int56     ',lwup_rad_spec_in(:,:,56), pcols, lchnk)
    call outfld('FDLC_int56     ',lwdown_rad_spec_in(:,:,56), pcols, lchnk)
    call outfld('FULC_int57     ',lwup_rad_spec_in(:,:,57), pcols, lchnk)
    call outfld('FDLC_int57     ',lwdown_rad_spec_in(:,:,57), pcols, lchnk)
    call outfld('FULC_int58     ',lwup_rad_spec_in(:,:,58), pcols, lchnk)
    call outfld('FDLC_int58     ',lwdown_rad_spec_in(:,:,58), pcols, lchnk)
    call outfld('FULC_int59     ',lwup_rad_spec_in(:,:,59), pcols, lchnk)
    call outfld('FDLC_int59     ',lwdown_rad_spec_in(:,:,59), pcols, lchnk)
    call outfld('FULC_int60     ',lwup_rad_spec_in(:,:,60), pcols, lchnk)
    call outfld('FDLC_int60     ',lwdown_rad_spec_in(:,:,60), pcols, lchnk)
    call outfld('FULC_int61     ',lwup_rad_spec_in(:,:,61), pcols, lchnk)
    call outfld('FDLC_int61     ',lwdown_rad_spec_in(:,:,61), pcols, lchnk)
    call outfld('FULC_int62     ',lwup_rad_spec_in(:,:,62), pcols, lchnk)
    call outfld('FDLC_int62     ',lwdown_rad_spec_in(:,:,62), pcols, lchnk)
    call outfld('FULC_int63     ',lwup_rad_spec_in(:,:,63), pcols, lchnk)
    call outfld('FDLC_int63     ',lwdown_rad_spec_in(:,:,63), pcols, lchnk)
    call outfld('FULC_int64     ',lwup_rad_spec_in(:,:,64), pcols, lchnk)
    call outfld('FDLC_int64     ',lwdown_rad_spec_in(:,:,64), pcols, lchnk)
    call outfld('FULC_int65     ',lwup_rad_spec_in(:,:,65), pcols, lchnk)
    call outfld('FDLC_int65     ',lwdown_rad_spec_in(:,:,65), pcols, lchnk)
    call outfld('FULC_int66     ',lwup_rad_spec_in(:,:,66), pcols, lchnk)
    call outfld('FDLC_int66     ',lwdown_rad_spec_in(:,:,66), pcols, lchnk)
    call outfld('FULC_int67     ',lwup_rad_spec_in(:,:,67), pcols, lchnk)
    call outfld('FDLC_int67     ',lwdown_rad_spec_in(:,:,67), pcols, lchnk)
    call outfld('FULC_int68     ',lwup_rad_spec_in(:,:,68), pcols, lchnk)
    call outfld('FDLC_int68     ',lwdown_rad_spec_in(:,:,68), pcols, lchnk)

    call outfld('FULC_toa_int01     ',lwup_rad_spec_in(:,1,1), pcols, lchnk)
    call outfld('FDLC_toa_int01     ',lwdown_rad_spec_in(:,1,1), pcols, lchnk)
    call outfld('FULC_toa_int02     ',lwup_rad_spec_in(:,1,2), pcols, lchnk)
    call outfld('FDLC_toa_int02     ',lwdown_rad_spec_in(:,1,2), pcols, lchnk)
    call outfld('FULC_toa_int03     ',lwup_rad_spec_in(:,1,3), pcols, lchnk)
    call outfld('FDLC_toa_int03     ',lwdown_rad_spec_in(:,1,3), pcols, lchnk)
    call outfld('FULC_toa_int04     ',lwup_rad_spec_in(:,1,4), pcols, lchnk)
    call outfld('FDLC_toa_int04     ',lwdown_rad_spec_in(:,1,4), pcols, lchnk)
    call outfld('FULC_toa_int05     ',lwup_rad_spec_in(:,1,5), pcols, lchnk)
    call outfld('FDLC_toa_int05     ',lwdown_rad_spec_in(:,1,5), pcols, lchnk)
    call outfld('FULC_toa_int06     ',lwup_rad_spec_in(:,1,6), pcols, lchnk)
    call outfld('FDLC_toa_int06     ',lwdown_rad_spec_in(:,1,6), pcols, lchnk)
    call outfld('FULC_toa_int07     ',lwup_rad_spec_in(:,1,7), pcols, lchnk)
    call outfld('FDLC_toa_int07     ',lwdown_rad_spec_in(:,1,7), pcols, lchnk)
    call outfld('FULC_toa_int08     ',lwup_rad_spec_in(:,1,8), pcols, lchnk)
    call outfld('FDLC_toa_int08     ',lwdown_rad_spec_in(:,1,8), pcols, lchnk)
    call outfld('FULC_toa_int09     ',lwup_rad_spec_in(:,1,9), pcols, lchnk)
    call outfld('FDLC_toa_int09     ',lwdown_rad_spec_in(:,1,9), pcols, lchnk)
    call outfld('FULC_toa_int10     ',lwup_rad_spec_in(:,1,10), pcols, lchnk)
    call outfld('FDLC_toa_int10     ',lwdown_rad_spec_in(:,1,10), pcols, lchnk)
    call outfld('FULC_toa_int11     ',lwup_rad_spec_in(:,1,11), pcols, lchnk)
    call outfld('FDLC_toa_int11     ',lwdown_rad_spec_in(:,1,11), pcols, lchnk)
    call outfld('FULC_toa_int12     ',lwup_rad_spec_in(:,1,12), pcols, lchnk)
    call outfld('FDLC_toa_int12     ',lwdown_rad_spec_in(:,1,12), pcols, lchnk)
    call outfld('FULC_toa_int13     ',lwup_rad_spec_in(:,1,13), pcols, lchnk)
    call outfld('FDLC_toa_int13     ',lwdown_rad_spec_in(:,1,13), pcols, lchnk)
    call outfld('FULC_toa_int14     ',lwup_rad_spec_in(:,1,14), pcols, lchnk)
    call outfld('FDLC_toa_int14     ',lwdown_rad_spec_in(:,1,14), pcols, lchnk)
    call outfld('FULC_toa_int15     ',lwup_rad_spec_in(:,1,15), pcols, lchnk)
    call outfld('FDLC_toa_int15     ',lwdown_rad_spec_in(:,1,15), pcols, lchnk)
    call outfld('FULC_toa_int16     ',lwup_rad_spec_in(:,1,16), pcols, lchnk)
    call outfld('FDLC_toa_int16     ',lwdown_rad_spec_in(:,1,16), pcols, lchnk)
    call outfld('FULC_toa_int17     ',lwup_rad_spec_in(:,1,17), pcols, lchnk)
    call outfld('FDLC_toa_int17     ',lwdown_rad_spec_in(:,1,17), pcols, lchnk)
    call outfld('FULC_toa_int18     ',lwup_rad_spec_in(:,1,18), pcols, lchnk)
    call outfld('FDLC_toa_int18     ',lwdown_rad_spec_in(:,1,18), pcols, lchnk)
    call outfld('FULC_toa_int19     ',lwup_rad_spec_in(:,1,19), pcols, lchnk)
    call outfld('FDLC_toa_int19     ',lwdown_rad_spec_in(:,1,19), pcols, lchnk)
    call outfld('FULC_toa_int20     ',lwup_rad_spec_in(:,1,20), pcols, lchnk)
    call outfld('FDLC_toa_int20     ',lwdown_rad_spec_in(:,1,20), pcols, lchnk)
    call outfld('FULC_toa_int21     ',lwup_rad_spec_in(:,1,21), pcols, lchnk)
    call outfld('FDLC_toa_int21     ',lwdown_rad_spec_in(:,1,21), pcols, lchnk)
    call outfld('FULC_toa_int22     ',lwup_rad_spec_in(:,1,22), pcols, lchnk)
    call outfld('FDLC_toa_int22     ',lwdown_rad_spec_in(:,1,22), pcols, lchnk)
    call outfld('FULC_toa_int23     ',lwup_rad_spec_in(:,1,23), pcols, lchnk)
    call outfld('FDLC_toa_int23     ',lwdown_rad_spec_in(:,1,23), pcols, lchnk)
    call outfld('FULC_toa_int24     ',lwup_rad_spec_in(:,1,24), pcols, lchnk)
    call outfld('FDLC_toa_int24     ',lwdown_rad_spec_in(:,1,24), pcols, lchnk)
    call outfld('FULC_toa_int25     ',lwup_rad_spec_in(:,1,25), pcols, lchnk)
    call outfld('FDLC_toa_int25     ',lwdown_rad_spec_in(:,1,25), pcols, lchnk)
    call outfld('FULC_toa_int26     ',lwup_rad_spec_in(:,1,26), pcols, lchnk)
    call outfld('FDLC_toa_int26     ',lwdown_rad_spec_in(:,1,26), pcols, lchnk)
    call outfld('FULC_toa_int27     ',lwup_rad_spec_in(:,1,27), pcols, lchnk)
    call outfld('FDLC_toa_int27     ',lwdown_rad_spec_in(:,1,27), pcols, lchnk)
    call outfld('FULC_toa_int28     ',lwup_rad_spec_in(:,1,28), pcols, lchnk)
    call outfld('FDLC_toa_int28     ',lwdown_rad_spec_in(:,1,28), pcols, lchnk)
    call outfld('FULC_toa_int29     ',lwup_rad_spec_in(:,1,29), pcols, lchnk)
    call outfld('FDLC_toa_int29     ',lwdown_rad_spec_in(:,1,29), pcols, lchnk)
    call outfld('FULC_toa_int30     ',lwup_rad_spec_in(:,1,30), pcols, lchnk)
    call outfld('FDLC_toa_int30     ',lwdown_rad_spec_in(:,1,30), pcols, lchnk)
    call outfld('FULC_toa_int31     ',lwup_rad_spec_in(:,1,31), pcols, lchnk)
    call outfld('FDLC_toa_int31     ',lwdown_rad_spec_in(:,1,31), pcols, lchnk)
    call outfld('FULC_toa_int32     ',lwup_rad_spec_in(:,1,32), pcols, lchnk)
    call outfld('FDLC_toa_int32     ',lwdown_rad_spec_in(:,1,32), pcols, lchnk)
    call outfld('FULC_toa_int33     ',lwup_rad_spec_in(:,1,33), pcols, lchnk)
    call outfld('FDLC_toa_int33     ',lwdown_rad_spec_in(:,1,33), pcols, lchnk)
    call outfld('FULC_toa_int34     ',lwup_rad_spec_in(:,1,34), pcols, lchnk)
    call outfld('FDLC_toa_int34     ',lwdown_rad_spec_in(:,1,34), pcols, lchnk)
    call outfld('FULC_toa_int35     ',lwup_rad_spec_in(:,1,35), pcols, lchnk)
    call outfld('FDLC_toa_int35     ',lwdown_rad_spec_in(:,1,35), pcols, lchnk)
    call outfld('FULC_toa_int36     ',lwup_rad_spec_in(:,1,36), pcols, lchnk)
    call outfld('FDLC_toa_int36     ',lwdown_rad_spec_in(:,1,36), pcols, lchnk)
    call outfld('FULC_toa_int37     ',lwup_rad_spec_in(:,1,37), pcols, lchnk)
    call outfld('FDLC_toa_int37     ',lwdown_rad_spec_in(:,1,37), pcols, lchnk)
    call outfld('FULC_toa_int38     ',lwup_rad_spec_in(:,1,38), pcols, lchnk)
    call outfld('FDLC_toa_int38     ',lwdown_rad_spec_in(:,1,38), pcols, lchnk)
    call outfld('FULC_toa_int39     ',lwup_rad_spec_in(:,1,39), pcols, lchnk)
    call outfld('FDLC_toa_int39     ',lwdown_rad_spec_in(:,1,39), pcols, lchnk)
    call outfld('FULC_toa_int40     ',lwup_rad_spec_in(:,1,40), pcols, lchnk)
    call outfld('FDLC_toa_int40     ',lwdown_rad_spec_in(:,1,40), pcols, lchnk)
    call outfld('FULC_toa_int41     ',lwup_rad_spec_in(:,1,41), pcols, lchnk)
    call outfld('FDLC_toa_int41     ',lwdown_rad_spec_in(:,1,41), pcols, lchnk)
    call outfld('FULC_toa_int42     ',lwup_rad_spec_in(:,1,42), pcols, lchnk)
    call outfld('FDLC_toa_int42     ',lwdown_rad_spec_in(:,1,42), pcols, lchnk)
    call outfld('FULC_toa_int43     ',lwup_rad_spec_in(:,1,43), pcols, lchnk)
    call outfld('FDLC_toa_int43     ',lwdown_rad_spec_in(:,1,43), pcols, lchnk)
    call outfld('FULC_toa_int44     ',lwup_rad_spec_in(:,1,44), pcols, lchnk)
    call outfld('FDLC_toa_int44     ',lwdown_rad_spec_in(:,1,44), pcols, lchnk)
    call outfld('FULC_toa_int45     ',lwup_rad_spec_in(:,1,45), pcols, lchnk)
    call outfld('FDLC_toa_int45     ',lwdown_rad_spec_in(:,1,45), pcols, lchnk)
    call outfld('FULC_toa_int46     ',lwup_rad_spec_in(:,1,46), pcols, lchnk)
    call outfld('FDLC_toa_int46     ',lwdown_rad_spec_in(:,1,46), pcols, lchnk)
    call outfld('FULC_toa_int47     ',lwup_rad_spec_in(:,1,47), pcols, lchnk)
    call outfld('FDLC_toa_int47     ',lwdown_rad_spec_in(:,1,47), pcols, lchnk)
    call outfld('FULC_toa_int48     ',lwup_rad_spec_in(:,1,48), pcols, lchnk)
    call outfld('FDLC_toa_int48     ',lwdown_rad_spec_in(:,1,48), pcols, lchnk)
    call outfld('FULC_toa_int49     ',lwup_rad_spec_in(:,1,49), pcols, lchnk)
    call outfld('FDLC_toa_int49     ',lwdown_rad_spec_in(:,1,49), pcols, lchnk)
    call outfld('FULC_toa_int50     ',lwup_rad_spec_in(:,1,50), pcols, lchnk)
    call outfld('FDLC_toa_int50     ',lwdown_rad_spec_in(:,1,50), pcols, lchnk)
    call outfld('FULC_toa_int51     ',lwup_rad_spec_in(:,1,51), pcols, lchnk)
    call outfld('FDLC_toa_int51     ',lwdown_rad_spec_in(:,1,51), pcols, lchnk)
    call outfld('FULC_toa_int52     ',lwup_rad_spec_in(:,1,52), pcols, lchnk)
    call outfld('FDLC_toa_int52     ',lwdown_rad_spec_in(:,1,52), pcols, lchnk)
    call outfld('FULC_toa_int53     ',lwup_rad_spec_in(:,1,53), pcols, lchnk)
    call outfld('FDLC_toa_int53     ',lwdown_rad_spec_in(:,1,53), pcols, lchnk)
    call outfld('FULC_toa_int54     ',lwup_rad_spec_in(:,1,54), pcols, lchnk)
    call outfld('FDLC_toa_int54     ',lwdown_rad_spec_in(:,1,54), pcols, lchnk)
    call outfld('FULC_toa_int55     ',lwup_rad_spec_in(:,1,55), pcols, lchnk)
    call outfld('FDLC_toa_int55     ',lwdown_rad_spec_in(:,1,55), pcols, lchnk)
    call outfld('FULC_toa_int56     ',lwup_rad_spec_in(:,1,56), pcols, lchnk)
    call outfld('FDLC_toa_int56     ',lwdown_rad_spec_in(:,1,56), pcols, lchnk)
    call outfld('FULC_toa_int57     ',lwup_rad_spec_in(:,1,57), pcols, lchnk)
    call outfld('FDLC_toa_int57     ',lwdown_rad_spec_in(:,1,57), pcols, lchnk)
    call outfld('FULC_toa_int58     ',lwup_rad_spec_in(:,1,58), pcols, lchnk)
    call outfld('FDLC_toa_int58     ',lwdown_rad_spec_in(:,1,58), pcols, lchnk)
    call outfld('FULC_toa_int59     ',lwup_rad_spec_in(:,1,59), pcols, lchnk)
    call outfld('FDLC_toa_int59     ',lwdown_rad_spec_in(:,1,59), pcols, lchnk)
    call outfld('FULC_toa_int60     ',lwup_rad_spec_in(:,1,60), pcols, lchnk)
    call outfld('FDLC_toa_int60     ',lwdown_rad_spec_in(:,1,60), pcols, lchnk)
    call outfld('FULC_toa_int61     ',lwup_rad_spec_in(:,1,61), pcols, lchnk)
    call outfld('FDLC_toa_int61     ',lwdown_rad_spec_in(:,1,61), pcols, lchnk)
    call outfld('FULC_toa_int62     ',lwup_rad_spec_in(:,1,62), pcols, lchnk)
    call outfld('FDLC_toa_int62     ',lwdown_rad_spec_in(:,1,62), pcols, lchnk)
    call outfld('FULC_toa_int63     ',lwup_rad_spec_in(:,1,63), pcols, lchnk)
    call outfld('FDLC_toa_int63     ',lwdown_rad_spec_in(:,1,63), pcols, lchnk)
    call outfld('FULC_toa_int64     ',lwup_rad_spec_in(:,1,64), pcols, lchnk)
    call outfld('FDLC_toa_int64     ',lwdown_rad_spec_in(:,1,64), pcols, lchnk)
    call outfld('FULC_toa_int65     ',lwup_rad_spec_in(:,1,65), pcols, lchnk)
    call outfld('FDLC_toa_int65     ',lwdown_rad_spec_in(:,1,65), pcols, lchnk)
    call outfld('FULC_toa_int66     ',lwup_rad_spec_in(:,1,66), pcols, lchnk)
    call outfld('FDLC_toa_int66     ',lwdown_rad_spec_in(:,1,66), pcols, lchnk)
    call outfld('FULC_toa_int67     ',lwup_rad_spec_in(:,1,67), pcols, lchnk)
    call outfld('FDLC_toa_int67     ',lwdown_rad_spec_in(:,1,67), pcols, lchnk)
    call outfld('FULC_toa_int68     ',lwup_rad_spec_in(:,1,68), pcols, lchnk)
    call outfld('FDLC_toa_int68     ',lwdown_rad_spec_in(:,1,68), pcols, lchnk)

end subroutine outfld_spectral_flux_clearsky
end subroutine radiation_tend

!===============================================================================

subroutine radiation_output_sw(lchnk, ncol, icall, rd, pbuf, cam_out)

   ! Dump shortwave radiation information to history buffer.

   integer ,               intent(in) :: lchnk
   integer,                intent(in) :: ncol
   integer,                intent(in) :: icall
   type(rad_out_t),        intent(in) :: rd
   type(physics_buffer_desc), pointer :: pbuf(:)
   type(cam_out_t),        intent(in) :: cam_out

   ! local variables
   real(r8), pointer :: qrs(:,:)
   real(r8), pointer :: fsnt(:)
   real(r8), pointer :: fsns(:)
   real(r8), pointer :: fsds(:)

   real(r8) :: ftem(pcols)
   !----------------------------------------------------------------------------

   call pbuf_get_field(pbuf, qrs_idx,  qrs)
   call pbuf_get_field(pbuf, fsnt_idx, fsnt)
   call pbuf_get_field(pbuf, fsns_idx, fsns)
   call pbuf_get_field(pbuf, fsds_idx, fsds)

!!$   call outfld('SOLIN'//diag(icall),    rd%solin,      pcols, lchnk)
!!$
!!$   call outfld('QRS'//diag(icall),      qrs(:ncol,:)/cpair,     ncol, lchnk)
!!$   call outfld('QRSC'//diag(icall),     rd%qrsc(:ncol,:)/cpair, ncol, lchnk)
!!$
!!$   call outfld('FSNT'//diag(icall),     fsnt,          pcols, lchnk)
!!$   call outfld('FSNTC'//diag(icall),    rd%fsntc,      pcols, lchnk)
!!$   call outfld('FSNTOA'//diag(icall),   rd%fsntoa,     pcols, lchnk)
!!$   call outfld('FSNTOAC'//diag(icall),  rd%fsntoac,    pcols, lchnk)
!!$
!!$   ftem(:ncol) = rd%fsntoa(:ncol) - rd%fsntoac(:ncol)
!!$   call outfld('SWCF'//diag(icall),     ftem,          pcols, lchnk)
!!$
!!$   call outfld('FSUTOA'//diag(icall),   rd%fsutoa,     pcols, lchnk)
!!$
!!$   call outfld('FSNIRTOA'//diag(icall), rd%fsnirt,     pcols, lchnk)
!!$   call outfld('FSNRTOAC'//diag(icall), rd%fsnrtc,     pcols, lchnk)
!!$   call outfld('FSNRTOAS'//diag(icall), rd%fsnirtsq,   pcols, lchnk)
!!$
!!$   call outfld('FSN200'//diag(icall),   rd%fsn200,     pcols, lchnk)
!!$   call outfld('FSN200C'//diag(icall),  rd%fsn200c,    pcols, lchnk)
!!$
!!$   call outfld('FSNR'//diag(icall),     rd%fsnr,       pcols, lchnk)
!!$
!!$   call outfld('SOLS'//diag(icall),     cam_out%sols,  pcols, lchnk)
!!$   call outfld('SOLL'//diag(icall),     cam_out%soll,  pcols, lchnk)
!!$   call outfld('SOLSD'//diag(icall),    cam_out%solsd, pcols, lchnk)
!!$   call outfld('SOLLD'//diag(icall),    cam_out%solld, pcols, lchnk)
!!$
!!$   call outfld('FSNS'//diag(icall),     fsns,          pcols, lchnk)
!!$   call outfld('FSNSC'//diag(icall),    rd%fsnsc,      pcols, lchnk)
!!$
!!$   call outfld('FSDS'//diag(icall),     fsds,          pcols, lchnk)
!!$   call outfld('FSDSC'//diag(icall),    rd%fsdsc,      pcols, lchnk)

end subroutine radiation_output_sw


!===============================================================================

subroutine radiation_output_cld(lchnk, ncol, rd)

   ! Dump shortwave cloud optics information to history buffer.

   integer ,        intent(in) :: lchnk
   integer,         intent(in) :: ncol
   type(rad_out_t), intent(in) :: rd
   !----------------------------------------------------------------------------

!!$   call outfld('TOT_CLD_VISTAU',  rd%tot_cld_vistau,  pcols, lchnk)
!!$   call outfld('TOT_ICLD_VISTAU', rd%tot_icld_vistau, pcols, lchnk)
!!$   call outfld('LIQ_ICLD_VISTAU', rd%liq_icld_vistau, pcols, lchnk)
!!$   call outfld('ICE_ICLD_VISTAU', rd%ice_icld_vistau, pcols, lchnk)
!!$   if (cldfsnow_idx > 0) then
!!$      call outfld('SNOW_ICLD_VISTAU', rd%snow_icld_vistau, pcols, lchnk)
!!$   endif
end subroutine radiation_output_cld

!===============================================================================

subroutine radiation_output_lw(lchnk, ncol, icall, rd, pbuf, cam_out, freqclr, flntclr)

   ! Dump longwave radiation information to history buffer

   integer,                intent(in) :: lchnk
   integer,                intent(in) :: ncol
   integer,                intent(in) :: icall  ! icall=0 for climate diagnostics
   type(rad_out_t),        intent(in) :: rd
   type(physics_buffer_desc), pointer :: pbuf(:)
   type(cam_out_t),        intent(in) :: cam_out
   real(r8),               intent(in) :: freqclr(pcols)
   real(r8),               intent(in) :: flntclr(pcols)

   ! local variables
   real(r8), pointer :: qrl(:,:)
   real(r8), pointer :: flnt(:)
   real(r8), pointer :: flns(:)

   real(r8) :: ftem(pcols)
   !----------------------------------------------------------------------------

   call pbuf_get_field(pbuf, qrl_idx,  qrl)
   call pbuf_get_field(pbuf, flnt_idx, flnt)
   call pbuf_get_field(pbuf, flns_idx, flns)

   call outfld('QRL'//diag(icall),     qrl(:ncol,:)/cpair,     ncol, lchnk)
!!$   call outfld('QRLC'//diag(icall),    rd%qrlc(:ncol,:)/cpair, ncol, lchnk)
!!$
!!$   call outfld('FLNT'//diag(icall),    flnt,          pcols, lchnk)
!!$   call outfld('FLNTC'//diag(icall),   rd%flntc,      pcols, lchnk)
!!$
!!$   call outfld('FREQCLR'//diag(icall), freqclr,       pcols, lchnk)
!!$   call outfld('FLNTCLR'//diag(icall), flntclr,       pcols, lchnk)
!!$
!!$   call outfld('FLUT'//diag(icall),    rd%flut,       pcols, lchnk)
!!$   call outfld('FLUTC'//diag(icall),   rd%flutc,      pcols, lchnk)
!!$
!!$   ftem(:ncol) = rd%flutc(:ncol) - rd%flut(:ncol)
!!$   call outfld('LWCF'//diag(icall),    ftem,          pcols, lchnk)
!!$
!!$   call outfld('FLN200'//diag(icall),  rd%fln200,     pcols, lchnk)
!!$   call outfld('FLN200C'//diag(icall), rd%fln200c,    pcols, lchnk)
!!$
!!$   call outfld('FLNR'//diag(icall),    rd%flnr,       pcols, lchnk)
!!$
!!$   call outfld('FLNS'//diag(icall),    flns,          pcols, lchnk)
!!$   call outfld('FLNSC'//diag(icall),   rd%flnsc,      pcols, lchnk)
!!$
!!$   call outfld('FLDS'//diag(icall),    cam_out%flwds, pcols, lchnk)
!!$   call outfld('FLDSC'//diag(icall),   rd%fldsc,      pcols, lchnk)

end subroutine radiation_output_lw

!===============================================================================

subroutine calc_col_mean(state, mmr_pointer, mean_value)

   ! Compute the column mean mass mixing ratio.

   type(physics_state),        intent(in)  :: state
   real(r8), dimension(:,:),   pointer     :: mmr_pointer  ! mass mixing ratio (lev)
   real(r8), dimension(pcols), intent(out) :: mean_value   ! column mean mmr

   integer  :: i, k, ncol
   real(r8) :: ptot(pcols)
   !-----------------------------------------------------------------------

   ncol         = state%ncol
   mean_value   = 0.0_r8
   ptot         = 0.0_r8

   do k=1,pver
      do i=1,ncol
         mean_value(i) = mean_value(i) + mmr_pointer(i,k)*state%pdeldry(i,k)
         ptot(i)         = ptot(i) + state%pdeldry(i,k)
      end do
   end do
   do i=1,ncol
      mean_value(i) = mean_value(i) / ptot(i)
   end do

end subroutine calc_col_mean

subroutine init_kcoeff

!------------------------------------------------------------------------
!
! Purpose:  Initialize k coefficient data from input file.
!
!------------------------------------------------------------------------

#if ( defined SPMD)
  use mpishorthand
#endif

    use ioFileMod, only: getfil
    use cam_pio_utils, only: cam_pio_openfile


    implicit none
    include 'netcdf.inc'

!------------------------------------------------------------------------
!
! Local Variables
!
    type(file_desc_t) :: ncid
    integer :: gid
    integer :: pid
    integer :: tid
    integer :: wid
    integer :: nid
    integer :: keff_id
    character(len=256) :: locfn
    character(len=256) :: filename
    integer :: ierr

!------------------------------------------------------------------------
!
! Start Code
!

    if ( masterproc ) then
      write (6, '(2x, a)') '_______________________________________________________'
      write (6, '(2x, a)') '_________ initializing gas absorption coeffs __________'
      write (6, '(2x, a)') '_______________________________________________________'
    endif

    ! Load K coefficients
    call getfil(trim(k_h2o_file), locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, k_h2o)
    call pio_closefile(ncid)

    call getfil(trim(k_co2_file), locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, k_co2)
    call pio_closefile(ncid)

    call getfil(trim(k_ch4_file), locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, k_ch4)
    call pio_closefile(ncid)

    call getfil(trim(k_c2h6_file), locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, k_c2h6)
    call pio_closefile(ncid)

    call getfil(trim(k_o3_file), locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, k_o3)
    call pio_closefile(ncid)

    call getfil(trim(k_o2_file), locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, k_o2)
    call pio_closefile(ncid)

    !! Load mtckd h2o continuum
    call getfil(trim(kh2o_mtckd_file), locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'KSELF',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, kh2oself_mtckd)
    call pio_closefile(ncid)

    call getfil(trim(kh2o_mtckd_file), locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'KFRGN',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, kh2ofrgn_mtckd)
    call pio_closefile(ncid)
    !! mtckd

    ! Load K coefficients, for n2n2 continuum
    call getfil(trim(kn2n2cia_file ), locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'sigma',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, kn2n2)
    call pio_closefile(ncid)

    ! Load K coefficients, for n2h2 continuum
    call getfil(trim(kn2h2cia_file ), locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'sigma',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, kn2h2)
    call pio_closefile(ncid)

    ! Load K coefficients, for h2h2 continuum
    call getfil(trim(kh2h2cia_file ), locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'sigma',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, kh2h2)
    call pio_closefile(ncid)

    ! Load K coefficients, for co2co2 lw continuum
    call getfil(trim(kco2co2cia_lw_file ), locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'sigma',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, kco2co2_lw)
    call pio_closefile(ncid)

    ! Load K coefficients, for co2co2 sw continuum
    call getfil(trim(kco2co2cia_sw_file ), locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'sigma',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, kco2co2_sw)
    call pio_closefile(ncid)

    ! Load K coefficients, for co2h2 continuum
    call getfil(trim(kco2h2cia_file ), locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'sigma',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, kco2h2)
    call pio_closefile(ncid)

    ! Load K coefficients, for co2h4 continuum
    call getfil(trim(kco2ch4cia_file ), locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'sigma',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, kco2ch4)
    call pio_closefile(ncid)

    call getfil(ko2o2cia_file, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'sigma',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, ko2o2)
    call pio_closefile(ncid)

    call getfil(ko2n2cia_file, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'sigma',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, ko2n2)
    call pio_closefile(ncid)

    call getfil(ko2co2cia_file, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'sigma',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, ko2co2)
    call pio_closefile(ncid)



! broadcast optical constants to all nodes
#if ( defined SPMD )
    call mpibcast(k_h2o,  ntot_wavlnrng*ngauss_8gpt*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k_co2,  ntot_wavlnrng*ngauss_8gpt*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k_ch4,  ntot_wavlnrng*ngauss_8gpt*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k_c2h6, ntot_wavlnrng*ngauss_8gpt*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k_o2,   ntot_wavlnrng*ngauss_8gpt*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k_o3,   ntot_wavlnrng*ngauss_8gpt*kc_npress*kc_ntemp, mpir8, 0, mpicom)

    call mpibcast(kh2oself_mtckd, ngauss_8gpt*ntot_wavlnrng*kmtckd_ntemp, mpir8, 0, mpicom)
    call mpibcast(kh2ofrgn_mtckd, ngauss_8gpt*ntot_wavlnrng*kmtckd_ntemp, mpir8, 0, mpicom)

    call mpibcast(kn2n2, ntot_wavlnrng*kn2n2_ntemp, mpir8, 0, mpicom)
    call mpibcast(kn2h2, ntot_wavlnrng*kn2h2_ntemp, mpir8, 0, mpicom)
    call mpibcast(kh2h2, ntot_wavlnrng*kh2h2_ntemp, mpir8, 0, mpicom)

    call mpibcast(kco2co2_sw, ntot_wavlnrng*kco2co2_sw_ntemp, mpir8, 0, mpicom)
    call mpibcast(kco2co2_lw, ntot_wavlnrng*kco2co2_lw_ntemp, mpir8, 0, mpicom)
    call mpibcast(kco2h2, ntot_wavlnrng*kco2h2_ntemp, mpir8, 0, mpicom)
    call mpibcast(kco2ch4, ntot_wavlnrng*kco2ch4_ntemp, mpir8, 0, mpicom)

    call mpibcast(ko2o2,  ntot_wavlnrng*ko2o2_ntemp,  mpir8, 0, mpicom)
    call mpibcast(ko2n2,  ntot_wavlnrng*ko2n2_ntemp,  mpir8, 0, mpicom)
    call mpibcast(ko2co2, ntot_wavlnrng*ko2co2_ntemp, mpir8, 0, mpicom)


#endif


  end subroutine init_kcoeff


!============================================================================

  subroutine init_solar

!------------------------------------------------------------------------
!
! Purpose:  Initialize solar data from input file.
!
!------------------------------------------------------------------------
!
#if ( defined SPMD)
  use mpishorthand
#endif

    use ioFileMod, only: getfil
    use cam_pio_utils, only: cam_pio_openfile
    use pio,  only: pio_inq_varid, pio_get_var, pio_closefile, pio_nowrite,  &
                    file_desc_t, var_desc_t, pio_inq_dimid, pio_inquire_dimension
    implicit none
    include 'netcdf.inc'

!------------------------------------------------------------------------
!
! Local Variables
!
    type(file_desc_t) :: ncid
    character(len=256) :: locfn
    integer :: solarflux_id
    integer :: S0_id
    integer :: ierr

    if (masterproc) then
        write(6,*) "INITIALIZING SOLAR SPECTRAL FILE"
    endif

    ! Load solar data
    call getfil(solar_file, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'S0',   S0_id)
    ierr =  pio_get_var(ncid, S0_id, S0)
    ierr =  pio_inq_varid(ncid, 'solarflux',   solarflux_id)
    ierr =  pio_get_var(ncid, solarflux_id, solarflux)
    call pio_closefile(ncid)


#if ( defined SPMD )
    call mpibcast(S0, 1, mpir8, 0, mpicom)
    call mpibcast(solarflux, ntot_wavlnrng, mpir8, 0, mpicom)
#endif

  end subroutine init_solar


!============================================================================

  subroutine init_cldopts

!------------------------------------------------------------------------
!
! Purpose:  Initialize the cloud optical constants from input file.
!
!------------------------------------------------------------------------

#if ( defined SPMD)
  use mpishorthand
#endif

    use ioFileMod, only: getfil
    use cam_pio_utils, only: cam_pio_openfile
    use pio,  only: pio_inq_varid, pio_get_var, pio_closefile, pio_nowrite,  &
                    file_desc_t, var_desc_t, pio_inq_dimid, pio_inquire_dimension
    implicit none
    include 'netcdf.inc'

!------------------------------------------------------------------------
!
! Local Variables

    type(file_desc_t) :: ncid
    integer :: bin_id
    integer :: wav_id
    integer :: ncldopt_lbins
    integer :: ncldopt_lwavs
    integer :: ncldopt_ibins
    integer :: ncldopt_iwavs
    integer :: q_id
    integer :: w_id
    integer :: g_id
    character(len=256) :: filename
    character(len=256) :: locfn
    integer :: ierr

!------------------------------------------------------------------------
!
! Start Code
!
    if (masterproc) then
      write(6,*) "CLDOPTS: INITIALIZING CLOUD OPTICAL PROPERTIES"
    endif

    ! Load K water cloud optics file
    call getfil(trim(cldoptsL_file), locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
!!    ierr =  pio_inq_dimid(ncid, 'rel_bins',   bin_id)
!!    ierr =  pio_inquire_dimension(ncid, bin_id, len=ncldopt_lbins)
    ierr =  pio_inq_varid(ncid, 'Qext_liq',   q_id)
    ierr =  pio_get_var(ncid, q_id, Qcldliq)
    ierr =  pio_inq_varid(ncid, 'W_liq',   w_id)
    ierr =  pio_get_var(ncid, w_id, Wcldliq)
    ierr =  pio_inq_varid(ncid, 'G_liq',   g_id)
    ierr =  pio_get_var(ncid, g_id, Gcldliq)
    call pio_closefile(ncid)


    ! Load ice cloud optics file
    call getfil(trim(cldoptsI_file), locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
!!    ierr =  pio_inq_dimid(ncid, 'rei_bins',   bin_id)
!!    ierr =  pio_inquire_dimension(ncid, bin_id, len=ncldopt_ibins)
    ierr =  pio_inq_varid(ncid, 'Qext_ice',   q_id)
    ierr =  pio_get_var(ncid, q_id, Qcldice)
    ierr =  pio_inq_varid(ncid, 'W_ice',   w_id)
    ierr =  pio_get_var(ncid, w_id, Wcldice)
    ierr =  pio_inq_varid(ncid, 'G_ice',   g_id)
    ierr =  pio_get_var(ncid, g_id, Gcldice)
    call pio_closefile(ncid)

! broadcast water cloud optical constants to all nodes
#if ( defined SPMD )
    call mpibcast(Qcldliq, nrel*ntot_wavlnrng, mpir8, 0, mpicom)
    call mpibcast(Wcldliq, nrel*ntot_wavlnrng, mpir8, 0, mpicom)
    call mpibcast(Gcldliq, nrel*ntot_wavlnrng, mpir8, 0, mpicom)
#endif

! broadcast ice cloud optical constants to all nodes
#if ( defined SPMD )
    call mpibcast(Qcldice, nrei*ntot_wavlnrng, mpir8, 0, mpicom)
    call mpibcast(Wcldice, nrei*ntot_wavlnrng, mpir8, 0, mpicom)
    call mpibcast(Gcldice, nrei*ntot_wavlnrng, mpir8, 0, mpicom)
#endif

  end subroutine init_cldopts



end module radiation
