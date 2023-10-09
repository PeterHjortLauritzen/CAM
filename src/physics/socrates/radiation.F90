module radiation

   ! stub module

   use shr_kind_mod,        only: r8=>shr_kind_r8, cl=>shr_kind_cl
   use ppgrid,              only: pcols, pver
   use camsrfexch,          only: cam_in_t,cam_out_t
   use physics_types,       only: physics_state, physics_tend, physics_ptend
   use physics_buffer,      only: physics_buffer_desc, pbuf_add_field, pbuf_get_field, dtype_r8
   use cam_history,         only: addfld, add_default, horiz_only, outfld, hist_fld_active
   use pio,                 only: file_desc_t, var_desc_t,               &
                               pio_int, pio_double, pio_noerr,        &
                               pio_seterrorhandling, pio_bcast_error, &
                               pio_inq_varid, pio_def_var,            &
                               pio_put_var, pio_get_var, pio_put_att


   implicit none
   private
   save

   public :: &
      radiation_do,             &
      radiation_readnl,         &
      radiation_init,           &
      radiation_register,       & 
      radiation_tend,           &
      radiation_define_restart, &! define variables for restart
      radiation_write_restart,  &! write variables to restart
      radiation_read_restart,   &! read variables from restart
      radheat_tend              ! <-- does not need to be public

   real(r8), public, protected :: nextsw_cday = -1._r8 ! future radiation calday for surface models

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

   ! physics buffer indices
   integer :: qrs_idx  = 0
   integer :: qrl_idx    = 0
!========================================================================================
contains
!========================================================================================
function radiation_do(op, timestep)

   ! Returns true if the specified operation is done this timestep.

   character(len=*), intent(in) :: op             ! name of operation
   integer, intent(in), optional:: timestep
   logical                      :: radiation_do   ! return value
   !---------------------------------------------------------------------------

   radiation_do = .false.

end function radiation_do


subroutine radiation_readnl(nlfile)

   ! this stub can be called, but does nothing

   character(len=*), intent(in) :: nlfile

end subroutine radiation_readnl

!========================================================================================

subroutine radiation_register()

   ! Register radiation fields in the physics buffer

   call pbuf_add_field('QRS' , 'global',dtype_r8,(/pcols,pver/), qrs_idx) ! shortwave radiative heating rate
   call pbuf_add_field('QRL' , 'global',dtype_r8,(/pcols,pver/), qrl_idx) ! longwave  radiative heating rate

   ! call pbuf_add_field('FSDS' , 'global',dtype_r8,(/pcols/), fsds_idx) ! Surface solar downward flux
   ! call pbuf_add_field('FSNS' , 'global',dtype_r8,(/pcols/), fsns_idx) ! Surface net shortwave flux
   ! call pbuf_add_field('FSNT' , 'global',dtype_r8,(/pcols/), fsnt_idx) ! Top-of-model net shortwave flux

   ! call pbuf_add_field('FLNS' , 'global',dtype_r8,(/pcols/), flns_idx) ! Surface net longwave flux
   ! call pbuf_add_field('FLNT' , 'global',dtype_r8,(/pcols/), flnt_idx) ! Top-of-model net longwave flux

end subroutine radiation_register

!========================================================================================

subroutine radiation_init(pbuf2d)

   ! arguments
   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   ! 
   ! initialize any additional radiation modules needed here
   !

   !
   ! provide output fields
   !
   call addfld('QRS',  (/ 'lev' /), 'A', 'K/s', 'Solar heating rate', sampling_seq='rad_lwsw')
   call addfld('QRL',  (/ 'lev' /), 'A', 'K/s', 'Longwave heating rate', sampling_seq='rad_lwsw')

   call add_default('QRS',   1, ' ')                                                  
   call add_default('QRL',   1, ' ')                                                  
end subroutine radiation_init

!========================================================================================

subroutine radiation_tend( &
   state, ptend, pbuf, cam_out, cam_in, net_flx, rd_out)

   ! Returns true if the specified operation is done this timestep.

   ! arguments
   type(physics_state), intent(in)    :: state
   type(physics_ptend), intent(out)   :: ptend            ! Package tendencies
   type(physics_buffer_desc), pointer :: pbuf(:)
   type(cam_in_t),      intent(in)    :: cam_in
   type(cam_out_t),     intent(inout) :: cam_out
   real(r8),            intent(out)        :: net_flx(pcols)
   type(rad_out_t), target, optional, intent(out) :: rd_out

   ! character(len=*), intent(in) :: op             ! name of operation
   ! integer, intent(in), optional:: timestep
   integer           :: lchnk         ! chunk identifier
   integer           :: ncol          ! number of atmospheric columns
   real(r8), pointer :: qrs(:,:)      ! shortwave radiative heating rate
   real(r8), pointer :: qrl(:,:)      ! longwave  radiative heating rate


   lchnk = state%lchnk
   ncol  = state%ncol

   call pbuf_get_field(pbuf, qrs_idx, qrs)
   call pbuf_get_field(pbuf, qrl_idx, qrl)

   !
   ! TODO: Get optical properties
   !

   !
   ! TODO: Radiation Calculation
   !
   qrs(:,:) = -1.0_r8
   qrl(:,:) =  1.0_r8

   call outfld('QRS', qrs(:ncol,:), ncol, lchnk)
   call outfld('QRL', qrl(:ncol,:), ncol, lchnk)

   ! apply to tendency
   call radheat_tend(state, pbuf, ptend, qrl, qrs)
   
end subroutine radiation_tend

!========================================================================================

subroutine radheat_tend(state, pbuf, ptend, qrl, qrs)
   use physics_types,      only: physics_ptend_init
   
   !-----------------------------------------------------------------------
   ! Compute net radiative heating from qrs and qrl, and the associated net
   ! boundary flux.
   ! See physics/cam/radheat.F90 for full version
   !-----------------------------------------------------------------------

   ! Arguments
   type(physics_state), intent(in)  :: state             ! Physics state variables
   type(physics_buffer_desc), pointer :: pbuf(:)
   type(physics_ptend), intent(out) :: ptend             ! indivdual parameterization tendencie
   real(r8),            intent(in)  :: qrl(pcols,pver)   ! longwave heating, (K/s)*cpair
   real(r8),            intent(in)  :: qrs(pcols,pver)   ! shortwave heating, (K/s)*cpair

   integer :: ncol

   ncol = state%ncol
   call physics_ptend_init(ptend,state%psetcols, 'radheat', ls=.true.)
   ptend%s(:ncol,:) = (qrs(:ncol,:) + qrl(:ncol,:))

end subroutine radheat_tend

!===============================================================================

subroutine radiation_define_restart(file)

   ! define variables to be written to restart file

   ! arguments
   type(file_desc_t), intent(inout) :: file

   ! local variables
   integer :: ierr
   !----------------------------------------------------------------------------

!!$   call pio_seterrorhandling(File, PIO_BCAST_ERROR)
!!$
!!$   ierr = pio_def_var(File, 'nextsw_cday', pio_double, nextsw_cday_desc)
!!$   ierr = pio_put_att(File, nextsw_cday_desc, 'long_name', 'future radiation calday for surface models')
!!$   if (docosp) then
!!$      ierr = pio_def_var(File, 'cosp_cnt_init', pio_int, cospcnt_desc)
!!$   end if

end subroutine radiation_define_restart

!===============================================================================

subroutine radiation_write_restart(file)

   ! write variables to restart file

   ! arguments
   type(file_desc_t), intent(inout) :: file

   ! local variables
   integer :: ierr
   !----------------------------------------------------------------------------

!!$   ierr = pio_put_var(File, nextsw_cday_desc, (/ nextsw_cday /))
!!$   if (docosp) then
!!$      ierr = pio_put_var(File, cospcnt_desc, (/cosp_cnt(begchunk)/))
!!$   end if

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

!!$   if (docosp) then
!!$      call pio_seterrorhandling(File, PIO_BCAST_ERROR, err_handling)
!!$      ierr = pio_inq_varid(File, 'cosp_cnt_init', vardesc)
!!$      call pio_seterrorhandling(File, err_handling)
!!$      if (ierr /= PIO_NOERR) then
!!$         cosp_cnt_init = 0
!!$      else
!!$         ierr = pio_get_var(File, vardesc, cosp_cnt_init)
!!$      end if
!!$   end if
!!$
!!$   ierr = pio_inq_varid(File, 'nextsw_cday', vardesc)
!!$   ierr = pio_get_var(File, vardesc, temp_var)
!!$   nextsw_cday = temp_var

end subroutine radiation_read_restart

!===============================================================================

end module radiation

