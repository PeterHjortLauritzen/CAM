module radiation

! stub module

use shr_kind_mod,        only: r8=>shr_kind_r8, cl=>shr_kind_cl
use ppgrid,              only: pcols, pver
use physics_types,      only: physics_state, physics_ptend
use physics_buffer,     only: physics_buffer_desc
use camsrfexch,         only: cam_out_t, cam_in_t


implicit none
private
save

public :: &
   radiation_readnl,         &
   radiation_register,         &
   radiation_tend,         &
   rad_is_active,            &
   radiation_init,           &
   radiation_do

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

real(r8), public, protected :: nextsw_cday = -1._r8 ! future radiation calday for surface models

!========================================================================================
contains
!========================================================================================

function rad_is_active()
  !-----------------------------------------------------------------------
  logical :: rad_is_active
  !-----------------------------------------------------------------------
  rad_is_active = .false.
end function rad_is_active

!================================================================================================

subroutine radiation_readnl(nlfile)

   ! this stub can be called, but does nothing

   character(len=*), intent(in) :: nlfile

end subroutine radiation_readnl

!========================================================================================

subroutine radiation_init(pbuf2d)
   ! arguments
   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   ! this stub can be called, but does nothing

end subroutine radiation_init
!================================================================================================

subroutine radiation_register

   ! this stub can be called, but does nothing

end subroutine radiation_register

!========================================================================================

function radiation_do(op, timestep)

   ! Returns true if the specified operation is done this timestep.

   character(len=*), intent(in) :: op             ! name of operation
   integer, intent(in), optional:: timestep
   logical                      :: radiation_do   ! return value
   !---------------------------------------------------------------------------

   radiation_do = .false.

end function radiation_do

!========================================================================================
subroutine radiation_tend( state, ptend, pbuf, cam_out, cam_in, net_flx, rd_out)

  ! Arguments
  type(physics_state), intent(in), target :: state
  type(physics_ptend), intent(out)        :: ptend

  type(physics_buffer_desc), pointer      :: pbuf(:)
  type(cam_out_t),     intent(inout)      :: cam_out
  type(cam_in_t),      intent(in)         :: cam_in
  real(r8),            intent(out)        :: net_flx(pcols)

  type(rad_out_t), target, optional, intent(out) :: rd_out


  ! this stub can be called, but does nothing

end subroutine radiation_tend

end module radiation
