module radconstants

! provide stubs to allow building with no radiation scheme active

use shr_kind_mod,   only: r8 => shr_kind_r8
use cam_abortutils, only: endrun

implicit none
private
save

integer, parameter, public :: nswbands = 1
integer, parameter, public :: nlwbands = 1
integer, parameter, public :: idx_sw_diag = 1
integer, parameter, public :: idx_lw_diag = 1
integer, parameter, public :: idx_nir_diag = 1
integer, parameter, public :: idx_uv_diag = 1

public :: rad_gas_index
public :: get_lw_spectral_boundaries, get_sw_spectral_boundaries

integer, public, parameter :: gasnamelength = 5
integer, public, parameter :: nradgas = 2
character(len=gasnamelength), public, parameter :: gaslist(nradgas) &
   = (/'CO2  ','H2O  '/)

!========================================================================================
contains
!========================================================================================

integer function rad_gas_index(gasname)
   ! return the index in the gaslist array of the specified gasname

   character(len=*),intent(in) :: gasname
   integer :: igas

   rad_gas_index = -1
   do igas = 1, nradgas
      if (trim(gaslist(igas)).eq.trim(gasname)) then
         rad_gas_index = igas
         return
      endif
   enddo
   call endrun ("rad_gas_index: can not find gas with name "//gasname)
end function rad_gas_index

!------------------------------------------------------------------------------

subroutine get_lw_spectral_boundaries(low_boundaries, high_boundaries, units)
   ! stub should not be called

   real(r8), intent(out) :: low_boundaries(nlwbands), high_boundaries(nlwbands)
   character(*), intent(in) :: units ! requested units

   call endrun('get_lw_spectral_boundaries: ERROR: this is a stub')

end subroutine get_lw_spectral_boundaries

!------------------------------------------------------------------------------

subroutine get_sw_spectral_boundaries(low_boundaries, high_boundaries, units)
   ! stub should not be called

   real(r8), intent(out) :: low_boundaries(nswbands), high_boundaries(nswbands)
   character(*), intent(in) :: units ! requested units

   call endrun('get_sw_spectral_boundaries: ERROR: this is a stub')

end subroutine get_sw_spectral_boundaries

!------------------------------------------------------------------------------

end module radconstants
