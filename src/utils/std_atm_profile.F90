module std_atm_profile
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------

use shr_kind_mod,        only: r8 => shr_kind_r8
use cam_logfile,         only: iulog
use cam_abortutils,      only: endrun
use physconst,           only: pi

implicit none
private
save

public :: &
   std_atm_pres,   & ! compute pressure given height
   std_atm_height, & ! compute height given pressure
   std_atm_temp      ! compute temperature given height

! Parameters
  integer,  parameter :: namp = 6 
  ! coefficients for 5-degree polynomial fit to vertical structure
  real(r8), parameter :: amp0(namp) = &
     (/ 1.90e2_r8, -2.14e0_r8 ,   4.68e-2_r8,-4.52e-4_r8  ,  1.52e-6_r8, -5.86e-10_r8/) 
  real(r8), parameter :: amp1(namp) = &
       (/ -3.61e1_r8 , 5.27e-1_r8,  1.53e-2_r8,-3.45e-4_r8  ,  2.29e-6_r8, -5.02e-9_r8/)
  real(r8), parameter :: amp2(namp) = &
       (/ -8.36e0_r8,  2.22e-1_r8, -7.06e-3_r8, 9.44e-5_r8  , -4.82e-7_r8,  8.46e-10_r8/)
  real(r8), parameter :: amp3(namp) = &
       (/ -3.70e-1_r8, 3.33e-1_r8 , -1.74e-2_r8,2.53e-4_r8  , -1.42e-6_r8,  2.80e-9_r8/)
  real(r8), parameter :: amp4(namp) = &
       (/ 1.72e0_r8, -3.29e-2_r8,   -4.64e-4_r8, -2.80e-6_r8,  1.24e-7_r8, -4.84e-10_r8/)
  !
  ! http://www-mars.lmd.jussieu.fr/mars/info_web/MCD4.3_ddd.pdf
  !
  real(r8), parameter :: scale_height_ref =  10.E3_r8 ! [m]
  real(r8), parameter :: psurf_ref        = 610.0_r8  ! reference surface pressure     
!=========================================================================================
CONTAINS
!=========================================================================================

subroutine std_atm_pres(height, pstd)
    
   ! arguments
   real(r8), intent(in)  :: height(:) ! height above sea level in meters
   real(r8), intent(out) :: pstd(:)   ! std pressure in Pa
    
   integer :: k, nlev
   character(len=*), parameter :: routine = 'std_atm_pres'
   !----------------------------------------------------------------------------
    
   nlev = size(height)
   do k = 1, nlev
     pstd(k) = psurf_ref*exp(-height(k)/scale_height_ref)
   end do

end subroutine std_atm_pres

!=========================================================================================

subroutine std_atm_height(pstd, height)
    
   ! arguments
   real(r8), intent(in)   :: pstd(:)   ! std pressure in Pa
   real(r8), intent(out)  :: height(:) ! height above sea level in meters
    
   integer :: k, nlev
   character(len=*), parameter :: routine = 'std_atm_height'
   !----------------------------------------------------------------------------
    
   nlev = size(height)
   do k = 1, nlev
     height(k) = -scale_height_ref*log((pstd(k)/psurf_ref))     
   end do
end subroutine std_atm_height

!=========================================================================================

subroutine std_atm_temp(height, lat, temp)
    
   ! arguments
   real(r8), intent(in)   :: height(:) ! std pressure in Pa
   real(r8), intent(in)   :: lat       ! std pressure in Pa
   real(r8), intent(out)  :: temp(:)   ! temperature
    
   ! local vars
   integer :: k, nlev
   real(r8):: amp(5), z, z2, z3, z4, z5
   real(r8):: y, cosy, cos2y, cos3y, cos4y
   character(len=*), parameter :: routine = 'std_atm_temp'
   !----------------------------------------------------------------------------
    
   nlev = size(height)
   
   y     = 2._r8*(lat + 0.5_r8*pi) ! Fourier expansion done on domain 0 to 2pi
   cosy  = cos(      y)
   cos2y = cos(2._r8*y)
   cos3y = cos(3._r8*y)
   cos4y = cos(4._r8*y)
      
   do k = 1, nlev
     z = height(k)*1.0E-3_r8!convert to km's
     z = min(z,178.0_r8)    !isothermal above 178km
     z2 = z*z
     z3 = z*z2
     z4 = z2*z2
     z5 = z2*z3
     
     amp(1) = amp0(1) + amp0(2)*z + amp0(3)*z2 + amp0(4)*z3 + amp0(5)*z4 + amp0(6)*z5
     amp(2) = amp1(1) + amp1(2)*z + amp1(3)*z2 + amp1(4)*z3 + amp1(5)*z4 + amp1(6)*z5
     amp(3) = amp2(1) + amp2(2)*z + amp2(3)*z2 + amp2(4)*z3 + amp2(5)*z4 + amp2(6)*z5
     amp(4) = amp3(1) + amp3(2)*z + amp3(3)*z2 + amp3(4)*z3 + amp3(5)*z4 + amp3(6)*z5
     amp(5) = amp4(1) + amp4(2)*z + amp4(3)*z2 + amp4(4)*z3 + amp4(5)*z4 + amp4(6)*z5

     temp(k) = amp(1) + amp(2)*cosy  + amp(3)*cos2y &
                      + amp(4)*cos3y + amp(5)*cos4y
   end do
end subroutine std_atm_temp

end module std_atm_profile
