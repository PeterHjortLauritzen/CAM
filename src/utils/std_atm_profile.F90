module std_atm_profile

#ifdef planet_mars
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

#else

!-------------------------------------------------------------------------------
!
! The barometric formula for U.S. Standard Atmosphere is valid up to 86 km.
! see https://en.wikipedia.org/wiki/Barometric_formula.
!
! N.B.  The extension above 86 km is using data from Hanli.  It is not complete
!       since the hardcoded parameter (c1) needs adjustment above 86 km.
!
!-------------------------------------------------------------------------------

use shr_kind_mod,        only: r8 => shr_kind_r8
use cam_logfile,         only: iulog
use cam_abortutils,      only: endrun

implicit none
private
save

public :: &
   std_atm_pres,   & ! compute pressure given height
   std_atm_height, & ! compute height given pressure
   std_atm_temp      ! compute temperature given height

! Parameters for barometric formula for U.S. Standard Atmosphere.

integer, parameter  :: nreg = 15  ! number of regions

real(r8), parameter :: hb(nreg) = & ! height at bottom of layer (m)
     (/0.0_r8, 1.1e4_r8, 2.0e4_r8, 3.2e4_r8, 4.7e4_r8, 5.1e4_r8, 7.1e4_r8, 8.6e4_r8, &
     9.1e4_r8, 1.1e5_r8, 1.2e5_r8, 1.5e5_r8, 2.0e5_r8, 3.0e5_r8, 7.e5_r8/)

real(r8), parameter :: pb(nreg) = & ! standard pressure (Pa)
     (/101325._r8, 22632.1_r8, 5474.89_r8, 868.02_r8, 110.91_r8, 66.94_r8, 3.96_r8, 3.7e-1_r8,  &
     1.5e-1_r8, 7.1e-3_r8, 2.5e-3_r8, 4.5e-4_r8, 8.47e-5_r8, 8.77e-6_r8, 3.19e-8_r8/)

real(r8), parameter :: tb(nreg) = & ! standard temperature (K)
     (/288.15_r8, 216.65_r8, 216.65_r8, 228.65_r8, 270.65_r8, 270.65_r8, 214.65_r8, 186.87_r8,  &
     186.87_r8, 240._r8, 360._r8, 634.39_r8, 854.56_r8, 976.01_r8, 1.e3_r8/)

real(r8), parameter :: lb(nreg) = & ! temperature lapse rate (K/m)
     (/-0.0065_r8, 0.0_r8, 0.001_r8, 0.0028_r8, 0.0_r8, -0.0028_r8, -0.001852_r8, 0.0_r8,       &
     2.796e-3_r8, 0.012_r8, 9.15e-3_r8, 4.4e-3_r8, 1.21e-3_r8, 6.e-5_r8, 0.0_r8/)

real(r8), parameter :: rg = 8.3144598_r8 ! universal gas constant (J/mol/K)
real(r8), parameter :: g0 = 9.80665_r8   ! gravitational acceleration (m/s^2)
real(r8), parameter :: mw = 0.0289644_r8 ! molar mass of dry air (kg/mol)
real(r8), parameter :: c1 = g0*mw/rg
  
!=========================================================================================
CONTAINS
!=========================================================================================

subroutine std_atm_pres(height, pstd)
    
   ! arguments
   real(r8), intent(in)  :: height(:) ! height above sea level in meters
   real(r8), intent(out) :: pstd(:)   ! std pressure in Pa
    
   integer :: i, ii, k, nlev
   character(len=*), parameter :: routine = 'std_atm_pres'
   !----------------------------------------------------------------------------
    
   nlev = size(height)
   do k = 1, nlev
      if (height(k) < 0.0_r8) then
         ! Extrapolate below mean sea level using troposphere lapse rate.
         ii = 1
      else
         ! find region containing height
         find_region: do i = nreg, 1, -1
            if (height(k) >= hb(i)) then
               ii = i
               exit find_region
            end if
         end do find_region
      end if
      
      if (lb(ii) /= 0._r8) then
         pstd(k) = pb(ii) * ( tb(ii) / (tb(ii) + lb(ii)*(height(k) - hb(ii)) ) )**(c1/lb(ii))
      else
         pstd(k) = pb(ii) * exp( -c1*(height(k) - hb(ii))/tb(ii) )
      end if
      
   end do

end subroutine std_atm_pres

!=========================================================================================

subroutine std_atm_height(pstd, height)
    
   ! arguments
   real(r8), intent(in)   :: pstd(:)   ! std pressure in Pa
   real(r8), intent(out)  :: height(:) ! height above sea level in meters
    
   integer :: i, ii, k, nlev
   logical :: found_region
   character(len=*), parameter :: routine = 'std_atm_height'
   !----------------------------------------------------------------------------
    
   nlev = size(height)
   do k = 1, nlev
      
      if (pstd(k) <= pb(nreg)) then
         ii = nreg
      else if (pstd(k) > pb(1)) then
         ii = 1
      else
         ! find region containing pressure
         find_region: do i = 2, nreg
            if (pstd(k) > pb(i)) then
               ii = i - 1
               exit find_region
            end if
         end do find_region
      end if

      if (lb(ii) /= 0._r8) then
         height(k) = hb(ii) + (tb(ii)/lb(ii)) * ( (pb(ii)/pstd(k))**(lb(ii)/c1) - 1._r8 )
      else
         height(k) = hb(ii) + (tb(ii)/c1)*log(pb(ii)/pstd(k))
      end if
   end do

end subroutine std_atm_height

!=========================================================================================

subroutine std_atm_temp(height, temp)
    
   ! arguments
   real(r8), intent(in)   :: height(:) ! std pressure in Pa
   real(r8), intent(out)  :: temp(:)   ! temperature
    
   ! local vars
   integer :: i, ii, k, nlev
   character(len=*), parameter :: routine = 'std_atm_temp'
   !----------------------------------------------------------------------------
    
   nlev = size(height)
   do k = 1, nlev
      if (height(k) < 0.0_r8) then
         ii = 1
      else
         ! find region containing height
         find_region: do i = nreg, 1, -1
            if (height(k) >= hb(i)) then
               ii = i
               exit find_region
            end if
         end do find_region
      end if

      if (lb(ii) /= 0._r8) then
         temp(k) = tb(ii) + lb(ii)*(height(k) - hb(ii))
      else
         temp(k) = tb(ii)
      end if
      
   end do

end subroutine std_atm_temp

#endif

end module std_atm_profile
