program main_select
  use, intrinsic :: iso_fortran_env
  implicit none

  integer, parameter :: dp = REAL64

  integer :: meth
  character(len=100) :: wl_path, sp, nk_path
  character(len=100) :: RTtab_path

  integer :: iint
  real(dp) :: r_med, sigma, rmin, rmax

  meth = 2


  !! Wavelength file path
  wl_path = 'wavelengths_RM.txt'

  !! Species name
  sp = 'Mg2SiO4_amorph'

  !! Path to nk data
  nk_path = 'nk/'//trim(sp)//'[s].dat'

  !! Rosseland mean radius-temperature table
  RTtab_path = 'RTtable.txt'

  !! lognormal distribution properties
  r_med = 1.0_dp
  sigma = 1.5_dp
  rmin = 1e-3_dp
  rmax = 100.0_dp

  iint = 100

  select case(meth)
  case(1)
    call Rosseland_clouds(RTtab_path,wl_path,sp,nk_path)
  case(2)
    call lognorm_clouds(wl_path,sp,nk_path,r_med,sigma,rmin,rmax,iint)
  case default
    print*, 'Invalid method integer selected', meth
  end select

end program main_select