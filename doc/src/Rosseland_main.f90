program Rosseland_main
  use Rosseland_data_mod
  implicit none

  integer :: meth
  character(len=100) :: wl_path, sp, nk_path, RT_path

  integer :: iint
  real(dp) :: N0, r_med, sigma, rmin, rmax

  integer :: u_nml

  namelist /Rosseland_clouds_nml/ &
    & meth, wl_path, sp, nk_path, RT_path, &
    & N0, r_med, sigma, rmin, rmax, iint


  !! Read input variables from main namelist
  open(newunit=u_nml, file='Rosseland_clouds.nml', status='old', action='read')
  read(u_nml, nml=Rosseland_clouds_nml)
  close(u_nml)

  !print*, meth, wl_path, sp, nk_path, RT_path
  !print*, r_med, sigma, rmin, rmax, iint


  call read_wavelength_table(wl_path)
  call read_RT_table(RT_path)
  call read_nk_data(nk_path)
  call interp_nk()

  select case(meth)
  case(1)
    call Rosseland_single(sp)
  case(2)
    call Rosseland_lognorm(sp, N0, sigma, rmin, rmax, iint)
  case(3)
    call Rosseland_reff(sp, sigma) 
  case(4)
    call Spectral_single(sp, r_med)
  case(5)
    call Spectral_lognorm(sp, N0, r_med, sigma, rmin, rmax, iint)       
  case(6)
    call Spectral_reff(sp, r_med, sigma)
  case(7)
    call Planck_single(sp)
  case(8)
    call Planck_lognorm(sp, N0, sigma, rmin, rmax, iint)
  case(9)
    call Planck_reff(sp, sigma) 
  case default
    print*, 'Invalid method integer selected', meth
  end select

end program Rosseland_main
