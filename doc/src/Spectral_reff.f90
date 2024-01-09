subroutine Spectral_reff(sp, r_med, sigma)
  use Rosseland_data_mod
  use lxmie_mod, only : lxmie
  implicit none

  real(dp), intent(in) :: r_med, sigma
  character(len=100), intent(in) :: sp

  integer :: l, u_reff
  complex(dp) :: ri
  real(dp) :: x, xsec, q_ext, q_sca, q_abs, g, r_eff

  real(dp), allocatable, dimension(:) :: kext_l, a_l, g_l

  ! Perform Mie calculations
  allocate(kext_l(nwl), a_l(nwl), g_l(nwl))

  ! Find effective radius from lognormal parameters
  r_eff = r_med * exp(5.0_dp/2.0_dp * log(sigma)**2)

  do l = 1, nwl

    print*, l, nwl, real(wl(l))

    ! Complex refractive index
    ri = cmplx(n_int(l),-k_int(l),dp)

    ! Size parameter
    x = (twopi * r_eff)/wl(l)

    ! cross section
    xsec = pi * (r_eff*1e-4_dp)**2

    ! Call Mie routine for this size and wavelength
    call lxmie(ri, x, q_ext, q_sca, q_abs, g)

    kext_l(l) = xsec * q_ext ! Units of kext_l [cm2]
    a_l(l) = q_sca/q_ext
    g_l(l) = g

  end do

  ! Output table for each wavelanth
  open(newunit=u_reff, file='results_Spectral_reff/'//trim(sp)//'_reff.txt',action='readwrite')
  write(u_reff,*) r_med, r_eff, sigma, nwl
  do l = 1, nwl
    write(u_reff,*) wl(l), kext_l(l), a_l(l), g_l(l)
    call flush(u_reff)
  end do

  close(u_reff)

end subroutine Spectral_reff