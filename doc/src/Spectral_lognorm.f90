subroutine Spectral_lognorm(sp, N0, r_med, sigma, rmin, rmax, iint)
  use Rosseland_data_mod
  use lxmie_mod, only : lxmie
  implicit none

  character(len=100), intent(in) :: sp
  integer, intent(in) :: iint
  real(dp), intent(in) :: N0, r_med, sigma, rmin, rmax

  integer :: l, i, u_ln
  complex(dp) :: ri
  real(dp) :: lamin, lamax
  real(dp) :: x, xsec, q_ext, q_sca, q_abs, g

  real(dp), allocatable, dimension(:) :: kext_r, a_r, g_r, nd_r, r_r
  real(dp), allocatable, dimension(:) :: kext_l, a_l, g_l

  ! Perform Mie calculations
  allocate(kext_l(nwl), a_l(nwl), g_l(nwl))

  !! allocate work arrays for each log-normal size
  allocate(kext_r(iint), a_r(iint), g_r(iint), r_r(iint), nd_r(iint))

  !! Calculate sizes to integrate lognormal function from
  lamin = log10(rmin)
  lamax = log10(rmax)
  do i = 1, iint
    r_r(i) = 10.0_dp**((lamax-lamin) * real(i-1,kind=dp) / real(iint-1,kind=dp) + lamin)
    !print*, i, r_r(i)
  end do

  do i = 1, iint
    !! Number density assuming size r_r(i) at median size r_med
    nd_r(i) = N0/(r_r(i)*1e-4_dp*sqrt(2.0_dp*pi)*log(sigma)) *  &
    & exp(-(log(r_r(i)/r_med)**2)/(2.0_dp*log(sigma)**2))
    nd_r(i) = max(nd_r(i), 1e-30_dp)
  end do

  do l = 1, nwl

    print*, l, nwl, real(wl(l))

    ! Complex refractive index
    ri = cmplx(n_int(l),-k_int(l),dp)

    do i = 1, iint

      ! Size parameter
      x = (twopi * r_r(i))/wl(l)

      ! cross section
      xsec = pi * (r_r(i)*1e-4_dp)**2

      ! Call Mie routine for this size and wavelength
      call lxmie(ri, x, q_ext, q_sca, q_abs, g)

      kext_r(i) = xsec * q_ext * nd_r(i)
      a_r(i) = xsec * q_sca * nd_r(i)
      g_r(i) = max(g,1.0e-12_dp) * xsec * q_sca * nd_r(i)

    end do

    !! Integrated normalised log-normal properties at wavelength l
    kext_l(l) = trapz(r_r(:)*1e-4_dp,kext_r(:))
    a_l(l) = trapz(r_r(:)*1e-4_dp,a_r(:))
    g_l(l) = trapz(r_r(:)*1e-4_dp,g_r(:))

    !! Now calculate ssa and g

    ! g is scattering opacity weighted by g divided by scattering opacity
    g_l(l) = g_l(l)/a_l(l)

    ! ssa is scattering opacity divided by extinction opacity
    a_l(l) = a_l(l)/kext_l(l)

  end do

  ! Output table for each wavelanth
  open(newunit=u_ln, file='results_Spectral_lognorm/'//trim(sp)//'_lognorm.txt',action='readwrite')
  write(u_ln,*) r_med, sigma, rmin, rmax, iint, nwl
  do l = 1, nwl
    write(u_ln,*) wl(l), kext_l(l), a_l(l), g_l(l)
    call flush(u_ln)
  end do

  close(u_ln)

end subroutine Spectral_lognorm