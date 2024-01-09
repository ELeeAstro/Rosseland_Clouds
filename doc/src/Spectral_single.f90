subroutine Spectral_single(sp, r_med)
  use Rosseland_data_mod
  use lxmie_mod, only : lxmie
  implicit none

  real(dp), intent(in) :: r_med
  character(len=100), intent(in) :: sp

  integer :: l, u_rmed
  complex(dp) :: ri
  real(dp) :: x, xsec, q_ext, q_sca, q_abs, g

  real(dp), allocatable, dimension(:) :: kext_l, a_l, g_l

  ! Perform Mie calculations
  allocate(kext_l(nwl), a_l(nwl), g_l(nwl))

  do l = 1, nwl

    print*, l, nwl, real(wl(l))

    ! Complex refractive index
    ri = cmplx(n_int(l),-k_int(l),dp)

    ! Size parameter
    x = (twopi * r_med)/wl(l)

    ! cross section
    xsec = pi * (r_med*1e-4_dp)**2

    ! Call Mie routine for this size and wavelength
    call lxmie(ri, x, q_ext, q_sca, q_abs, g)

    kext_l(l) = xsec * q_ext ! Units of kext_l [cm2]
    a_l(l) = q_sca/q_ext
    g_l(l) = g

  end do

  ! Output table for each wavelanth
  open(newunit=u_rmed, file='results_Spectral_single/'//trim(sp)//'_single.txt',action='readwrite')
  write(u_rmed,*) r_med, nwl
  do l = 1, nwl
    write(u_rmed,*) wl(l), kext_l(l), a_l(l), g_l(l)
    call flush(u_rmed)
  end do

  close(u_rmed)

end subroutine Spectral_single