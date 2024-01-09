subroutine Planck_reff(sp, sigma)
  use Rosseland_data_mod
  use lxmie_mod, only : lxmie
  implicit none

  real(dp), intent(in) :: sigma
  character(len=100), intent(in) :: sp

  integer :: aa, tt, l, u_k, u_a, u_g

  real(dp), dimension(:), allocatable :: kext_l, a_l, g_l
  real(dp), dimension(:,:), allocatable :: Pl_kext, Pl_a, Pl_g

  complex(dp) :: ri
  real(dp) :: x, xsec, q_ext, q_sca, q_abs, g
  real(dp), allocatable, dimension(:) :: reff


  !! Perform Mie calculations


  !! allocate work arrays for each wavelength
  allocate(kext_l(nwl), a_l(nwl), g_l(nwl))

  !! allocate end Rosseland mean arrays
  allocate(Pl_kext(na,nT),Pl_a(na,nT),Pl_g(na,nT))

  !! Allocate reff array
  allocate(reff(na))


  do aa = 1, na

    ! Effective radius for log-normal
    reff(aa) = a(aa) * exp(5.0_dp/2.0_dp * log(sigma)**2)

    print*, a(aa), reff(aa), aa, na

    ! Cross section
    xsec = pi * (reff(aa)* 1e-4_dp)**2

    do tt = 1, nt

      do l = 1, nwl

        ! Optical constant
        ri = cmplx(n_int(l),-k_int(l),dp)

        ! Size parameter
        x = (twopi * reff(aa))/wl(l)

        call lxmie(ri, x, q_ext, q_sca, q_abs, g)

        kext_l(l) = xsec * q_ext
        a_l(l) = xsec * q_sca
        g_l(l) = max(g,1.0e-12_dp) * xsec * q_sca

        !print*, a(aa), wl(l), x, q_ext, q_sca, g

      end do

      call Planck_mean(nwl, wl(:), T(tt), kext_l(:), Pl_kext(aa,tt))
      call Planck_mean(nwl, wl(:), T(tt), a_l(:), Pl_a(aa,tt))
      call Planck_mean(nwl, wl(:), T(tt), g_l(:), Pl_g(aa,tt))

      !! Now calculate ssa and g

      ! g is scattering opacity weighted by g divided by scattering opacity
      Pl_g(aa,tt) = Pl_g(aa,tt)/Pl_a(aa,tt)

      ! ssa is scattering opacity divided by extinction opacity
      Pl_a(aa,tt) = Pl_a(aa,tt)/Pl_kext(aa,tt)

    end do
  end do


  ! Output table - go radius outer loop, temperature inner loop
  ! Output cross section, ssa and g

  open(newunit=u_k, file='results_Planck_reff/'//trim(sp)//'_kext.txt',action='readwrite')
  write(u_k,*) na, nT, sigma
  write(u_k,*) a(:)
  write(u_k,*) reff(:)
  write(u_k,*) T(:)
  do aa = 1, na
    write(u_k,*) (real(Pl_kext(aa,tt)), tt = 1, nT)
    call flush(u_k)
  end do

  close(u_k)

  open(newunit=u_a, file='results_Planck_reff/'//trim(sp)//'_a.txt',action='readwrite')
  write(u_a,*) na, nT, sigma
  write(u_a,*) a(:)
  write(u_a,*) reff(:)
  write(u_a,*) T(:)
  do aa = 1, na
    write(u_a,*) (real(Pl_a(aa,tt)), tt = 1, nT)
    call flush(u_a)
  end do

  close(u_a)

  open(newunit=u_g, file='results_Planck_reff/'//trim(sp)//'_g.txt',action='readwrite')
  write(u_g,*) na, nT, sigma
  write(u_g,*) a(:)
  write(u_g,*) reff(:)
  write(u_g,*) T(:)
  do aa = 1, na
    write(u_g,*) (real(Pl_g(aa,tt)), tt = 1, nT)
    call flush(u_g)
  end do

  close(u_g)


end subroutine Planck_reff
