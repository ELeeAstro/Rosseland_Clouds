subroutine Planck_lognorm(sp, N0, sigma, rmin, rmax, iint)
  use Rosseland_data_mod
  use lxmie_mod, only : lxmie
  implicit none

  character(len=100), intent(in) :: sp
  integer, intent(in) :: iint
  real(dp), intent(in) :: N0, sigma, rmin, rmax

  integer :: i, aa, tt, l, u_k, u_a, u_g

  real(dp), dimension(:), allocatable :: kext_r, a_r, g_r, nd_r, r_r
  real(dp), dimension(:), allocatable :: kext_l, a_l, g_l
  real(dp), dimension(:, :), allocatable :: Pl_kext, Pl_a, Pl_g

  complex(dp) :: ri
  real(dp) :: x, xsec, q_ext, q_sca, q_abs, g
  real(dp) :: lamin, lamax


  !! Perform Mie calculations

  !! allocate work arrays for each log-normal size
  allocate(kext_r(iint), a_r(iint), g_r(iint), r_r(iint), nd_r(iint)) 

  !! allocate work arrays for each wavelength
  allocate(kext_l(nwl), a_l(nwl), g_l(nwl))

  !! allocate end Planck mean arrays
  allocate(Pl_kext(na,nT),Pl_a(na,nT),Pl_g(na,nT))

  !! Calculate sizes to integrate lognormal function from
  lamin = log10(rmin)
  lamax = log10(rmax)
  do i = 1, iint
    r_r(i) = 10.0_dp**((lamax-lamin) * real(i-1,kind=dp) / real(iint-1,kind=dp) + lamin)
    !print*, i, r_r(i)
  end do

  do aa = 1, na

    print*, a(aa), aa, na

    do tt = 1, nt

      do l = 1, nwl

        ! Optical constant
        ri = cmplx(n_int(l),-k_int(l),dp)

        do i = 1, iint

          !! Number density assuming size r_r(i) at median size a(aa)
          nd_r(i) = N0/(r_r(i)*1e-4_dp*sqrt(2.0_dp*pi)*log(sigma)) *  &
            & exp(-(log(r_r(i)/(a(aa)))**2)/(2.0_dp*log(sigma)**2))

          ! Limiter for underflow
          nd_r(i) = max(nd_r(i), 1e-30_dp)

          ! Cross section
          xsec = pi * (r_r(i)* 1e-4_dp)**2

          ! Size parameter
          x = (twopi * r_r(i))/wl(l)

          call lxmie(ri, x, q_ext, q_sca, q_abs, g)

          kext_r(i) = xsec * q_ext * nd_r(i)
          a_r(i) = xsec * q_sca * nd_r(i)
          g_r(i) = max(g,1.0e-12_dp) * xsec * q_sca * nd_r(i)
        
        end do 

        !! Integrated normalised log-normal properties at wavelength l
        kext_l(l) = trapz(r_r(:)*1e-4_dp,kext_r(:))
        a_l(l) = trapz(r_r(:)*1e-4_dp,a_r(:))
        g_l(l) = trapz(r_r(:)*1e-4_dp,g_r(:))
        !print*, a(aa), wl(l), x, q_ext, q_sca, g

      end do

      call Planck_mean(nwl, wl(:),T(tt),kext_l(:),Pl_kext(aa,tt))
      call Planck_mean(nwl, wl(:),T(tt),a_l(:),Pl_a(aa,tt))
      call Planck_mean(nwl, wl(:),T(tt),g_l(:),Pl_g(aa,tt))


      !! Now calculate ssa and g

      ! g is scattering opacity weighted by g divided by scattering opacity
      Pl_g(aa,tt) = Pl_g(aa,tt)/Pl_a(aa,tt)

      ! ssa is scattering opacity divided by extinction opacity
      Pl_a(aa,tt) = Pl_a(aa,tt)/Pl_kext(aa,tt)

    end do
  end do


  ! Output table - go radius outer loop, temperature inner loop
  ! Output cross section, ssa and g

  open(newunit=u_k, file='results_Planck_lognorm/'//trim(sp)//'_kext.txt',action='readwrite')
  write(u_k,*) na, nT, sigma
  write(u_k,*) a(:)
  write(u_k,*) T(:)
  do aa = 1, na
    write(u_k,*) (real(Pl_kext(aa,tt)), tt = 1, nT)
    call flush(u_k)
  end do

  close(u_k)

  open(newunit=u_a, file='results_Planck_lognorm/'//trim(sp)//'_a.txt',action='readwrite')
  write(u_a,*) na, nT, sigma
  write(u_a,*) a(:)
  write(u_a,*) T(:)
  do aa = 1, na
    write(u_a,*) (real(Pl_a(aa,tt)), tt = 1, nT)
    call flush(u_a)
  end do

  close(u_a)

  open(newunit=u_g, file='results_Planck_lognorm/'//trim(sp)//'_g.txt',action='readwrite')
  write(u_g,*) na, nT, sigma
  write(u_g,*) a(:)
  write(u_g,*) T(:)
  do aa = 1, na
    write(u_g,*) (real(Pl_g(aa,tt)), tt = 1, nT)
    call flush(u_g)
  end do

  close(u_g)


end subroutine Planck_lognorm
