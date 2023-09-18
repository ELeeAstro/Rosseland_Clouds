subroutine lognorm_clouds(wl_path,sp,nk_path,r_med,sigma,rmin,rmax,iint)  
  use, intrinsic :: iso_fortran_env
  use ieee_arithmetic
  use lxmie_mod, only : lxmie
  implicit none

  integer, parameter :: dp = REAL64

  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi

  integer :: iint
  real(dp), intent(in) :: r_med, sigma, rmin, rmax
  character(len=100), intent(in) ::  wl_path, sp, nk_path

  logical :: c_flag
  integer :: r, l, n
  integer :: nlines, nwl
  integer :: uln, uwl, unk

  real(dp), dimension(:), allocatable :: rad, wl
  real(dp), dimension(:), allocatable :: wl_ori, n_ori, k_ori, n_int, k_int

  real(dp), dimension(:), allocatable :: kext_l, a_l, g_l, nd_dist
  real(dp), dimension(:), allocatable :: ln_knorm, ln_a, ln_g

  complex(dp) :: ri
  real(dp) :: x, q_ext, q_sca, q_abs, g

  real(dp) :: lrmin, lrmax, xsec

  ! Read in wavelength grid
  open(newunit=uwl, file=trim(wl_path),action='read')
  read(uwl,*) nwl
  allocate(wl(nwl))
  do l = 1, nwl
    read(uwl,*) wl(l)
    !print*, l, wl(l)
  end do

  ! Read in n,k constants for species - Note comment out/in lines for 'MgSiO3_2' lines (Xianyu's original data)
  open(newunit=unk, file=trim(nk_path),action='read')
  print*, trim(nk_path)
  read(unk,*) nlines, c_flag
  allocate(wl_ori(nlines),n_ori(nlines),k_ori(nlines))
  read(unk,*) ; read(unk,*); read(unk,*); read(unk,*)
  do n = 1, nlines
  !do n = nlines, 1, -1 !'MgSiO3_2' - comment in/out
   read(unk,*) wl_ori(n),n_ori(n),k_ori(n)
    !wl_ori(n) = 1.0_dp/wl_ori(n) * 1e4 ! 'MgSiO3_2' - comment in/out
    n_ori(n) = max(0.0_dp,n_ori(n))
    k_ori(n) = max(0.0_dp,k_ori(n))
    !print*,  n, wl_ori(n),n_ori(n),k_ori(n)
  end do

  ! Interpolate n,k constants to wavelength grid
  allocate(n_int(nwl),k_int(nwl))
  call interp_nk(c_flag, nwl, nlines, wl, wl_ori, n_ori, k_ori, n_int, k_int)

  ! Perform Mie calculations
  allocate(ln_knorm(nwl),ln_a(nwl),ln_g(nwl))
  allocate(kext_l(iint), a_l(iint), g_l(iint))
  allocate(nd_dist(iint))

  ! Find radius grid
  allocate(rad(iint))
  lrmin = log10(rmin)
  lrmax = log10(rmax)
  do r = 1, iint
    rad(r) = 10.0_dp**((lrmax-lrmin) * real(r-1,dp) / real(iint-1,dp) + lrmin)
  end do

  do l = 1, nwl

    print*, l, nwl, real(wl(l))

    ! Complex refractive index
    ri = cmplx(n_int(l),-k_int(l),dp)

    do r = 1, iint

      ! Normalised distribution [um-1]
      nd_dist(r) = (1.0_dp  / (rad(r) * sqrt(twopi) * log(sigma))) * &
          & exp(-(log(rad(r)/r_med))**2/(2.0_dp * log(sigma)**2))

      if ((ieee_is_nan(nd_dist(r)) .eqv. .True.) .or. (ieee_is_finite(nd_dist(r)) .eqv. .False.)) then
        nd_dist(r) = 1.0e-99_dp
      end if

      ! Limiter for very low numbers
      nd_dist(r) = max(nd_dist(r),1.0e-99_dp)

      ! Size parameter
      x = (twopi * rad(r))/wl(l)

      ! Particle cross section [um2]
      xsec = pi * rad(r)**2

      ! Call Mie routine for this size and wavelength
      call lxmie(ri, x, q_ext, q_sca, q_abs, g)

      !print*, rad(r), wl(l), q_ext, q_sca, g
      kext_l(r) = q_ext * xsec * nd_dist(r) ! Units of kext_l [um]
      a_l(r) = q_sca/q_ext
      g_l(r) = g

    end do

    ! Integrate to find integrated properties for this wavelength
    ! Total extinction = integral over all sizes - units of ln_knorm [um2]
    ln_knorm(l) = trapz(rad(:),kext_l(:))
    ! effective SSA = integral for k_sca / k_ext =  integral(k_ext * a) / k_ext
    ln_a(l) = trapz(rad(:),kext_l(:) * a_l(:)) ! Store intermediate result for cl_out_g
    ! effective g = itergral for k_sca * g / k_sca = integeral(k_sca * g) / k_sca
    ln_g(l) = trapz(rad(:),kext_l(:) * a_l(:) * g_l(:))/ln_a(l)

    ln_a(l) = ln_a(l) / ln_knorm(l) ! Albedo is scattering/extinction

    ! Convert units of ln_knorm to [cm2]
    ! So when multipled by N0 [cm-3] and divided by density [g cm-3] cloud opacity becomes cm2 g-1
    ln_knorm(l) = ln_knorm(l) * 1.0e-8_dp

  end do

  ! Output table for each wavelanth
  open(newunit=uln, file='results_lognorm/'//trim(sp)//'_ln.txt',action='readwrite')
  write(uln,*) r_med, sigma, rmin, rmax, iint, nwl
  do l = 1, nwl
    write(uln,*) wl(l), ln_knorm(l), ln_a(l), ln_g(l)
    call flush(uln)
  end do

contains

  pure function trapz(xx, yy) result(rr)
    !! Calculates the integral of an array y with respect to x using the trapezoid
    !! approximation. Note that the mesh spacing of x does not have to be uniform.
    real(dp), intent(in)  :: xx(:)         !! Variable x
    real(dp), intent(in)  :: yy(size(xx))   !! Function y(x)
    real(dp)              :: rr            !! Integral ∫y(x)·dx

    ! Integrate using the trapezoidal rule
    associate(n => size(xx))
      rr = sum((yy(1+1:n-0) + yy(1+0:n-1))*(xx(1+1:n-0) - xx(1+0:n-1)))/2.0_dp
    end associate
  end function trapz

end subroutine lognorm_clouds