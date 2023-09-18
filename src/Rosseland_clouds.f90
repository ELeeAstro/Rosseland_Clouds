subroutine Rosseland_clouds(RTtab_path,wl_path,sp,nk_path)
  use, intrinsic :: iso_fortran_env
  use lxmie_mod, only : lxmie
  implicit none

  integer, parameter :: dp = REAL64

  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi

  character(len=100), intent(in) :: RTtab_path, wl_path, sp, nk_path

  logical :: c_flag
  integer :: r, t, l, n
  integer :: nr, nt, nlines, nwl
  integer :: uQext, uQsca, ugg, uin, uwl, unk

  real(dp), dimension(:), allocatable :: rad, temp, wl
  real(dp), dimension(:), allocatable :: wl_ori, n_ori, k_ori, n_int, k_int

  real(dp), dimension(:), allocatable :: Qext_l, Qsca_l, gg_l
  real(dp), dimension(:, :), allocatable :: Ross_Qext, Ross_Qsca, Ross_gg

  complex(dp) :: ri
  real(dp) :: x, q_ext, q_sca, q_abs, g

  ! Read in temperature and grain size grid - NOTE: different ones depending on T and a grid
  open(newunit=uin, file=trim(RTtab_path),action='read')
  read(uin,*) nr, nt
  allocate(temp(nt),rad(nr))
  read(uin,*) (rad(r),r=1,nr)
  read(uin,*) (temp(t),t=1,nt)

  print*, rad(:)
  print*, temp(:)

  print*, nr, nt

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
  allocate(Ross_Qext(nr,nt),Ross_Qsca(nr,nt),Ross_gg(nr,nt))
  allocate(Qext_l(nwl), Qsca_l(nwl), gg_l(nwl))

  do r = 1, nr
    print*, r, nr, real(rad(r))
    do t = 1, nt
      do l = 1, nwl

        ri = cmplx(n_int(l),-k_int(l),dp)
        x = (twopi * rad(r))/wl(l)

        call lxmie(ri, x, q_ext, q_sca, q_abs, g)

        Qext_l(l) = q_ext
        Qsca_l(l) = q_sca
        gg_l(l) = max(g,1.0e-12_dp)
        gg_l(l) = min(gg_l(l), 1.0_dp)

        !print*, rad(r), wl(l), q_ext, q_sca, g

      end do

      call Ross_mean(nwl, wl(:),temp(t),Qext_l(:),Ross_Qext(r,t))
      call Ross_mean(nwl, wl(:),temp(t),Qsca_l(:),Ross_Qsca(r,t))
      call Ross_mean(nwl, wl(:),temp(t),gg_l(:),Ross_gg(r,t))

    end do
  end do


  ! Output table - go radius outer loop, temperature inner loop
  ! Output Qext - Qsca and g

  open(newunit=uQext, file='results_Ross/'//trim(sp)//'_rosselandMean_qext.txt',action='readwrite')
  write(uQext,*) nr, nt
  do r = 1, nr
    write(uQext,*) (real(Ross_Qext(r,t)), t = 1, nt)
    call flush(uQext)
  end do

  open(newunit=uQsca, file='results_Ross/'//trim(sp)//'_rosselandMean_qscat.txt',action='readwrite')
  write(uQsca,*) nr, nt
  do r = 1, nr
    write(uQsca,*) (real(Ross_Qsca(r,t)), t = 1, nt)
    call flush(uQsca)
  end do

  open(newunit=ugg, file='results_Ross/'//trim(sp)//'_rosselandMean_gg.txt',action='readwrite')
  write(ugg,*) nr, nt
  do r = 1, nr
    write(ugg,*) (real(Ross_gg(r,t)), t = 1, nt)
    call flush(ugg)
  end do


end subroutine Rosseland_clouds

subroutine Ross_mean(nwl, wl, temp, Vl, Vr)
  use, intrinsic :: iso_fortran_env
  implicit none

  integer, parameter :: dp = REAL64

  real(dp), parameter :: hp = 6.62607015e-27_dp ! erg s - Planck's constant
  real(dp), parameter :: c_s = 2.99792458e10_dp ! cm s^-1 - Vacuum speed of light
  real(dp), parameter :: kb = 1.380649e-16_dp ! erg K^-1 - Boltzmann's constant

  integer, intent(in) :: nwl
  real(dp), dimension(nwl), intent(in) :: wl, Vl
  real(dp), intent(in) :: temp

  real(dp), intent(out) :: Vr

  integer :: l
  real(dp) :: top, bot,  xx
  real(dp), dimension(nwl) :: dBdT, wl_cm

  !! Subroutine calculates Rosseland mean weighted value
  !! Returns Rosseland mean to all array values

  do l = 1, nwl

    wl_cm(l) = wl(l) * 1e-4_dp
    xx = (hp * c_s) / (wl_cm(l) * kb * temp)
    top = 2.0_dp * hp**2 * c_s**3 * exp(xx)
    bot = wl_cm(l)**6 * kb * temp**2 * (exp(xx) - 1.0_dp)**2
    dBdT(l) = top/bot
  end do

  top = trapz(wl_cm(:),(1.0_dp/Vl(:))*dBdT(:))
  bot = trapz(wl_cm(:),dBdT(:))

  Vr = 1.0_dp / (top/bot)

contains

  pure function trapz(x, y) result(r)
    !! Calculates the integral of an array y with respect to x using the trapezoid
    !! approximation. Note that the mesh spacing of x does not have to be uniform.
    real(dp), intent(in)  :: x(:)         !! Variable x
    real(dp), intent(in)  :: y(size(x))   !! Function y(x)
    real(dp)              :: r            !! Integral ∫y(x)·dx

    ! Integrate using the trapezoidal rule
    associate(n => size(x))
      r = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2.0_dp
    end associate
  end function trapz

end subroutine Ross_mean
