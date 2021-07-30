program Rosseland_clouds
  use, intrinsic :: iso_fortran_env
  use mie_routines
  implicit none

  integer, parameter :: dp = REAL64

  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi

  character(len=20) :: sp
  logical :: c_flag
  integer :: r, t, l, n
  integer :: nr, nt, nlines, nwl
  integer :: uQext, uQsca, ugg, uin, uwl, unk

  real(dp), dimension(:), allocatable :: rad, temp, wl
  real(dp), dimension(:), allocatable :: wl_ori, n_ori, k_ori, n_int, k_int

  real(dp), dimension(:), allocatable :: Qext_l, Qsca_l, gg_l
  real(dp), dimension(:, :), allocatable :: Ross_Qext, Ross_Qsca, Ross_gg

  integer, parameter :: rnang = 1
  integer :: rier
  complex(dp) :: ri
  real(dp) :: x, rQext, rQsca, rQabs, rQbk, rQpr, ralbedo, rg
  complex(dp), dimension(rnang) :: rSA1, rSA2
  logical, parameter :: rdoSA = .False.

  ! Species name - 'MgSiO3_2' is Xianyu's MgSiO3 data
  !! Change to the species name here and recompile
  !! will auto read in the nk constants in the nk directory
  sp = 'MgSiO3'

  ! Read in temperature and grain size grid
  open(newunit=uin, file='rosselandMean_RTtable.txt',action='read')
  read(uin,*) nr, nt
  allocate(temp(nt),rad(nr))
  read(uin,*) (rad(r),r=1,nr)
  read(uin,*) (temp(t),t=1,nt)

  !print*, rad(:)
  !print*, temp(:)

  ! Read in wavelength grid
  open(newunit=uwl, file='wavelengths.txt',action='read')
  read(uwl,*) nwl
  allocate(wl(nwl))
  do l = 1, nwl
    read(uwl,*) wl(l)
    !print*, l, wl(l)
  end do

  ! Read in n,k constants for species - Note comment out/in lines for 'MgSiO3_2' lines (Xianyu's original data)
  open(newunit=unk, file='nk/'//trim(sp)//'[s].dat',action='read')
  print*, 'nk/'//trim(sp)//'[s].dat'
  read(unk,*) nlines, c_flag
  allocate(wl_ori(nlines),n_ori(nlines),k_ori(nlines))
  read(unk,*) ; read(unk,*); read(unk,*); read(unk,*)
  do n = 1, nlines
  !do n = nlines, 1, -1 !'MgSiO3_2' - comment in/out
   read(unk,*) wl_ori(n),n_ori(n),k_ori(n)
    !wl_ori(n) = 1.0_dp/wl_ori(n) * 1e4 ! 'MgSiO3_2' - comment in/out
    n_ori(n) = max(0.0_dp,n_ori(n))
    k_ori(n) = max(0.0_dp,k_ori(n))
    print*,  n, wl_ori(n),n_ori(n),k_ori(n)
  end do


  ! Interpolate n,k constants to wavelength grid
  allocate(n_int(nwl),k_int(nwl))
  call interp_nk(c_flag, nwl, nlines, wl, wl_ori, n_ori, k_ori, n_int, k_int)

  ! Perform Mie calculations
  allocate(Ross_Qext(nr,nt),Ross_Qsca(nr,nt),Ross_gg(nr,nt))
  allocate(Qext_l(nwl), Qsca_l(nwl), gg_l(nwl))

  do r = 1, nr
    print*, r, nr, real(rad(r)*1.0e6_dp)
    do t = 1, nt
      do l = 1, nwl

        ! Real and complex n,k constants
        ri = cmplx(n_int(l),k_int(l),dp)
        ! Size parameter
        x = (twopi * rad(r) * 1.0e6_dp)/wl(l)

        ! Mie theory solver
        call shexqnn2(ri, x, rQext, rQsca, rQabs, rQbk, rQpr, ralbedo, rg, &
                     & rier, rSA1, rSA2, rdoSA, rnang)

        ! Save results to arrays
        Qext_l(l) = rQext
        Qsca_l(l) = rQsca
        gg_l(l) = max(rg,1.0e-12_dp)
        gg_l(l) = min(gg_l(l), 1.0_dp)

        !print*, rad(r), wl(l), rQext, rQsca, rg

      end do

      ! Rosseland mean weighting calculation (1 = trapezium rule, 2 = simple sum rule)
      call Ross_mean(nwl, wl(:),temp(t),Qext_l(:),Ross_Qext(r,t),1)
      call Ross_mean(nwl, wl(:),temp(t),Qsca_l(:),Ross_Qsca(r,t),1)
      call Ross_mean(nwl, wl(:),temp(t),gg_l(:),Ross_gg(r,t),1)
    end do
  end do


  ! Output table - go radius outer loop, temperature inner loop
  ! Output Qext - Qsca and g

  open(newunit=uQext, file='results/'//trim(sp)//'_rosselandMean_qext.txt',action='readwrite')
  write(uQext,*) nr, nt
  do r = 1, nr
    write(uQext,*) (real(Ross_Qext(r,t)), t = 1, nt)
  end do

  open(newunit=uQsca, file='results/'//trim(sp)//'_rosselandMean_qscat.txt',action='readwrite')
  write(uQsca,*) nr, nt
  do r = 1, nr
    write(uQsca,*) (real(Ross_Qsca(r,t)), t = 1, nt)
  end do

  open(newunit=ugg, file='results/'//trim(sp)//'_rosselandMean_gg.txt',action='readwrite')
  write(ugg,*) nr, nt
  do r = 1, nr
    write(ugg,*) (real(Ross_gg(r,t)), t = 1, nt)
  end do


end program Rosseland_clouds

subroutine interp_nk(conducting, nwl, n_lines, wl, wl_work, n_work, k_work, no, ko)
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  logical, intent(in) :: conducting
  integer, intent(in) :: nwl, n_lines
  real(dp), dimension(n_lines), intent(in) :: wl_work, n_work, k_work
  real(dp), dimension(nwl), intent(in) :: wl
  real(dp), dimension(nwl), intent(out) :: no, ko

  integer :: l, l1
  real(dp) :: fac

  do l = 1, nwl
    ! If required wavelength is less than availible data - keep constant
    if (wl(l) < wl_work(1)) then
      no(l) = n_work(1)
      ko(l) = k_work(1)
      no(l) = max(no(l),0.0_dp)
      ko(l) = max(ko(l),0.0_dp)
    ! If required wavelength is greater than availible data - extrapolate
    ! Non conducting: n is constant - k is linear decreasing
    ! Conducting: n and k are log-log extrapolated
    else if (wl(l) > wl_work(n_lines)) then
      if (conducting .eqv. .False.) then
        no(l) = n_work(n_lines)
        ko(l) = k_work(n_lines)*wl_work(n_lines)/wl(l)
        no(l) = max(no(l),0.0_dp)
        ko(l) = max(ko(l),0.0_dp)
      else if (conducting .eqv. .True.) then
        do l1 = n_lines,1,-1
          if (wl_work(l1) < 0.7_dp*wl_work(n_lines)) then  ! data can be noisy, so
            exit                            ! it's safer to use larger
          end if                            ! region to get the slope
        enddo
        fac = log(wl(l)/wl_work(n_lines))/log(wl_work(l1)/wl_work(n_lines))
        no(l) = exp(log(n_work(n_lines)) &
          & + fac*log(n_work(l1)/n_work(n_lines)))
        ko(l) = exp(log(k_work(n_lines)) &
          & + fac*log(k_work(l1)/k_work(n_lines)))
        no(l) = max(no(l),0.0_dp)
        ko(l) = max(ko(l),0.0_dp)
      end if
    ! Data is availible in the required wavelength range - log-log interpolation
    else
      ! Loop across work arrays untill straddle point is point then interpolate
      do l1 = 1, n_lines - 1
        if (wl(l) >= wl_work(l1) .and. wl(l) <= wl_work(l1+1)) then
          fac = log(wl(l)/wl_work(l1))/log(wl_work(l1+1)/wl_work(l1))
          no(l) = exp(log(n_work(l1)) &
            & + fac*log(n_work(l1+1)/n_work(l1)))
          if (k_work(l1) <= 0.0_dp .or. k_work(l1+1) <= 0.0_dp) then
            ko(l) = 0.0_dp
          else
            ko(l) = exp(log(k_work(l1)) &
              & + fac*log(k_work(l1+1)/k_work(l1)))
          end if
          no(l) = max(no(l),0.0_dp)
          ko(l) = max(ko(l),0.0_dp)
          exit
        end if
      end do
    end if

  end do

end subroutine interp_nk


subroutine Ross_mean(nwl, wl, temp, Vl, Vr, idx)
  use, intrinsic :: iso_fortran_env
  implicit none

  integer, parameter :: dp = REAL64
  real(kind=dp), parameter :: hp = 6.62607015e-27_dp ! erg s - Planck's constant
  real(kind=dp), parameter :: c_s = 2.99792458e10_dp ! cm s^-1 - Vacuum speed of light
  real(kind=dp), parameter :: kb = 1.380649e-16_dp ! erg K^-1 - Boltzmann's constant

  integer, intent(in) :: nwl, idx
  real(dp), dimension(nwl), intent(in) :: wl, Vl
  real(dp), intent(in) :: temp
  real(dp), intent(out) :: Vr

  integer :: l
  real(dp) :: top, bot,  x, Ross_w, Vw
  real(dp), dimension(nwl) :: dBdT, wl_cm

  !! Subroutine calculates Rosseland mean weighted value
  !! Returns Rosseland mean to all array values

  Ross_w = 0.0_dp
  Vw = 0.0_dp

  do l = 1, nwl
    wl_cm(l) = wl(l) * 1e-4_dp
    x = (hp * c_s) / (wl_cm(l) * kb * temp)
    top = 2.0_dp * hp**2 * c_s**3 * exp(x)
    bot = wl_cm(l)**6 * kb * temp**2 * (exp(x - 1.0_dp))**2
    dBdT(l) = top/bot

    Vw = Vw + 1.0_dp/Vl(l) * dBdT(l)
    Ross_w = Ross_w + dBdT(l)

  end do

  if (idx == 2) then
    Vr = 1.0_dp / (Vw/Ross_w)
    return
  end if

  top = trapz(wl(:),1.0_dp/Vl(:)*dBdT(:))
  bot = trapz(wl(:),dBdT(:))

  if (idx == 1)then
    Vr = 1.0_dp / (top/bot)
  end if

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
