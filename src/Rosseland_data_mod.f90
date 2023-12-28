module Rosseland_data_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  integer, parameter :: dp = REAL64

  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi

  real(dp), parameter :: hp = 6.62607015e-27_dp ! erg s - Planck's constant
  real(dp), parameter :: c_s = 2.99792458e10_dp ! cm s^-1 - Vacuum speed of light
  real(dp), parameter :: kb = 1.380649e-16_dp ! erg K^-1 - Boltzmann's constant

  integer :: nwl, na, nT
  real(dp), allocatable, dimension(:) :: a, T, wl

  logical :: c_flag
  integer :: nlines
  real(dp), allocatable, dimension(:) :: wl_ori, n_ori, k_ori, n_int, k_int

contains

  subroutine Ross_mean(nwl, wl, temp, Vl, Vr)
    implicit none

    integer, intent(in) :: nwl
    real(dp), dimension(nwl), intent(in) :: wl, Vl
    real(dp), intent(in) :: temp

    real(dp), intent(out) :: Vr

    integer :: l
    real(dp) :: top, bot,  xx, expx
    real(dp), dimension(nwl) :: dBdT, wl_cm

    !! Subroutine calculates Rosseland mean weighted value
    !! Returns Rosseland mean to all array values

    do l = 1, nwl

      wl_cm(l) = wl(l) * 1e-4_dp
      xx = (hp * c_s) / (wl_cm(l) * kb * temp)
      xx = min(xx, 35.0_dp) ! Avoid overfloating
      expx = exp(xx)
      top = 2.0_dp * hp**2 * c_s**3 * expx
      bot = wl_cm(l)**6 * kb * temp**2 * (expx - 1.0_dp)**2
      dBdT(l) = top/bot
      
    end do

    top = trapz(wl_cm(:),(1.0_dp/Vl(:))*dBdT(:))
    bot = trapz(wl_cm(:),dBdT(:))

    Vr = 1.0_dp / (top/bot)

  end subroutine Ross_mean

  subroutine read_RT_table(RT_path)
    implicit none

    character(len=100), intent(in) :: RT_path 

    integer :: uRT

    print*, 'Reading: ', trim(RT_path)
    open(newunit=uRT, file=trim(RT_path),action='read')
    read(uRT,*) na, nT
    allocate(a(na),T(nT))
    read(uRT,*) a(:)
    read(uRT,*) T(:)
    close(uRT)

  end subroutine read_RT_table

  subroutine read_nk_data(nk_path)
    implicit none

    character(len=100), intent(in) :: nk_path  

    integer :: unk, n

    ! Read in n,k constants for species - Note comment out/in lines for 'MgSiO3_2' lines (Xianyu's original data)
    print*, 'Reading: ', trim(nk_path)
    open(newunit=unk, file=trim(nk_path),action='read')
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

    close(unk)

  end subroutine read_nk_data

  subroutine read_wavelength_table(wl_path)
    implicit none

    character(len=100), intent(in) :: wl_path
    
    integer :: uwl, l

    ! Read in wavelength grid
    print*, 'Reading: ', trim(wl_path)
    open(newunit=uwl, file=trim(wl_path),action='read')
    read(uwl,*) nwl
    allocate(wl(nwl))
    do l = 1, nwl
      read(uwl,*) wl(l)
      !print*, l, wl(l)
    end do

    close(uwl)

  end subroutine read_wavelength_table

  subroutine interp_nk()
    implicit none

    integer :: l, l1
    real(dp) :: fac

    ! Interpolate n,k constants to wavelength grid
    allocate(n_int(nwl),k_int(nwl))

    do l = 1, nwl
      ! If required wavelength is less than availible data - keep constant
      if (wl(l) < wl_ori(1)) then
        n_int(l) = n_ori(1)
        k_int(l) = k_ori(1)
        n_int(l) = max(n_int(l),0.0_dp)
        k_int(l) = max(k_int(l),0.0_dp)
      ! If required wavelength is greater than availible data - extrapolate
      ! Non conducting: n is constant - k is linear decreasing
      ! Conducting: n and k are log-log extrapolated
      else if (wl(l) > wl_ori(nlines)) then
        if (c_flag .eqv. .False.) then
          n_int(l) = n_ori(nlines)
          k_int(l) = k_ori(nlines)*wl_ori(nlines)/wl(l)
          n_int(l) = max(n_int(l),0.0_dp)
          k_int(l) = max(k_int(l),0.0_dp)
        else if (c_flag .eqv. .True.) then
          do l1 = nlines,1,-1
            ! data can be noisy, so it's safer to use larger region to get the slope
            if (wl_ori(l1) < 0.7_dp*wl_ori(nlines)) then 
              exit                            
            end if                            
          enddo
          fac = log(wl(l)/wl_ori(nlines))/log(wl_ori(l1)/wl_ori(nlines))
          n_int(l) = exp(log(n_ori(nlines)) &
            & + fac*log(n_ori(l1)/n_ori(nlines)))
          k_int(l) = exp(log(k_ori(nlines)) &
            & + fac*log(k_ori(l1)/k_ori(nlines)))
          n_int(l) = max(n_int(l),0.0_dp)
          k_int(l) = max(k_int(l),0.0_dp)
        end if
      ! Data is availible in the required wavelength range - log-log interpolation
      else
        ! Loop across work arrays untill straddle point is point then interpolate
        do l1 = 1, nlines - 1
          if (wl(l) >= wl_ori(l1) .and. wl(l) <= wl_ori(l1+1)) then
            fac = log(wl(l)/wl_ori(l1))/log(wl_ori(l1+1)/wl_ori(l1))
            n_int(l) = exp(log(n_ori(l1)) &
              & + fac*log(n_ori(l1+1)/n_ori(l1)))
            if (k_ori(l1) <= 0.0_dp .or. k_ori(l1+1) <= 0.0_dp) then
              k_int(l) = 0.0_dp
            else
              k_int(l) = exp(log(k_ori(l1)) &
                & + fac*log(k_ori(l1+1)/k_ori(l1)))
            end if
            n_int(l) = max(n_int(l),0.0_dp)
            k_int(l) = max(k_int(l),0.0_dp)
            exit
          end if
        end do
      end if

    end do

  end subroutine interp_nk

  function trapz(xx, yy) result(rr)
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

end module Rosseland_data_mod