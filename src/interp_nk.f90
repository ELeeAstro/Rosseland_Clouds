subroutine interp_nk(conducting, nwl, n_lines, wl, wl_work, n_work, k_work, no, ko)
  use, intrinsic :: iso_fortran_env
  implicit none

  integer, parameter :: dp = REAL64

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