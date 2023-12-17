subroutine Rosseland_lognorm(sp)
  use Rosseland_data_mod
  use lxmie_mod, only : lxmie
  implicit none

  character(len=100), intent(in) :: sp

  integer :: aa, tt, l, u_k, u_a, u_g

  real(dp), dimension(:), allocatable :: kext_l, a_l, g_l
  real(dp), dimension(:, :), allocatable :: Ross_kext, Ross_a, Ross_g

  complex(dp) :: ri
  real(dp) :: x, xsec, q_ext, q_sca, q_abs, g


  !! Perform Mie calculations


  !! allocate work arrays for each wavelength
  allocate(kext_l(nwl), a_l(nwl), g_l(nwl))

  !! allocate end Rosseland mean arrays
  allocate(Ross_kext(na,nT),Ross_a(na,nT),Ross_g(na,nT))


  do aa = 1, na

    print*, a(aa), aa, na
    xsec = pi * (a(aa)* 1e-4_dp)**2

    do tt = 1, nt

      do l = 1, nwl

        ri = cmplx(n_int(l),-k_int(l),dp)
        x = (twopi * a(aa))/wl(l)

        call lxmie(ri, x, q_ext, q_sca, q_abs, g)

        kext_l(l) = xsec * q_ext
        a_l(l) = q_sca/q_ext
        g_l(l) = max(g,1.0e-12_dp)
        g_l(l) = min(g_l(l), 1.0_dp)

        !print*, a(aa), wl(l), x, q_ext, q_sca, g

      end do

      call Ross_mean(nwl, wl(:),T(tt),kext_l(:),Ross_kext(aa,tt))
      call Ross_mean(nwl, wl(:),T(tt),a_l(:),Ross_a(aa,tt))
      call Ross_mean(nwl, wl(:),T(tt),g_l(:),Ross_g(aa,tt))

    end do
  end do


  ! Output table - go radius outer loop, temperature inner loop
  ! Output Qext - Qsca and g

  open(newunit=u_k, file='results_Rosseland_lognorm/'//trim(sp)//'_kext.txt',action='readwrite')
  write(u_k,*) na, nT
  write(u_k,*) a(:)
  write(u_k,*) T(:)
  do aa = 1, na
    write(u_k,*) (real(Ross_kext(aa,tt)), tt = 1, nT)
    call flush(u_k)
  end do

  close(u_k)

  open(newunit=u_a, file='results_Rosseland_lognorm/'//trim(sp)//'_a.txt',action='readwrite')
  write(u_a,*) na, nt
  write(u_a,*) a(:)
  write(u_a,*) T(:)
  do aa = 1, na
    write(u_a,*) (real(Ross_a(aa,tt)), tt = 1, nT)
    call flush(u_a)
  end do

  close(u_a)

  open(newunit=u_g, file='results_Rosseland_lognorm/'//trim(sp)//'_g.txt',action='readwrite')
  write(u_g,*) na, nt
  write(u_g,*) a(:)
  write(u_g,*) T(:)
  do aa = 1, na
    write(u_g,*) (real(Ross_g(aa,tt)), tt = 1, nT)
    call flush(u_g)
  end do

  close(u_g)


end subroutine Rosseland_lognorm
