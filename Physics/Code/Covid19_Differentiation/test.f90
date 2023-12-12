  implicit none
  integer, parameter:: N_max = 1000
  double precision, parameter :: dx_min = 1.0d-2, dx_max=1.0d-1
  double precision, parameter :: x_min = -3.0, x_max = 3.0
  double precision f, df_dx, dx, x, a, actual
  double precision tmp1, tmp2, tmp3, tmp4
  integer n, f_loop

  open(file='df_dx_central.dat', unit=11)

    dx = .01
    do n=1, N_max
        x = x_min+(x_max-x_min)*float(n-1)/float(N_max-1) ! x at place where u want to find der
        !fourth order/ app
        tmp1 = f(x+2.0*dx)
        tmp2 = f(x+dx)
        tmp3 = f(x-dx)
        tmp4 = f(x-2.0*dx)
        df_dx = (-tmp1+8.0*tmp2-8.0*tmp3+tmp4)/(12.0*dx)
        ! analy
        actual = a(x)
        write(11, *) x, df_dx, actual
        enddo
    enddo

  close(unit=11)
end program