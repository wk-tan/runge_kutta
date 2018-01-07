!/////////////////////////////////////////////////////////////////////
module setup
  real, parameter :: tstep = 0.01
  integer, parameter :: nsteps = 4000
  integer, parameter :: n = 2
  real, dimension(1:n) :: x
  real :: t
end module setup
!/////////////////////////////////////////////////////////////////////
module parameter
  real, parameter :: a  = 1.0
  real, parameter :: b  = 2.0
  real, parameter :: c  = 2.0
  real, parameter :: d  = 1.0
end module parameter
!/////////////////////////////////////////////////////////////////////
program main
  use setup
  use parameter
  implicit none

  integer :: i

  x(1) = 5.0
  x(2) = 3.0
  
  open(6, file = 'output.txt')
  write(6,9) 0, x(1:n)
  do i = 1, nsteps
     t = i*tstep
     call rk4
     write(6,9) t, x(1:n)
  end do
  close(6)

  call gnuplot
  call system('gnuplot -p plot.gnu')
9 format(f10.2, 20f10.4)
end program main
!/////////////////////////////////////////////////////////////////////
subroutine rk4
  use setup
  implicit none

  real :: h
  real, dimension(1:n) :: k1, k2, k3, k4
  real, dimension(1:n) :: dxdt

  h = tstep/2.0

  call odesystem(t, n, x, dxdt)
  k1(1:n) = tstep*dxdt(1:n)

  call odesystem(t + h, n, x + k1/2.0, dxdt)
  k2(1:n) = tstep*dxdt(1:n)

  call odesystem(t + h, n, x + k2/2.0, dxdt)
  k3(1:n) = tstep*dxdt(1:n)

  call odesystem(t + h, n, x + k3, dxdt)
  k4(1:n) = tstep*dxdt(1:n)

  x(1:n) = x(1:n) + (k1(1:n) + 2.0*(k2(1:n) + k3(1:n)) + k4(1:n))/6.0
end subroutine rk4
!/////////////////////////////////////////////////////////////////////
subroutine odesystem(t, n, x, dxdt)
  use parameter
  implicit none
  real :: t
  integer :: n
  real, dimension(1:n) :: x, dxdt
  
  dxdt(1) = a*x(2)*x(1) - b*x(1)
  dxdt(2) = c*x(2) - d*x(2)*x(1)
end subroutine odesystem
!/////////////////////////////////////////////////////////////////////
subroutine gnuplot
  open(66, file = 'plot.gnu')
  write(66,*) "reset"
  write(66,*) 'set terminal pngcairo dashed enhanced size 480,360',   &
       &      "font 'Arial,12' fontscale 1.0"
  write(66,*) "set encoding utf8"
  write(66,*) "set output 'output.png'"
  write(66,*) 'set xlabel "{/Arial:Bold Time}"'
  write(66,*) 'set ylabel "{/Arial:Bold Variable}"'
  write(66,*) 'set key bmargin center horizontal height 1.0 ',        &
       &      "font 'Arial-Bold,10'"
  write(66,*) "plot 'output.txt' u 1:2 lt 1 lw 2.5 lc ", '"red" ',    &
       &      'title "x_1" w l, \'
  write(66,*) "     'output.txt' u 1:3 lt 1 lw 2.5 lc ", '"blue" ',   &
       &      'title "x_2" w l'
  close(66)
end subroutine gnuplot
!/////////////////////////////////////////////////////////////////////
