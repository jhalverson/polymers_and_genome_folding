! This F90 code computes g1(t) and g3(t) according to the
! Rouse model for either a linear chain or ring.

! Time values are stored in an array so that one can easily
! go between linear and logarithm scales. The units for
! all times are tau (i.e., dt is in tau).

! m is a continuous variable between 0 and N; when m=1/2 the
! linear and ring cases give the same results for g1. They
! always give the same result for g3.

program rouse_g1_g3

implicit none

integer*4 upper_bound
parameter (upper_bound = 100)
real*8 time(0:upper_bound), g1(0:upper_bound)

integer*4 N, m
integer*4 i, p, p2, pcut
real*8 pi, pi2
real*8 dt, t, zeta, b2, kT, D6, cos2, rsum
real*8 tau_linear, tau_ring
logical*4 linear, time_scale_linear

integer*4 start_exp, end_exp
real*8 step_exp

linear = .false.
time_scale_linear = .false.
N = 100
m = 0
zeta = 43.0D0
b2 = 2.71D0
kT = 1.0D0
D6 = 6.0D0 * kT / (N * zeta)

pcut = 10000

if(time_scale_linear) then
  ! linear time
  dt = 1.0D0
  do i = 0, upper_bound
    time(i) = i * dt
  enddo
else
  ! log scale
  start_exp = 3
  end_exp = 6
  step_exp = real(end_exp - start_exp) / upper_bound
  do i = 0, upper_bound
    time(i) = 10.0D0**(start_exp + step_exp * i)
  enddo
endif

pi = 4.0D0 * atan(1.0D0)
pi2 = pi**2
tau_linear = zeta * b2 * N**2 / (3.0D0 * pi2 * kT)
tau_ring = tau_linear / 4.0D0
write(*,*) 'tau_linear =', tau_linear
write(*,*) 'tau_ring =', tau_ring

do i = 0, upper_bound
  t = time(i)
  rsum = 0.0D0
  if(linear) then
    ! linear
    do p = 1, pcut
      p2 = p**2
      cos2 = (cos(pi * p * m / N))**2
      rsum = rsum + cos2 * (1.0D0 - exp(-p2 * t / tau_linear)) / p2
    enddo
  else
    ! ring
    do p = 2, pcut, 2
      p2 = p**2
      rsum = rsum + (1.0D0 - exp(-p2 * t / tau_linear)) / p2
    enddo
  endif
  g1(i) = D6 * t + (4.0D0 * N * b2 / pi2) * rsum
enddo

open(10, file='msd_rouse.dat')
write(10,*) '# time/tau  g1/sigma**2  g3/sigma**2'
do i = 0, upper_bound
  write(10,'(f15.6,f15.6,f15.6)') time(i), g1(i), D6 * time(i)
enddo
close(10)

end program
