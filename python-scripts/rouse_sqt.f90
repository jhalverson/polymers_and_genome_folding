! This F90 code computes S(q,t) according to the
! Rouse model for either a linear chain or ring. See
! Eqns. 36-39 in Tsolou et al., Macromolecules 43 (2010).

! Time values are stored in an array so that one can easily
! go between linear and logarithmic scales. The units for
! all times are tau (i.e., dt is in tau).

program rouse_sqt

implicit none

integer*4 upper_bound
parameter (upper_bound = 100)
real*8 time(0:upper_bound), Sqt(0:upper_bound)

integer*4 monomers, m, n
integer*4 i, p, p2, pcut
real*8 q, q26, phi
real*8 pi, pi2
real*8 dt, t, zeta, b2, kT, D6, fac
real*8 tau_linear
logical*4 linear, time_scale_linear

integer*4 start_exp, end_exp
character(20) outfile

linear = .false.
time_scale_linear = .false.
q = 0.2D0
monomers = 400

zeta = 43.0D0
b2 = 2.71D0
kT = 1.0D0
D6 = 6.0D0 * kT / (monomers * zeta)
pcut = 20

if(time_scale_linear) then
  ! linear time
  dt = 1.0D0
  do i = 0, upper_bound
    time(i) = i * dt
  enddo
else
  ! log scale
  start_exp = 3
  end_exp = 7
  do i = 0, upper_bound
    time(i) = 10.0D0**(start_exp + real(end_exp - start_exp) * i / upper_bound)
  enddo
endif

pi = 4.0D0 * atan(1.0D0)
pi2 = pi**2
q26 = q * q / 6.0D0
fac = 4.0D0 * monomers * b2 / pi2
tau_linear = zeta * b2 * monomers**2 / (3.0D0 * pi2 * kT)
if(linear) then
  write(*,*) 'linear q = ', q
  write(*,*) 'tau_linear =', tau_linear
else
  write(*,*) 'ring q = ', q
  write(*,*) 'tau_ring =', tau_linear / 4.0D0
endif

do i = 0, upper_bound
  t = time(i)
  Sqt(i) = 0.0D0
  if(linear) then

    ! linear
    do n = 1, monomers
      do m = 1, monomers
        phi = 0.0D0
        do p = 1, pcut
          p2 = p**2
          phi = phi + cos(p * pi * n / monomers) * cos(p * pi * m / monomers) * (1.0D0 - exp(-p2 * t / tau_linear)) / p2
        enddo
        phi = fac * phi + D6 * t + abs(n - m) * b2
        Sqt(i) = Sqt(i) + exp(-q26 * phi)
      enddo
    enddo

  else

    ! ring
    do n = 1, monomers
      do m = 1, monomers
        phi = 0.0D0
        do p = 2, pcut, 2
          p2 = p**2
          phi = phi + cos(p * pi * (n - m) / monomers) * (1.0D0 - exp(-p2 * t / tau_linear)) / p2
        enddo
        phi = fac * phi + D6 * t + abs(n - m) * (monomers - abs(n - m)) * b2 / monomers
        Sqt(i) = Sqt(i) + exp(-q26 * phi)
      enddo
    enddo

  endif
  Sqt(i) = Sqt(i) / monomers
enddo

outfile = 'Sqt_rouse_ring.dat'
if(linear) outfile = 'Sqt_rouse_linear.dat'

open(10, file=outfile)
if(linear) then
  write(10,*) '# linear, q = ', q
else
  write(10,*) '# ring, q = ', q
endif
write(10,*) '# pcut = ', pcut
write(10,*) '# time/tau  S(q,t)/S(q,0)  S(q,t)'
do i = 0, upper_bound
  write(10,'(f15.1,f15.6,f15.6)') time(i), Sqt(i) / Sqt(0), Sqt(i)
enddo
close(10)

end program
