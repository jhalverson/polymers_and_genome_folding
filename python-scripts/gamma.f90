! This F90 code computes probability that monomers i and j of the
! same chain are separated by nij bonds and within rc

! Large integers are needed so compile with the following:
! gfortran -fdefault-integer-8 -O3 gamma.f90

program gamma

implicit none

integer*4 chains, monomers, monomers_half, total_particles
parameter (chains = 200)
parameter (monomers = 200)
parameter (monomers_half = monomers / 2)
parameter (total_particles = chains * monomers)
real*8 x(total_particles), y(total_particles), z(total_particles)

integer*4 trials(monomers_half)
integer*4 successes(monomers_half)
integer*4 nij
real*8 Pij

integer*4 i, j, k
integer*4 p1, p2
integer*4 ix, iy, iz

real*8 rho
real*8 L, L2
real*8 rc, rcsq
real*8 rij, rijsq
real*8 x1, y1, z1
real*8 x2, y2, z2
real*8 xij, yij, zij

integer*8 f, istr, iend, incr
character(10) chrf10
character(80) prefix
character(80) flnm

rho = 0.85D0
L = (total_particles / rho)**(1.0 / 3.0)
L2 = 0.5D0 * L
rc = 2.0D0**(1.0D0 / 6.0D0)
rcsq = rc**2

trials = 0
successes = 0

prefix = '../../200/'
istr = 1000000000
iend = 2000000000
incr = 5000000

do f = istr, iend, incr

  write(chrf10, '(i10)') f
  flnm = trim(adjustl(prefix)) // trim(adjustl(chrf10)) // '.pos'
  write(*,*) 'Working on file ', flnm
 
  ! read in coordinates
  open(10, file = flnm)
  do i = 1, total_particles
    read(10,*) x(i), y(i), z(i)
  enddo
  close(10)
 
  if(f.eq.istr) then
 
    write(*,*) "density = ", total_particles / L**3
    write(*,*) "L = ", L, ", L2 = ", L2
    write(*,*) "total_particles = ", total_particles, "  monomers = ", monomers
    write(*,*) "r(1) = ", x(1), y(1), z(1)
    write(*,*) "r(end) = ", x(total_particles), y(total_particles), z(total_particles)
 
  endif

  do i = 1, chains

    do j = 1, monomers - 1

      do k = j + 1, monomers

        p1 = (i - 1) * monomers + j
        p2 = (i - 1) * monomers + k

        ! number of bonds separating p1 and p2
        nij = min(k - j, monomers - (k - j))

        ! translate p1 to central box
        x1 = x(p1)
        y1 = y(p1)
        z1 = z(p1)

        ix = int(x1 / L)
        iy = int(y1 / L)
        iz = int(z1 / L)

        if(x1.lt.0.0D0) ix = ix - 1
        if(y1.lt.0.0D0) iy = iy - 1
        if(z1.lt.0.0D0) iz = iz - 1

        x1 = x1 - ix * L
        y1 = y1 - iy * L
        z1 = z1 - iz * L

        ! translate p2 to central box
        x2 = x(p2)
        y2 = y(p2)
        z2 = z(p2)

        ix = int(x2 / L)
        iy = int(y2 / L)
        iz = int(z2 / L)

        if(x2.lt.0.0D0) ix = ix - 1
        if(y2.lt.0.0D0) iy = iy - 1
        if(z2.lt.0.0D0) iz = iz - 1

        x2 = x2 - ix * L
        y2 = y2 - iy * L
        z2 = z2 - iz * L

        ! compute minimum image distance
        xij = x1 - x2
        yij = y1 - y2
        zij = z1 - z2

        if(xij.gt. L2) xij = xij - L
        if(xij.lt.-L2) xij = xij + L
        if(yij.gt. L2) yij = yij - L
        if(yij.lt.-L2) yij = yij + L
        if(zij.gt. L2) zij = zij - L
        if(zij.lt.-L2) zij = zij + L

        rijsq = xij**2 + yij**2 + zij**2
        trials(nij) = trials(nij) + 2
        if(rijsq.le.rcsq) successes(nij) = successes(nij) + 2

      enddo

    enddo

  enddo
 
enddo

open(11, file="pij.dat")
do i = 1, monomers_half
  Pij = real(successes(i)) / trials(i)
  write(11,'(i4,i11,i11,f9.6)') i, trials(i), successes(i), Pij
enddo
close(11)

end program
