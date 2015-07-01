! This F90 code computes the internal distances of polymer
! chains.

! Large integers are needed so compile with the following:
! gfortran -fdefault-integer-8 -O3 internal_distances.f90

program internal_distances

implicit none

integer*4 chains, monomers, monomers_half, total_particles
parameter (chains = 260)
parameter (monomers = 200)
parameter (monomers_half = monomers / 2)
parameter (total_particles = chains * monomers)
integer*4 id(total_particles)
real*8 x(total_particles), y(total_particles), z(total_particles)
real*8 x_buf(total_particles), y_buf(total_particles), z_buf(total_particles)
real*8 xp1, yp1, zp1, rsq
real*8 dsq(monomers_half), dsq_dsqave(monomers_half)
integer*8 ct(monomers_half)

integer*4 i, j, k, p1, p2, s

integer*8 f, istr, iend, incr
character(10) chrf10
character(80) prefix
character(80) flnm

character(80) chrf80
integer*8 timestep
integer*4 np
real*8 xmin, ymin, zmin
real*8 xmax, ymax, zmax
real*8 dx, dy, dz
real*8 rho, drho

ct = 0
dsq = 0.0D0

prefix = '../../rings10_250.'
istr = 1000000
iend = 2000000000
incr = 1000000

do f = istr, iend, incr
 
  write(chrf10, '(i10)') f
  flnm = trim(adjustl(prefix)) // trim(adjustl(chrf10))
  write(*,*) 'Working on file ', flnm

  open(10, file = flnm)
  read(10,*) chrf80
  read(10,*) timestep
  read(10,*) chrf80
  read(10,*) np
  read(10,*) chrf80
  read(10,*) xmin, xmax
  read(10,*) ymin, ymax
  read(10,*) zmin, zmax
  read(10,*) chrf80
  do i = 1, total_particles
    read(10,*) id(i), x_buf(i), y_buf(i), z_buf(i)
  enddo
  close(10)

  ! consistency checks
  if(timestep.ne.f) write(*,*) 'WARNING: Mismatch between timestep and file name: ', timestep, f
  if(np.ne.total_particles) write(*,*) 'WARNING: Mismatch between np and total_particles.'
  rho = np / ((xmax - xmin) * (ymax - ymin) * (zmax - zmin))
  drho = abs(1.0D0 - rho / 0.85D0)
  if(drho.gt.0.01D0) write(*,*) 'WARNING: Mismatch between densities: ', drho

  ! sort by id
  do i = 1, total_particles
    j = id(i)
    x(j) = x_buf(i)
    y(j) = y_buf(i)
    z(j) = z_buf(i)
  enddo

  ! loop over only first 10 chains (which are rings)
  do i = 1, 10
    do j = 1, monomers - 1
      p1 = (i - 1) * monomers + j
      xp1 = x(p1)
      yp1 = y(p1)
      zp1 = z(p1)
      do k = j + 1, monomers
        p2 = (i - 1) * monomers + k
        rsq = (xp1 - x(p2))**2 + (yp1 - y(p2))**2 + (zp1 - z(p2))**2
        s = k - j
        if(s.gt.monomers_half) s = monomers - s
        dsq(s) = dsq(s) + rsq
        ct(s) = ct(s) + 1
      enddo
    enddo
  enddo

enddo

do i = 1, monomers_half
  dsq(i) = dsq(i) / ct(i)
  dsq_dsqave(i) = dsq(i) * (1.0 / i + 1.0 / (monomers - i))
enddo

open(10, file = 'dsq.dat')
write(10,*) '# s  dsq(s)/sigma^2  dsq(s)/<dsq>'
do i = 1, monomers_half
  write(10,'(i5,f8.3,f7.3)') i, dsq(i), dsq_dsqave(i)
enddo
close(10)

end program
