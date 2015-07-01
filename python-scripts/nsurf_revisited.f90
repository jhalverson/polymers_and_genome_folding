program nsurf_revisited

! The number of interchain interactions per chain was computed
! previously. Here we compute the number of monomers per chain
! that have at least one interaction with another chain where
! an interaction means a separation distance of less than rc

! Large integers are needed so compile with the following:
! gfortran -fdefault-integer-8 -O3 nsurf_revisited.f90

implicit none

integer*4 chains, monomers, total_particles
parameter (chains = 200)
parameter (monomers = 100)
parameter (total_particles = chains * monomers)
real*8 x(total_particles), y(total_particles), z(total_particles)

integer*4 i, j, k, p
integer*4 ichain, jchain
integer*4 ix, iy, iz
integer*4 nsurf, nsurf_total
integer*4 files

integer*8 f, istr, iend, incr
character(10) chrf10
character(80) prefix
character(80) flnm

real*8 L, L_half
real*8 rijsq, rc, rcsq
real*8 xi, yi, zi
real*8 xj, yj, zj
real*8 xij, yij, zij

logical*4 surface_monomer

! cutoff
rc = (10.0D0 / 10.0D0) * 2.0D0**(1.0D0/6.0D0)
rcsq = rc**2

L = (total_particles / 0.85D0)**(1.0 / 3.0)
L_half = 0.5D0 * L

prefix = '/data/pckr120/halvers/rings_data/200_0/stripped/'
prefix = '/data/pckr117/halvers/100_0/stripped/'
istr = 1000000000
iend = 2000000000
incr = 100000000

nsurf_total = 0
files = 0

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
    write(*,*) "L = ", L, ", L_half = ", L_half, ", rcsq = ", rcsq, ", rc = ", rc
    write(*,*) "total_particles = ", total_particles, "  monomers = ", monomers
    write(*,*) "r(1) = ", x(1), y(1), z(1)
    write(*,*) "r(end) = ", x(total_particles), y(total_particles), z(total_particles)
 
  endif

  nsurf = 0
  do k = 1, chains

    do i = 1, monomers

      p = (k - 1) * monomers + i

      xi = x(p)
      yi = y(p)
      zi = z(p)

      ix = int(xi / L)
      iy = int(yi / L)
      iz = int(zi / L)

      if(xi.lt.0.0D0) ix = ix - 1
      if(yi.lt.0.0D0) iy = iy - 1
      if(zi.lt.0.0D0) iz = iz - 1

      xi = xi - ix * L
      yi = yi - iy * L
      zi = zi - iz * L

      surface_monomer = .false.

      do j = 1, total_particles

        jchain = (j - 1) / monomers + 1

        if(k.ne.jchain) then

          xj = x(j)
          yj = y(j)
          zj = z(j)

          ix = int(xj / L)
          iy = int(yj / L)
          iz = int(zj / L)

          if(xj.lt.0.0D0) ix = ix - 1
          if(yj.lt.0.0D0) iy = iy - 1
          if(zj.lt.0.0D0) iz = iz - 1
          
          xj = xj - ix * L
          yj = yj - iy * L
          zj = zj - iz * L

          xij = xi - xj
          yij = yi - yj
          zij = zi - zj

          if(xij.gt. L_half) xij = xij - L
          if(xij.lt.-L_half) xij = xij + L
          if(yij.gt. L_half) yij = yij - L
          if(yij.lt.-L_half) yij = yij + L
          if(zij.gt. L_half) zij = zij - L
          if(zij.lt.-L_half) zij = zij + L

          rijsq = xij**2 + yij**2 + zij**2
          if(rijsq.le.rcsq) then
            surface_monomer = .true.
            exit
          endif

        endif

      enddo

      if(surface_monomer) nsurf = nsurf + 1

    enddo

  enddo

  files = files + 1
  nsurf_total = nsurf_total + nsurf
 
  write(*,*) "n_surf = ", real(nsurf) / chains
 
enddo

write(*,*) "final nsurf = ", real(nsurf_total) / (chains * files)
 
end program
