! This F90 code computes the ring-ring, ring-linear, and
! linear-linear COM-COM rdf.

! Large integers are needed so compile with the following:
! gfortran -fdefault-integer-8 -O3 rdf_com_blend.f90

program rdf_com_blend

implicit none

integer*4 M_rings, M_linear, chains, monomers, total_particles, ibinmax
real*8 L, dr
parameter (M_rings = 200)
parameter (M_linear = 26)
parameter (monomers = 200)
parameter (chains = M_rings + M_linear)
parameter (total_particles = chains * monomers)
parameter (L = 37.60450D0)
parameter (dr = 0.5D0)
parameter (ibinmax = 0.5D0 * L / dr + 1)
integer*4 id(total_particles)
real*8 x(total_particles), y(total_particles), z(total_particles)
real*8 x_buf(total_particles), y_buf(total_particles), z_buf(total_particles)

integer*4 i, j, k, p, files
integer*8 f, istr, iend, incr

character(10) chrf10
character(80) prefix
character(80) flnm

integer*4 ix, iy, iz
real*8 xcm, ycm, zcm
real*8 xi, yi, zi
real*8 xj, yj, zj
real*8 xij, yij, zij
real*8 rijsq, rij
real*8 L_half, L2sq

integer*4 ibin
integer*4 ihist_ring_ring(ibinmax)
integer*4 ihist_linear_linear(ibinmax)
integer*4 ihist_ring_linear(ibinmax)
real*8 pi, fac, rinner, router
real*8 rhist(ibinmax)

logical*4 lammps

character(80) chrf80
integer*8 timestep
integer*4 np
real*8 xmin, ymin, zmin
real*8 xmax, ymax, zmax
real*8 dx, dy, dz
real*8 rho, drho

L_half = 0.5D0 * L
L2sq = L_half**2

files = 0
ihist_ring_ring = 0
ihist_linear_linear = 0
ihist_ring_linear = 0

istr = 1000000000
iend = 4000000000
incr = 100000

prefix = '../configs/rings200_26_glob.'
lammps = .true.

do f = istr, iend, incr

  if(lammps) then

    ! LAMMPS format
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

  else
 
    ! pos format
    write(chrf10, '(i10)') f
    flnm = trim(adjustl(prefix)) // trim(adjustl(chrf10)) // '.pos'
    write(*,*) 'Working on file ', flnm

    open(10, file = flnm)
    do i = 1, total_particles
      read(10,*) x(i), y(i), z(i)
    enddo
    close(10)

  endif

  if(f.eq.istr) then
 
    write(*,*) "computed rho = ", total_particles / L**3
    write(*,*) "L = ", L, ", L_half = ", L_half
    write(*,*) "M_rings = ", M_rings, "  M_linear = ", M_linear
    write(*,*) "total_particles = ", total_particles, "  monomers = ", monomers
    write(*,*) "r(1) = ", x(1), y(1), z(1)
    write(*,*) "r(end) = ", x(total_particles), y(total_particles), z(total_particles)
 
  endif

  do i = 1, chains - 1

    xcm = 0.0D0
    ycm = 0.0D0
    zcm = 0.0D0

    do k = 1, monomers

      p = (i - 1) * monomers + k
      xcm = xcm + x(p)
      ycm = ycm + y(p)
      zcm = zcm + z(p)

    enddo

    xi = xcm / monomers
    yi = ycm / monomers
    zi = zcm / monomers

    ix = int(xi / L)
    iy = int(yi / L)
    iz = int(zi / L)

    if(xi.lt.0.0D0) ix = ix - 1
    if(yi.lt.0.0D0) iy = iy - 1
    if(zi.lt.0.0D0) iz = iz - 1

    xi = xi - ix * L
    yi = yi - iy * L
    zi = zi - iz * L

    do j = i + 1, chains

      xcm = 0.0D0
      ycm = 0.0D0
      zcm = 0.0D0

      do k = 1, monomers

        p = (j - 1) * monomers + k
        xcm = xcm + x(p)
        ycm = ycm + y(p)
        zcm = zcm + z(p)

      enddo

      xj = xcm / monomers
      yj = ycm / monomers
      zj = zcm / monomers

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
      if(rijsq.le.L2sq) then

        rij = sqrt(rijsq)
        ibin = int(rij / dr) + 1
        if(ibin.gt.ibinmax.or.ibin.lt.1) write(*,*) "ERROR: ibin ibinmax rij", ibin, ibinmax, rij, i, j

        if(i.le.M_rings.and.j.le.M_rings) then
          ihist_ring_ring(ibin) = ihist_ring_ring(ibin) + 2
        elseif(i.gt.M_rings.and.j.gt.M_rings) then
          ihist_linear_linear(ibin) = ihist_linear_linear(ibin) + 2
        else
          ihist_ring_linear(ibin) = ihist_ring_linear(ibin) + 1
        endif

      endif

    enddo

  enddo
 
  files = files + 1
 
enddo

pi = 4.0D0 * atan(1.0D0)
if(M_rings.gt.1) then

  ! normalize
  fac = 4.0 * pi / (3.0D0 * L**3)
  do i = 1, ibinmax
    rinner = (i - 1) * dr
    router = rinner + dr
    rhist(i) = ihist_ring_ring(i) / (fac * (router**3 - rinner**3) *  M_rings * (M_rings - 1) * files)
  enddo
 
  ! write to file
  open(11, file="rdf_ring_ring.dat")
  do i = 1, ibinmax - 1
    write(11,'(f8.3,f8.3,f8.3,f8.3)') (i - 1) * dr, i * dr, i * dr - 0.5 * dr, rhist(i)
  enddo
  close(11)

endif

if(M_linear.gt.1) then

  ! normalize
  fac = 4.0 * pi / (3.0D0 * L**3)
  do i = 1, ibinmax
    rinner = (i - 1) * dr
    router = rinner + dr
    rhist(i) = ihist_linear_linear(i) / (fac * (router**3 - rinner**3) * M_linear * (M_linear - 1) * files)
  enddo

  ! write to file
  open(11, file="rdf_linear_linear.dat")
  do i = 1, ibinmax - 1
    write(11,'(f8.3,f8.3,f8.3,f8.3)') (i - 1) * dr, i * dr, i * dr - 0.5 * dr, rhist(i)
  enddo
  close(11)

endif

if(M_rings.gt.0.and.M_linear.gt.0) then

  ! normalize
  fac = 4.0 * pi * (M_rings / L**3) / 3.0D0
  do i = 1, ibinmax
    rinner = (i - 1) * dr
    router = rinner + dr
    rhist(i) = ihist_ring_linear(i) / (fac * (router**3 - rinner**3)) / (M_linear * files)
  enddo

  ! write to file
  open(11, file="rdf_ring_linear.dat")
  do i = 1, ibinmax - 1
    write(11,'(f8.3,f8.3,f8.3,f8.3)') (i - 1) * dr, i * dr, i * dr - 0.5 * dr, rhist(i)
  enddo
  close(11)

endif

end program
