! This F90 code computes the self-excluded density profile
! or self-excluded RDF of a polymer melt.

! Large integers are needed so compile with the following:
! gfortran -fdefault-integer-8 -O3 ex_rdf.f90

program ex_rdf

implicit none

integer*4 chains, monomers, total_particles, ibin, ibinmax
real*8 L, dr
parameter (chains = 200)
parameter (monomers = 200)
parameter (total_particles = chains * monomers)
parameter (L = 36.103310D0)
parameter (dr = 0.25D0)
parameter (ibinmax = 0.5 * L / dr + 1)
integer*4 ihist(ibinmax)
real*8 rhist(ibinmax)
real*8 x(total_particles), y(total_particles), z(total_particles)

integer*4 i, j, p, jchain, files
integer*4 ix, iy, iz

real*8 rij, rijsq, L2sq, L2
real*8 xi, yi, zi
real*8 xj, yj, zj
real*8 xij, yij, zij
real*8 xcm, ycm, zcm
real*8 rho, rho_, pi, fac, rinner, router

integer*8 f, istr, iend, incr
character(10) chrf10
character(80) prefix
character(80) flnm

! variables for finding bin with highest density
integer*4 density_bins
parameter (density_bins = 1000)
integer*4 ct(density_bins)
real*8 xcen(density_bins), ycen(density_bins), zcen(density_bins)
real*8 xmin, ymin, zmin
real*8 xmax, ymax, zmax
real*8 xtmp, ytmp, ztmp
real*8 db
integer*4 nx, ny, nz
integer*4 num_max, num_max_ave, j_max
integer*4 dbins
real*8 dist_from_mass_center

! use mass center of central chain or center of bin with max density
logical*4 mass_center

! output is rdf or density
logical*4 rdf

rho_ = abs(1.0D0 - (total_particles / L**3) / 0.85D0)
if(rho_.gt.0.01D0) write(*,*) 'WARNING: Mismatch between densities: ', rho_, rho

! bin size for density calculation (all lengths in sigma)
db = 2.5D0

L2 = 0.5D0 * L
L2sq = L2**2
rho = 0.85D0

files = 0
ihist = 0
num_max_ave = 0
dist_from_mass_center = 0.0D0

prefix = '../../200/'
istr = 500000
iend = 2000000000
incr = 50000000

mass_center = .false.
rdf = .false.

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

    xcm = 0.0D0
    ycm = 0.0D0
    zcm = 0.0D0

    do j = 1, monomers

      p = (i - 1) * monomers + j
      xcm = xcm + x(p)
      ycm = ycm + y(p)
      zcm = zcm + z(p)

    enddo

    xi = xcm / monomers
    yi = ycm / monomers
    zi = zcm / monomers

    if(mass_center) then

      ix = int(xi / L)
      iy = int(yi / L)
      iz = int(zi / L)

      if(xi.lt.0.0D0) ix = ix - 1
      if(yi.lt.0.0D0) iy = iy - 1
      if(zi.lt.0.0D0) iz = iz - 1

      xi = xi - ix * L
      yi = yi - iy * L
      zi = zi - iz * L

    else

      xmin = 1e9
      ymin = 1e9
      zmin = 1e9
      xmax = -1e9
      ymax = -1e9
      zmax = -1e9

      do j = 1, monomers

        p = (i - 1) * monomers + j
        if(x(p) < xmin) xmin = x(p)
        if(y(p) < ymin) ymin = y(p)
        if(z(p) < zmin) zmin = z(p)

        if(x(p) > xmax) xmax = x(p)
        if(y(p) > ymax) ymax = y(p)
        if(z(p) > zmax) zmax = z(p)

      enddo

      nx = int((xmax - xmin) / db) + 1
      ny = int((ymax - ymin) / db) + 1
      nz = int((zmax - zmin) / db) + 1
      dbins = nx * ny * nz

      if(dbins.gt.density_bins) write(*,*) "ERROR: density_bins too small: ", dbins, density_bins

      ! define bin centers
      do ix = 1, nx
        do iy = 1, ny
          do iz = 1, nz
            xtmp = (ix - 1) * db + 0.5D0 * db
            ytmp = (iy - 1) * db + 0.5D0 * db
            ztmp = (iz - 1) * db + 0.5D0 * db
            ibin = 1 + int(xtmp / db) + int(ytmp / db) * nx + int(ztmp / db) * nx * ny
            xcen(ibin) = xmin + xtmp
            ycen(ibin) = ymin + ytmp
            zcen(ibin) = zmin + ztmp
          enddo
        enddo
      enddo

      ct = 0
      do j = 1, monomers

        p = (i - 1) * monomers + j
        ibin = 1 + int((x(p) - xmin) / db) + &
                   int((y(p) - ymin) / db) * nx + &
                   int((z(p) - zmin) / db) * nx * ny
        ct(ibin) = ct(ibin) + 1

      enddo

      if(sum(ct).ne.monomers) write(*,*) "ERROR: incorrect sum", sum(ct), monomers

      ! find bin with maximum density
      num_max = 0
      do j = 1, dbins
        if(ct(j) > num_max) then
          num_max = ct(j)
          j_max = j
        endif
      enddo

      num_max_ave = num_max_ave + num_max
      dist_from_mass_center = dist_from_mass_center + sqrt((xi - xcen(j_max))**2 + &
                                                           (yi - ycen(j_max))**2 + &
                                                           (zi - zcen(j_max))**2)

      ! assign and shift center
      xi = xcen(j_max)
      yi = ycen(j_max)
      zi = zcen(j_max)

      ix = int(xi / L)
      iy = int(yi / L)
      iz = int(zi / L)

      if(xi.lt.0.0D0) ix = ix - 1
      if(yi.lt.0.0D0) iy = iy - 1
      if(zi.lt.0.0D0) iz = iz - 1

      xi = xi - ix * L
      yi = yi - iy * L
      zi = zi - iz * L

    endif

    do j = 1, total_particles

      jchain = (j - 1) / monomers + 1

      if(i.ne.jchain) then

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

        if(xij.gt. L2) xij = xij - L
        if(xij.lt.-L2) xij = xij + L
        if(yij.gt. L2) yij = yij - L
        if(yij.lt.-L2) yij = yij + L
        if(zij.gt. L2) zij = zij - L
        if(zij.lt.-L2) zij = zij + L

        rijsq = xij**2 + yij**2 + zij**2
        if(rijsq.le.L2sq) then

          rij = sqrt(rijsq)
          ibin = int(rij / dr) + 1
          if(ibin.gt.ibinmax) write(*,*) "ERROR: ibin"
          ihist(ibin) = ihist(ibin) + 1

        endif

      endif

    enddo

  enddo
 
  files = files + 1
 
enddo

pi = 4.0D0 * atan(1.0D0)
fac = 4.0D0 * pi * rho / 3.0D0
do i = 1, ibinmax
  rinner = (i - 1) * dr
  router = rinner + dr
  rhist(i) = ihist(i) / (fac * (router**3 - rinner**3)) / (chains * files)
enddo

if(.not.rdf) then
  do i = 1, ibinmax
    rhist(i) = rho * rhist(i)
  enddo
endif

open(11, file="rdfex.dat")
do i = 1, ibinmax - 1
  write(11,'(f8.3,f8.3,f8.3,f8.3)') (i - 1) * dr, i * dr, i * dr - 0.5D0 * dr, rhist(i)
enddo
close(11)

write(*,*) "<num_max> = ", real(num_max_ave) / (chains * files)
write(*,*) "<rho_max> = ", real(num_max_ave) / (chains * files * db**3)
write(*,*) "<|r_CM - r_max_density|> = ", dist_from_mass_center / (chains * files)

end program
