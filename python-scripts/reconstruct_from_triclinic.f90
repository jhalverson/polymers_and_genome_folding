! This F90 code reconstructs chains from LAMMPS sllod simulations where
! the shear direction is x and the velocity changes in the z-direction.
! It assumes a dump command with "id x y z ix iy iz" although the image
! flags are not used. Note the x, y, z coordinates may be slightly outside
! of the central triclinic box if the particles were written on a time
! step where the neighbor lists were not built.

! By setting the flags switch to false a dump of id xu yu zu may be
! used. The code now wraps the particles in the central box before
! attempting to reconstruct the chains.

! To reconstruct the chains the code starts with the first atom of
! each chain. The position of the second atom in the central cell and its
! sh images are considered as candidates with the selected position being
! that with the minimum separation from the previously accepted atom. The
! chains are then wrapped by their COM position to a orthorhombic box - it
! is understood that the image boxes in the z-direction are shifted by ix*xz.
! One could modify the code to maintain the angles if desired.

! This code was checked by reconstructing the chains and then computing
! the total potential energy and comparing that with the simulation value.
! It was necessary to use a triclinic box with the same xz value in the
! LAMMPS data file otherwise LAMMPS assumes an orthogonal box and the
! boundary conditions are wrong. This leads to two particles being very
! close and hence the simulation crashes immediately. Useful files:
! ~/research/clc/lammps_triclinic/dec23_triclinic_wrapped_to_orthogonal/configs

! Large integers are needed so compile with the following:
! gfortran -fdefault-integer-8 -O2 reconstruct_from_triclinic.f90

program reconstruct_from_triclinic

implicit none

integer*4 chains, monomers, total_particles
parameter (chains = 200)
parameter (monomers = 200)
parameter (total_particles = chains * monomers)
integer*4 id(total_particles), ifx(total_particles), ify(total_particles), ifz(total_particles)
real*8 x(total_particles), y(total_particles), z(total_particles)
real*8 x_buf(total_particles), y_buf(total_particles), z_buf(total_particles)

integer*4 i, j
integer*4 ix, iy, iz
integer*4 ux, uy, uz, sh
integer*4 p, p_next
integer*4 long_bonds
integer*4 xout, yout, zout
real*8 xcm, ycm, zcm
real*8 x_next, y_next, z_next
real*8 x_next_min, y_next_min, z_next_min
real*8 xi, yi, zi
real*8 xj, yj, zj
real*8 xij, yij, zij
real*8 dsq, dsq_min

integer*8 timestep
integer*4 np
real*8 xmin, ymin, zmin
real*8 xmax, ymax, zmax
real*8 xy, xz, yz
real*8 rho, drho
real*8 L, L_half

integer*8 f, istr, iend, incr
character(10) chrf10
character(80) chrf80
character(80) prefix
character(80) prefix_out
character(80) flnm
character(80) flnm_out

logical*4 checks
logical*4 flags

sh = 5
checks = .false.
flags = .false.

prefix = '../../sllod_M200N400_rings_G0.00001.'
prefix_out = 'sllod_M200N400_rings_G0.00001.'
istr = 40000000
iend = istr
incr = 100000

long_bonds = 0
do f = istr, iend, incr

  write(chrf10, '(i10)') f
  flnm = trim(adjustl(prefix)) // trim(adjustl(chrf10))
  flnm_out = trim(adjustl(prefix_out)) // trim(adjustl(chrf10))
  write(*,*) 'Working on file ', flnm

  open(10, file = flnm)
  read(10,*) chrf80
  read(10,*) timestep
  read(10,*) chrf80
  read(10,*) np
  read(10,*) chrf80
  read(10,*) xmin, xmax, xy
  read(10,*) ymin, ymax, xz
  read(10,*) zmin, zmax, yz
  read(10,*) chrf80
  if(flags) then
    do i = 1, total_particles
      read(10,*) id(i), x_buf(i), y_buf(i), z_buf(i), ifx(i), ify(i), ifz(i)
    enddo
  else
    ifx = 0
    ify = 0
    ifz = 0
    do i = 1, total_particles
      read(10,*) id(i), x_buf(i), y_buf(i), z_buf(i)
    enddo
  endif
  close(10)

  ! consistency checks
  write(*,*) 'r(1) = ', x_buf(1), y_buf(1), z_buf(1), ifx(1), ify(1), ifz(1)
  write(*,*) 'r(N) = ', x_buf(np), y_buf(np), z_buf(np), ifx(np), ify(np), ifz(np)
  if(timestep.ne.f) write(*,*) 'WARNING: Mismatch between timestep and file name: ', timestep, f
  if(np.ne.total_particles) write(*,*) 'WARNING: Mismatch between np and total_particles.'
  if(ymax - ymin.ne.zmax - zmin) write(*,*) 'WARNING: box size'
  rho = total_particles / ((xmax - xmin - abs(xz)) * (ymax - ymin) * (zmax - zmin))
  drho = abs(1.0D0 - rho / 0.85D0)
  if(drho.gt.0.001D0) write(*,*) 'WARNING: Density is not 0.85: ', rho
  write(*,*) 'tilt factors: ', xy, xz, yz

  L = ymax - ymin
  L_half = 0.5D0 * L
  write(*,*) 'L =  ', L

  ! sort by id (ignore flags at moment)
  do i = 1, total_particles
    j = id(i)
    x(j) = x_buf(i) + L * ifx(i) * 0
    y(j) = y_buf(i) + L * ify(i) * 0
    z(j) = z_buf(i) + L * ifz(i) * 0
  enddo

  ! wrap particles to central orthorhombic cell (xyz order matters)
  do i = 1, total_particles
    iy = int(y(i) / L)
    if(y(i).lt.0.0D0) iy = iy - 1
    y(i) = y(i) - iy * L

    iz = int(z(i) / L)
    if(z(i).lt.0.0D0) iz = iz - 1
    z(i) = z(i) - iz * L
    x(i) = x(i) - iz * xz

    ix = int(x(i) / L)
    if(x(i).lt.0.0D0) ix = ix - 1
    x(i) = x(i) - ix * L
  enddo

  ! check that particles are inside the central cell
  xout = 0
  yout = 0
  zout = 0
  do i = 1, total_particles
    !if(x(i).gt.L + z(i) * xz / L) xout = xout + 1
    if(x(i).gt.L) xout = xout + 1
    if(y(i).gt.L) yout = yout + 1
    if(z(i).gt.L) zout = zout + 1
    !if(x(i).lt.z(i) * xz / L) xout = xout + 1
    if(x(i).lt.0.0D0) xout = xout + 1
    if(y(i).lt.0.0D0) yout = yout + 1
    if(z(i).lt.0.0D0) zout = zout + 1
  enddo
  if(xout + yout + zout.ne.0) write(*,*) 'xout, yout, zout = ', xout, yout, zout

  x_buf = -1e10
  y_buf = -1e10
  z_buf = -1e10

  do i = 1, chains

    p = (i - 1) * monomers + 1
    x_buf(p) = x(p)
    y_buf(p) = y(p)
    z_buf(p) = z(p)

    do j = 2, monomers

      dsq_min = 1e10
      p_next = (i - 1) * monomers + j
      do ux = -sh, sh, 1
        do uy = -sh, sh, 1
          do uz = -sh, sh, 1
            x_next = ux * L + x(p_next) + uz * xz
            y_next = uy * L + y(p_next)
            z_next = uz * L + z(p_next)
            dsq = (x_buf(p_next - 1) - x_next)**2 + &
                  (y_buf(p_next - 1) - y_next)**2 + &
                  (z_buf(p_next - 1) - z_next)**2
            if(dsq < dsq_min) then
              x_next_min = x_next
              y_next_min = y_next
              z_next_min = z_next
              dsq_min = dsq
            endif
          enddo
        enddo
      enddo

      x_buf(p_next) = x_next_min
      y_buf(p_next) = y_next_min
      z_buf(p_next) = z_next_min
      if(dsq_min.gt.1.5D0.or.dsq_min.lt.0.64D0) then
        write(*,*) 'ERROR: long bond = ', dsq_min**0.5D0, ' May need to increase sh.'
        long_bonds = long_bonds + 1
      endif
    enddo

  enddo
  if(long_bonds.ne.0) write(*,*) "ERROR: long bonds = ", long_bonds

  ! wrap by COM (order is important)
  do i = 1, chains
    xcm = 0.0D0
    ycm = 0.0D0
    zcm = 0.0D0
    do j = 1, monomers
      p = (i - 1) * monomers + j
      xcm = xcm + x_buf(p)
      ycm = ycm + y_buf(p)
      zcm = zcm + z_buf(p)
    enddo
    xcm = xcm / monomers
    ycm = ycm / monomers
    zcm = zcm / monomers
    ! wrap in y
    iy = int(ycm / L)
    if(ycm.lt.0.0D0) iy = iy - 1
    do j = 1, monomers
      p = (i - 1) * monomers + j
      y_buf(p) = y_buf(p) - iy * L
    enddo
    ! wrap in z (implies changing x)
    iz = int(zcm / L)
    if(zcm.lt.0.0D0) iz = iz - 1
    do j = 1, monomers
      p = (i - 1) * monomers + j
      x_buf(p) = x_buf(p) - iz * xz
      z_buf(p) = z_buf(p) - iz * L
    enddo
    ! wrap in x
    xcm = xcm - xz * iz
    ix = int(xcm / L)
    if(xcm.lt.0.0D0) ix = ix - 1
    do j = 1, monomers
      p = (i - 1) * monomers + j
      x_buf(p) = x_buf(p) - ix * L
    enddo
  enddo

  if(checks) then
    ! check the min distance between all pairs
    dsq_min = 1e10
    do i = 1, total_particles - 1
      xi = x_buf(i)
      yi = y_buf(i)
      zi = z_buf(i)
      do j = i + 1, total_particles
        xij = xi - x_buf(j)
        yij = yi - y_buf(j)
        zij = zi - z_buf(j)
        dsq = xij**2 + yij**2 + zij**2
        if(dsq.lt.dsq_min) dsq_min = dsq
      enddo
    enddo
    write(*,*) 'min rij = ', dsq_min**0.5D0

    ! check the min distance between all pairs with minimum image convention
    ! only works for a cubic orthorhombic box, must be modified for triclinic
    dsq_min = 1e10
    do i = 1, total_particles - 1

      xi = x_buf(i)
      yi = y_buf(i)
      zi = z_buf(i)

      ix = int(xi / L)
      iy = int(yi / L)
      iz = int(zi / L)

      if(xi.lt.0.0D0) ix = ix - 1
      if(yi.lt.0.0D0) iy = iy - 1
      if(zi.lt.0.0D0) iz = iz - 1

      xi = xi - ix * L
      yi = yi - iy * L
      zi = zi - iz * L

      do j = i + 1, total_particles

        xj = x_buf(j)
        yj = y_buf(j)
        zj = z_buf(j)

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

        dsq = xij**2 + yij**2 + zij**2
        if(dsq.lt.dsq_min) dsq_min = dsq
      enddo
    enddo
    write(*,*) 'min rij = ', dsq_min**0.5D0
  endif

  ! write file with reconstructed chains
  flnm_out = 'recon_'//flnm_out
  open(10, file = flnm_out)
  write(10,*) 'ITEM: TIMESTEP'
  write(10,'(i12)') timestep
  write(10,*) 'ITEM: NUMBER OF ATOMS'
  write(10,'(i7)') np
  write(10,*) 'ITEM: BOX BOUNDS xy xz yz'
  write(10,'(f10.4,f10.4,f12.6)') 0.0000, L, xy
  write(10,'(f10.4,f10.4,f12.6)') 0.0000, L, xz
  write(10,'(f10.4,f10.4,f12.6)') 0.0000, L, yz
  write(10,*) 'ITEM: ATOMS id xcon ycon zcon'
  do i = 1, total_particles
    write(10,'(i7,f10.3,f10.3,f10.3)') i, x_buf(i), y_buf(i), z_buf(i)
  enddo
  close(10)

enddo

end program
