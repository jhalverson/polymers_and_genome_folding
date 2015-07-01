program rdfex

! This program computes an excluded monomer-monomer RDF of a
! polymer melt where the monomers of the central chain are
! ignored in the analysis. Distances are computed using the
! center-of-mass of the central chain and the monomers of the
! other chains in the system.

! This program requires three inputs:
! * number of monomers
! * file names of the configurations 
! * output filename

implicit none

integer*4 chains, monomers, np, np27, chains27
integer*4 ibinmax
real*8 L, L2sq, dr, rho
parameter (chains = 200)
parameter (monomers = 1600)
parameter (np = chains * monomers)
parameter (np27 = 27 * np)
parameter (chains27 = 27 * chains)
parameter (rho = 0.85)
parameter (L = (np / rho)**(1.0 / 3.0))
parameter (L2sq = (L / 2.0)**2)
parameter (dr = 0.25)
parameter (ibinmax = 0.5 * L / dr + 1)
integer*4 ihist(ibinmax)
real*8 rhist(ibinmax)
real*8 x(np27), y(np27), z(np27)

integer*4 i, j, k, m, f
integer*4 ix, iy, iz
integer*4 ipart, ict, ibin

character(9) chrf
character(13) flnm
integer*4 istr, iend, incr, files

real*8 rsq, rij
real*8 xcm, ycm, zcm
real*8 pi, fac, rinner, router

real*8 xmin, ymin, zmin
real*8 xmax, ymax, zmax
real*8 xb, yb, zb

logical*4 debug
debug = .false.

! initialize histogram
ihist = 0

istr = 360800000
iend = 362000000
incr = 100000
files = 0

! need Python script (block_strip.py) to prepare blockfiles

do f = istr, iend, incr

  write(chrf,'(i9)') f
  flnm = chrf(1:9) // '.pos'

  write(*,*) "File ", f, " of ", (iend - istr) / incr
  
  ! read in coordinates
  open(10, file = flnm)
  do i = 1, np
    read(10,*) x(i), y(i), z(i)
  enddo
  close(10)
  
  if(debug) then
  
    ! find smallest bounding box
    xmin = 1e10; ymin = xmin; zmin = xmin
    xmax = -1e10; ymax = xmax; zmax = xmax
    do i = 1, np
  
      if(x(i).lt.xmin) xmin = x(i)
      if(x(i).gt.xmax) xmax = x(i)
  
      if(y(i).lt.ymin) ymin = y(i)
      if(y(i).gt.ymax) ymax = y(i)
  
      if(z(i).lt.zmin) zmin = z(i)
      if(z(i).gt.zmax) zmax = z(i)
  
    enddo
  
    xb = xmax - xmin
    yb = ymax - ymin
    zb = zmax - zmin
  
    write(*,*) "Particles fit in a box of dimensions: "
    write(*,'(f8.3,f8.3,f8.3)') xb, yb, zb
    write(*,*) "Simulation cell is of dimensions    : "
    write(*,'(f8.3,f8.3,f8.3)') L, L, L
  
  endif
  
  ! wrap chains to central simulation box
  do i = 1, chains
  
    ! compute center-of-mass of molecule i
    xcm = 0.0
    ycm = 0.0
    zcm = 0.0
    do j = 1, monomers
      ipart = (i - 1) * monomers + j
      xcm = xcm + x(ipart)
      ycm = ycm + y(ipart)
      zcm = zcm + z(ipart)
    enddo
    xcm = xcm / monomers
    ycm = ycm / monomers
    zcm = zcm / monomers
  
    ! shift chains such that the central box is between [0, L]
    ix = xcm / L
    iy = ycm / L
    iz = zcm / L
  
    ! correct in case of negative positions
    if(xcm.lt.0) ix = ix - 1
    if(ycm.lt.0) iy = iy - 1
    if(zcm.lt.0) iz = iz - 1
  
    ! perform shift
    do j = 1, monomers
      ipart = (i - 1) * monomers + j
      x(ipart) = x(ipart) - ix * L 
      y(ipart) = y(ipart) - iy * L
      z(ipart) = z(ipart) - iz * L
    enddo
  
  enddo
  
  ! create 26 translations
  ict = 0
  ipart = np
  do i = -1, 1, 1
    do j = -1, 1, 1
      do k = -1, 1, 1
        if(i**2 + j**2 + k**2.ne.0) then
          ict = ict + 1
          do m = 1, np
            ipart = ipart + 1
            x(ipart) = x(ipart - ict * np) + i * L
            y(ipart) = y(ipart - ict * np) + j * L
            z(ipart) = z(ipart - ict * np) + k * L
          enddo
        endif
      enddo
    enddo
  enddo
  
  if(debug) then
  
    write(*,*) "L = ", L, " L2sq = ", L2sq, " dr = ", dr
    write(*,*) "np = ", np
    write(*,*) " ipart = ", ipart, " np27 = ", np27
    write(*,*) "ict = ", ict, " ibinmax = ", ibinmax
  
  endif
  
  ! loop over each central molecule
  do i = 1, chains
  
    ! compute center-of-mass of molecule i
    xcm = 0.0
    ycm = 0.0
    zcm = 0.0
    do m = 1, monomers
      ipart = (i - 1) * monomers + m
      xcm = xcm + x(ipart)
      ycm = ycm + y(ipart)
      zcm = zcm + z(ipart)
    enddo
    xcm = xcm / monomers
    ycm = ycm / monomers
    zcm = zcm / monomers
  
    ! loop over all monomers except for those of molecule i
    do j = 1, chains27
      do k = 1, monomers
        if(i.ne.j) then
          ipart = (j - 1) * monomers + k
          rsq = (x(ipart) - xcm)**2 + (y(ipart) - ycm)**2 + (z(ipart) - zcm)**2
          if(rsq.le.L2sq) then
            rij = sqrt(rsq)
            ibin = int(rij / dr) + 1
            if(ibin.gt.ibinmax) write(*,*) "ERROR: ibin"
            ihist(ibin) = ihist(ibin) + 1
          endif
        endif
      enddo
    enddo
  
  enddo
  
  files = files + 1
  
enddo
  
! normalize
pi = 4.0 * atan(1.0)
fac = 4.0 * pi * rho / 3.0
do i = 1, ibinmax
  rinner = (i - 1) * dr
  router = rinner + dr
  rhist(i) = ihist(i) / (fac * (router**3 - rinner**3)) / (chains * files)
enddo
  
! write to file
open(11, file="1600.out")
do i = 1, ibinmax - 1
  write(11,'(f8.3,f8.3,f8.3,f8.3)') (i - 1) * dr, i * dr, i * dr - 0.5 * dr, rhist(i)
enddo
close(11)
  
end program
