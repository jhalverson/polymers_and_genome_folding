program neighbors_second_kind

implicit none

integer*4 chains, monomers, total_particles
real*8 L, L_half
parameter (chains = 200)
parameter (monomers = 100)
parameter (total_particles = chains * monomers)
parameter (L = 28.6552162316D0)
!parameter (L = 36.1033101196D0)
!parameter (L = 45.4873203905D0)
!parameter (L = 57.3104324633D0)
!parameter (L = 72.2066202391D0)
real*8 x(total_particles), y(total_particles), z(total_particles)
real*8 xcm, ycm, zcm

integer*4 i, j, k, p, files, k1, k2
integer*4 ix, iy, iz

character(5) chrf5
character(6) chrf6
character(7) chrf7
character(8) chrf8
character(9) chrf9
character(10) chrf10
character(50) flnm
integer*4 f, istr, iend, incr

real*8 rijsq, rcsq
real*8 xi, yi, zi
real*8 xj, yj, zj
real*8 xij, yij, zij

! Rg2
!rcsq = 17.16D0
!rcsq = 30.77D0
!rcsq = 52.92D0
!rcsq = 87.61D0
!rcsq = 145.61D0

! Re2
rcsq = 50.76D0
!rcsq = 88.75D0
!rcsq = 149.41D0
!rcsq = 242.41D0
!rcsq = 401.73D0

L_half = 0.5D0 * L
k1 = 0
k2 = 0
files = 0

istr = 500000
iend = 2000000000
incr = 5000000

do f = istr, iend, incr

  if(f.lt.100000) then
     write(chrf5,'(i5)') f
     flnm = chrf5 // '.pos'
  elseif(f.lt.1000000) then
     write(chrf6,'(i6)') f
     flnm = chrf6 // '.pos'
  elseif(f.lt.10000000) then
     write(chrf7,'(i7)') f
     flnm = chrf7 // '.pos'
  elseif(f.lt.100000000) then
     write(chrf8,'(i8)') f
     flnm = chrf8 // '.pos'
  elseif(f.lt.1000000000) then
     write(chrf9,'(i9)') f
     flnm = chrf9 // '.pos'
  else
     write(chrf10,'(i10)') f
     flnm = chrf10 // '.pos'
  endif
  flnm = '../sqt/stripped/' // flnm
  
  ! read in coordinates
  open(10, file = flnm)
  do i = 1, total_particles
    read(10,*) x(i), y(i), z(i)
  enddo
  close(10)
  
  if(f.eq.istr) then
 
    write(*,*) "computed rho = ", total_particles / L**3
    write(*,*) "L = ", L, ", L_half = ", L_half
    write(*,*) "total_particles = ", total_particles, "  monomers = ", monomers
    write(*,*) "r(1) = ", x(1), y(1), z(1)
    write(*,*) "r(end) = ", x(total_particles), y(total_particles), z(total_particles)
 
  endif

  do i = 1, chains

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

    do j = 1, chains

      if(i.ne.j) then

        do k = 1, monomers

          p = (j - 1) * monomers + k
          xj = x(p)
          yj = y(p)
          zj = z(p)

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
            k2 = k2 + 1
            exit
          endif

        enddo

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
        if(rijsq.le.rcsq) then
          k1 = k1 + 1
        endif

      endif

    enddo

  enddo
 
  files = files + 1
  write(*,*) real(k1) / (files * chains), real(k2) / (files * chains)
 
enddo

end program
