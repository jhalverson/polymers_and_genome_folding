program inter_per_chain

implicit none

integer*8 chains, monomers, total_particles
parameter (chains = 200)
parameter (monomers = 200)
parameter (total_particles = chains * monomers)
real*8 x(total_particles), y(total_particles), z(total_particles)

integer*4 i, j
integer*4 ichain, jchain
integer*4 ix, iy, iz
integer*8 ct, intra, inter, sum_inter
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
real*8 U

! cutoff
rc = 2.0D0**(1.0D0 / 6.0D0)
rcsq = rc**2

L = (total_particles / 0.85D0)**(1.0 / 3.0)
L_half = 0.5D0 * L

prefix = '../../200/'
istr = 1000000000
iend = 2000000000
incr = 100000000

sum_inter = 0
files = 0

! need Python script (block_strip.py) to prepare blockfiles
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
    write(*,*) "L = ", L, ", L_half = ", L_half, ", rcsq = ", rcsq
    write(*,*) "total_particles = ", total_particles, "  monomers = ", monomers
    write(*,*) "r(1) = ", x(1), y(1), z(1)
    write(*,*) "r(end) = ", x(total_particles), y(total_particles), z(total_particles)
 
  endif

  U = 0.0D0
  ct = 0
  intra = 0
  inter = 0
  do i = 1, total_particles - 1

    ichain = (i - 1) / monomers

    xi = x(i)
    yi = y(i)
    zi = z(i)

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

      ct = ct + 1

      jchain = (j - 1) / monomers

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
        U = U + 4.0D0 * (1.0D0 / rijsq**6 - 1.0D0 / rijsq**3)
        if(ichain.eq.jchain) then
          intra = intra + 2
        else
          inter = inter + 2
        endif
      endif

    enddo

  enddo

  sum_inter = sum_inter + inter
  files = files + 1
 
  write(*,*) "avg number of interchain interactions per chain = ", real(inter) / chains
  write(*,*) "U = ", U
  write(*,*) "counts = ", ct, total_particles * (total_particles - 1) / 2
 
enddo

write(*,*) "final avg number of interchain interactions per chain = ", real(sum_inter) / (chains * files)
 
end program
